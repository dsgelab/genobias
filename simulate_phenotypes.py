#!/usr/bin/env python
# coding: utf-8
import hail as hl
import hail.expr.aggregators as agg
import hail.methods
import pandas as pd
from typing import *
import random
hl.init()

import requests
url = 'https://raw.githubusercontent.com/nikbaya/ldscsim/master/ldscsim.py'
r = requests.get(url).text
exec(r)

mt = hl.import_plink(bed='gs://risteys-data-transfer/simulations/ukb_imp_chr1_v3.bed',
    bim='gs://risteys-data-transfer/simulations/ukb_imp_chr1_v3.bim',
    fam='gs://risteys-data-transfer/simulations/ukb_imp_chr1_v3.fam',
    reference_genome='GRCh37')

for chrom in range(2, 23):
    mtT = hl.import_plink(bed='gs://risteys-data-transfer/simulations/ukb_imp_chr%s_v3.bed' % chrom,
                          bim='gs://risteys-data-transfer/simulations/ukb_imp_chr%s_v3.bim' % chrom,
                          fam='gs://risteys-data-transfer/simulations/ukb_imp_chr%s_v3.fam' % chrom,
                          reference_genome='GRCh37')
    
    mt = mt.union_rows(mtT)

tb = hl.import_table('gs://risteys-data-transfer/simulations/ukbb_pca_pops_rf.txt', impute=True)
tb = tb.annotate(s_str=hl.str(tb.s)).key_by('s_str')
mt = mt.annotate_cols(sample_info=tb[mt.s])

# Keep only unrelated individuals
tb2 = hl.import_table('gs://mattia-simulations/ukb31063.neale_gwas_samples.both_sexes.txt', impute=True)
tb2 = tb2.annotate(s_str=hl.str(tb2.s)).key_by('s_str')
mt = mt.semi_join_cols(tb2)

mt_eur = mt.filter_cols((mt.sample_info.pop == "EUR") & (hl.is_defined(mt.sample_info.pop)))
mt_eur = mt_eur.checkpoint('gs://mattia-simulations/mt_eur_checkpoint.mt', overwrite=True)

mt_eur = hl.variant_qc(mt_eur)
mt_eur = mt_eur.filter_rows((mt_eur.variant_qc.AF[1] >= 0.05) & (mt_eur.variant_qc.AF[1] <= 0.95)).add_col_index()
mt_eur = mt_eur.annotate_cols(s_index=mt_eur.col_idx).key_cols_by('s_index')

# Randomly assign sex (1 F, 0 M)
random.seed(123)
mt_eur = mt_eur.annotate_cols(sex=hl.cond(hl.rand_bool(0.5), 1, 0))
mt_eur = mt_eur.checkpoint('gs://mattia-simulations/mt_eur_qc_checkpoint.mt', overwrite=True)

sim1 = simulate_phenotypes(mt_eur, mt_eur.GT.n_alt_alleles(), h2=[0.476, 0.462], rg=0)
sim1.write('gs://mattia-simulations/simEUR350_mcv_height.mt', overwrite=True)
