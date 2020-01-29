#!/usr/bin/env python
# coding: utf-8
import hail as hl
import random
hl.init(log="/home/bordinki/hail.log")

# import requests
# url = 'https://raw.githubusercontent.com/nikbaya/ldscsim/master/ldscsim.py'
# r = requests.get(url).text
# exec(r)
# !!! Nik added a "exact_h2" flag to the calculate_phenotypes method, still working on so the function from his repo is
# not working

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

# # Keep only ukb31063 (~361xxx samples)
tb2 = hl.import_table('gs://mattia-simulations/ukb31063.neale_gwas_samples.both_sexes.txt', impute=True)
tb2 = tb2.annotate(s_str=hl.str(tb2.s)).key_by('s_str')
mt = mt.semi_join_cols(tb2).add_col_index()
mt = mt.annotate_cols(s_index=mt.col_idx).key_cols_by('s_index')

# Extract 350000 samples
random.seed(123)
indices = random.sample(range(mt.count_cols()), 350000)
mt = mt.choose_cols(list(indices))
# QC
mt = hl.variant_qc(mt)
mt = mt.filter_rows((mt.variant_qc.AF[1] >= 0.05) & (mt.variant_qc.AF[1] <= 0.95))
# Randomly assign sex (1 F, 0 M)
mt = mt.annotate_cols(sex=hl.cond(hl.rand_bool(0.5, seed=123), 1, 0))

# Simulate phenotypes
# sim = simulate_phenotypes(mt, mt.GT, h2=[0.1, 0.1, 0.1, 0.1, 0.3, 0.3, 0.3, 0.3], rg=[0]*28)
# ! Nik added a "exact_h2" flag to the calculate_phenotypes method, still working on so the function from his repo is
# not working

sim = hl.experimental.ldscsim.simulate_phenotypes(mt, mt.GT, h2=[0.1, 0.1, 0.1, 0.1, 0.3, 0.3, 0.3, 0.3], rg=[0]*28)
print(sim.describe())
sim.write('gs://mattia-simulations/simEUR350_2.mt', overwrite=True)
