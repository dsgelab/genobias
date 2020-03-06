#!/usr/bin/env python
# coding: utf-8
import hail as hl
import random
hl.init()

# # # REQUIRES HAIL v0.2.24


# # h2 and rgs specifications

# Simulate genetically unrelated phenotypes with h2 = 0.1, 0.3

# h2 = [0.1, 0.1, 0.1, 0.1, 0.3, 0.3, 0.3, 0.3]
# rg = [0]*28
# out = 'gs://mattia/mattia-simulations/simEUR350_2.mt'

# Simulate genetically correlated phenotypes with same h2 = 0.3 and different rgs + one uncorrelated

h2 = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3]
rg = [-0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5, 0] + 28*[0]
out = 'gs://mattia/mattia-simulations/simEUR350_correlated.mt'

# # Simulations
# Import UKBB data
mt = hl.import_plink(bed='gs://mattia/simulations/ukb_imp_chr1_v3.bed',
                     bim='gs://mattia/simulations/ukb_imp_chr1_v3.bim',
                     fam='gs://mattia/simulations/ukb_imp_chr1_v3.fam',
                     reference_genome='GRCh37')

for chrom in range(2, 23):
    mtT = hl.import_plink(bed='gs://mattia/simulations/ukb_imp_chr%s_v3.bed' % chrom,
                          bim='gs://mattia/simulations/ukb_imp_chr%s_v3.bim' % chrom,
                          fam='gs://mattia/simulations/ukb_imp_chr%s_v3.fam' % chrom,
                          reference_genome='GRCh37')
    mt = mt.union_rows(mtT)

# Import PCA
tb = hl.import_table('gs://mattia-simulations/ukbb_pca_pops_rf.txt', impute=True)
tb = tb.annotate(s_str=hl.str(tb.s)).key_by('s_str')
mt = mt.annotate_cols(sample_info=tb[mt.s])

# Keep only unrelated (~361k samples)
tb2 = hl.import_table('gs://mattia/mattia-simulations/ukb31063.neale_gwas_samples.both_sexes.txt', impute=True)
tb2 = tb2.annotate(s_str=hl.str(tb2.s)).key_by('s_str')
mt = mt.semi_join_cols(tb2).add_col_index()
mt = mt.annotate_cols(s_index=mt.col_idx).key_cols_by('s_index')

# Extract 350k samples + basic QC
random.seed(123)
indices = random.sample(range(mt.count_cols()), 350000)
mt = mt.choose_cols(list(indices))
mt = hl.variant_qc(mt)
mt = mt.filter_rows((mt.variant_qc.AF[1] >= 0.05) & (mt.variant_qc.AF[1] <= 0.95))

# Randomly assign sex (1 F, 0 M)
mt = mt.annotate_cols(sex=hl.cond(hl.rand_bool(0.5, seed=123), 1, 0))

# Simulate phenotypes and save MatrixTable
sim = hl.experimental.ldscsim.simulate_phenotypes(mt, mt.GT, h2=h2, rg=rg)
print(sim.describe())
sim.write(out, overwrite=True)
