#!/usr/bin/env python
# coding: utf-8
import hail as hl
import random
hl.init()

# # # REQUIRES HAIL v0.2.24

# # # Suggested gcloud cluster configuration:
# hailctl dataproc start hail-0.2.24 --num-preemptible-workers=36 --worker-machine-type=n1-highmem-8


# # h2 and rgs specifications
# Simulate genetically unrelated phenotypes with h2 = 0.1, 0.3

# h2 = [0.1, 0.1, 0.1, 0.1, 0.3, 0.3, 0.3, 0.3]
# rg = [0]*28
# out = 'gs://.../sim_350k_uncorrelated.mt'

# Simulate genetically correlated phenotypes with same h2 = 0.3 and different rgs + one uncorrelated

h2 = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3]
rg = [-0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5, 0] + 28*[0]
out = 'gs://.../sim_350k_correlated.mt'

# # Simulations
# Import UKBB hm3 genotype data
mt = hl.import_plink(bed='gs://.../ukb_imp_chr1_v3.bed',
                     bim='gs://.../ukb_imp_chr1_v3.bim',
                     fam='gs://.../ukb_imp_chr1_v3.fam',
                     reference_genome='GRCh37')

for chrom in range(2, 23):
    mtT = hl.import_plink(bed='gs://.../ukb_imp_chr%s_v3.bed' % chrom,
                          bim='gs://.../ukb_imp_chr%s_v3.bim' % chrom,
                          fam='gs://.../ukb_imp_chr%s_v3.fam' % chrom,
                          reference_genome='GRCh37')
    mt = mt.union_rows(mtT)

# Keep only unrelated (~361k samples)
tb2 = hl.import_table('gs://.../unrelated_samples.txt', impute=True)
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
