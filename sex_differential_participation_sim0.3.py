#!/usr/bin/env python
# coding: utf-8

import hail as hl
import numpy as np
from math import exp, log

hl.init()


def gwas(y, x, cov):
    g = hl.linear_regression_rows(y=y,
                                  x=x,
                                  covariates=cov,
                                  pass_through=['rsid'])
    return g


def export_gwas(g, fname):
    gann = g.annotate(A1=g.alleles[0],
                      A2=g.alleles[1]).key_by()
    gann.select('locus', 'A1', 'A2', 'rsid', 'n', 'beta', 'standard_error', 'p_value').export(fname)
    return


mt = hl.read_matrix_table('gs://mattia-simulations/simEUR350_2.mt')
out_bucket = 'gs://mattia-simulations/simulations_k_or/'

# Save mt cols to table and convert to pandas dataframe for sampling
# mt.y:
# y0:y3    h2 = 0.1
# y4:y7    h2 = 0.3

mt = mt.annotate_cols(y0=mt.y[4])
mt = mt.annotate_cols(y1=mt.y[5])
df = mt.cols().select('s', 'y0', 'y1', 'sex').key_by().to_pandas()


# OR = [1, 1.2, 1.5, 1.8, 2, 3]
OR = [1]
K = [-0.5, -0.3, 0, 0.3, 0.7, 1, 1.5]

for k in K:
    for o in OR:

        df['z'] = df['y0'] * log(o) + df['y1'] * log(o) + np.random.normal(0, 0.1, len(df.index))
        df.loc[(df['sex'] == 0), 'z'] = df.loc[(df['sex'] == 0), 'z'] * k
        df['prob'] = 0.5 if o == 1 else [1 / (1 + exp(-z)) for z in df['z']]
        df['select'] = np.random.binomial(n=1, p=df['prob'], size=len(df.index))

        samples_to_keep = set(df.loc[(df['select'] == 1), 's'])
        set_to_keep = hl.literal(samples_to_keep)
        mt_sampled = mt.filter_cols(set_to_keep.contains(mt['s']), keep=True)

        i = '_'+str(k)+'_'+str(o)

        # Export phenotypes
        mt_sampled.cols().select('s', 'sex', 'y0', 'y1').key_by().\
            export(out_bucket + 'phenotypes/pheno' + i + '.tsv')

        # ------------ sex -------------
        # GWAS of sex
        gwas_s = gwas(mt_sampled.sex, mt_sampled.GT.n_alt_alleles(), [1.0])
        fn = out_bucket + 'gwas/gwas_sex' + i + '.tsv'
        export_gwas(gwas_s, fn)

        # ------------ y0 --------------
        # GWAS of y0
        gwas_y0 = gwas(mt_sampled.y0, mt_sampled.GT.n_alt_alleles(), [1.0])
        fn = out_bucket + 'gwas/gwas_y0' + i + '.tsv'
        export_gwas(gwas_y0, fn)

        # GWAS of y0 adjusted for sex
        gwas_y0_adj = gwas(mt_sampled.y0, mt_sampled.GT.n_alt_alleles(), [1.0, mt_sampled.sex])
        fn = out_bucket + 'gwas/gwas_y0' + i + '_adj.tsv'
        export_gwas(gwas_y0_adj, fn)

        # GWAS of y0 for females
        mt_sampled_f = mt_sampled.filter_cols(mt_sampled.sex == 1)
        gwas_y0_f = gwas(mt_sampled_f.y0, mt_sampled_f.GT.n_alt_alleles(), [1.0])
        fn = out_bucket + 'gwas/gwas_y0' + i + '_f.tsv'
        export_gwas(gwas_y0_f, fn)

        # GWAS of y0 for males
        mt_sampled_m = mt_sampled.filter_cols(mt_sampled.sex == 0)
        gwas_y0_m = gwas(mt_sampled_m.y0, mt_sampled_m.GT.n_alt_alleles(), [1.0])
        fn = out_bucket + 'gwas/gwas_y0' + i + '_m.tsv'
        export_gwas(gwas_y0_m, fn)

        # ------------ y1 ---------------
        # GWAS of y1
        gwas_y1 = gwas(mt_sampled.y1, mt_sampled.GT.n_alt_alleles(), [1.0])
        fn = out_bucket + 'gwas/gwas_y1' + i + '.tsv'
        export_gwas(gwas_y1, fn)

        # GWAS of y1 adjusted for sex
        gwas_y1_adj = gwas(mt_sampled.y1, mt_sampled.GT.n_alt_alleles(), [1.0, mt_sampled.sex])
        fn = out_bucket + 'gwas/gwas_y1' + i + '_adj.tsv'
        export_gwas(gwas_y1_adj, fn)

        # GWAS of y1 for females
        mt_sampled_f = mt_sampled.filter_cols(mt_sampled.sex == 1)
        gwas_y1_f = gwas(mt_sampled_f.y1, mt_sampled_f.GT.n_alt_alleles(), [1.0])
        fn = out_bucket + 'gwas/gwas_y1' + i + '_f.tsv'
        export_gwas(gwas_y1_f, fn)

        # GWAS of y1 for males
        mt_sampled_m = mt_sampled.filter_cols(mt_sampled.sex == 0)
        gwas_y1_m = gwas(mt_sampled_m.y1, mt_sampled_m.GT.n_alt_alleles(), [1.0])
        fn = out_bucket + 'gwas/gwas_y1' + i + '_m.tsv'
        export_gwas(gwas_y1_m, fn)
