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


mt = hl.read_matrix_table('gs://mattia-simulations/simEUR350.mt')

output_bucket = 'gs://mattia-simulations/EUR350_scenario_1/'

# Subsampling with different prob for M and F and running GWAS
# Save mt cols to table and convert to pandas df for sampling
mt = mt.annotate_cols(y0=mt.y[0])
mt = mt.annotate_cols(y1=mt.y[1])
mt = mt.annotate_cols(y2=mt.y[2])
mt = mt.annotate_cols(y3=mt.y[3])
df = mt.cols().select('s', 'y0', 'y1', 'y2', 'y3', 'sex').key_by().to_pandas()

# Export phenotypes
mt.cols().select('s',
                 'sex',
                 'y0',
                 'y1',
                 'y2',
                 'y3').key_by().export(
                 output_bucket + 'phenotypes/pheno_0.tsv')


OR_y0 = [1.2, 1.5, 1.8, 2, 4]
OR_y1 = [1.2, 1.5, 1.8, 2, 4]

output_bucket = 'gs://mattia-simulations/EUR350_scenario_1/10/'
i = 0
for ory0 in OR_y0:
    for ory1 in OR_y1:
        i += 1

        # Simulation scenario 1: M and F picked based on y0 and y1, inversely
        df['z'] = df['y0'] * log(ory0) + df['y1'] * log(ory1)
        df.loc[(df['sex'] == 0), 'z'] = -df['z']
        df['prob'] = [1 / (1 + exp(-z)) for z in df['z']]
        df['select'] = np.random.binomial(n=1, p=df['prob'], size=len(df.index))

        samples_to_keep = set(df.loc[(df['select'] == 1), 's'])
        set_to_keep = hl.literal(samples_to_keep)
        mt_sampled = mt.filter_cols(set_to_keep.contains(mt['s']), keep=True)

        # Export phenotypes
        mt_sampled.cols().select('s',
                                 'sex',
                                 'y0',
                                 'y1',
                                 'y2',
                                 'y3').key_by().export(
                                 output_bucket + 'phenotypes/pheno_sample' + str(i) + '.tsv')

        # ------------ sex -------------
        # GWAS of sex
        gwas_s = gwas(mt_sampled.sex, mt_sampled.GT.n_alt_alleles(), [1.0])
        fn = output_bucket + 'gwas/gwas_sex_sample' + str(i) + '.tsv'
        export_gwas(gwas_s, fn)

        # ------------ y0 --------------
        # GWAS of y0
        gwas_y0 = gwas(mt_sampled.y[0], mt_sampled.GT.n_alt_alleles(), [1.0])
        fn = output_bucket + 'gwas/gwas_y0_sample' + str(i) + '.tsv'
        export_gwas(gwas_y0, fn)

        # GWAS of y0 adjusted for sex
        gwas_y0_adj = gwas(mt_sampled.y[0], mt_sampled.GT.n_alt_alleles(), [1.0, mt_sampled.sex])
        fn = output_bucket + 'gwas/gwas_y0_sample' + str(i) + '_adj.tsv'
        export_gwas(gwas_y0_adj, fn)

        # GWAS of y0 for females
        mt_sampled_f = mt_sampled.filter_cols(mt_sampled.sex == 1)
        gwas_y0_f = gwas(mt_sampled_f.y[0], mt_sampled_f.GT.n_alt_alleles(), [1.0])
        fn = output_bucket + 'gwas/gwas_y0_f_sample' + str(i) + '.tsv'
        export_gwas(gwas_y0_f, fn)

        # GWAS of y0 for males
        mt_sampled_m = mt_sampled.filter_cols(mt_sampled.sex == 0)
        gwas_y0_m = gwas(mt_sampled_m.y[0], mt_sampled_m.GT.n_alt_alleles(), [1.0])
        fn = output_bucket + 'gwas/gwas_y0_m_sample' + str(i) + '.tsv'
        export_gwas(gwas_y0_m, fn)

        # ------------ y1 ---------------
        # GWAS of y1
        gwas_y1 = gwas(mt_sampled.y[1], mt_sampled.GT.n_alt_alleles(), [1.0])
        fn = output_bucket + 'gwas/gwas_y1_sample' + str(i) + '.tsv'
        export_gwas(gwas_y1, fn)

        # GWAS of y1 adjusted for sex
        gwas_y1_adj = gwas(mt_sampled.y[1], mt_sampled.GT.n_alt_alleles(), [1.0, mt_sampled.sex])
        fn = output_bucket + 'gwas/gwas_y1_sample' + str(i) + '_adj.tsv'
        export_gwas(gwas_y1_adj, fn)

        # GWAS of y1 for females
        mt_sampled_f = mt_sampled.filter_cols(mt_sampled.sex == 1)
        gwas_y1_f = gwas(mt_sampled_f.y[1], mt_sampled_f.GT.n_alt_alleles(), [1.0])
        fn = output_bucket + 'gwas/gwas_y1_f_sample' + str(i) + '.tsv'
        export_gwas(gwas_y1_f, fn)

        # GWAS of y1 for males
        mt_sampled_m = mt_sampled.filter_cols(mt_sampled.sex == 0)
        gwas_y1_m = gwas(mt_sampled_m.y[1], mt_sampled_m.GT.n_alt_alleles(), [1.0])
        fn = output_bucket + 'gwas/gwas_y1_m_sample' + str(i) + '.tsv'
        export_gwas(gwas_y1_m, fn)


output_bucket = 'gs://mattia-simulations/EUR350_scenario_1/30/'
i = 0
for ory0 in OR_y0:
    for ory1 in OR_y1:
        i += 1

        # Simulation scenario 1: M and F picked based on y0 and y1, inversely
        df['z'] = df['y2'] * log(ory0) + df['y3'] * log(ory1)
        df.loc[(df['sex'] == 0), 'z'] = -df['z']
        df['prob'] = [1 / (1 + exp(-z)) for z in df['z']]
        df['select'] = np.random.binomial(n=1, p=df['prob'], size=len(df.index))

        samples_to_keep = set(df.loc[(df['select'] == 1), 's'])
        set_to_keep = hl.literal(samples_to_keep)
        mt_sampled = mt.filter_cols(set_to_keep.contains(mt['s']), keep=True)

        # Export phenotypes
        mt_sampled.cols().select('s',
                                 'sex',
                                 'y0',
                                 'y1',
                                 'y2',
                                 'y3').key_by().export(
            output_bucket + 'phenotypes/pheno_sample' + str(i) + '.tsv')

        # ------------ sex -------------
        # GWAS of sex
        gwas_s = gwas(mt_sampled.sex, mt_sampled.GT.n_alt_alleles(), [1.0])
        fn = output_bucket + 'gwas/gwas_sex_sample' + str(i) + '.tsv'
        export_gwas(gwas_s, fn)

        # ------------ y0 --------------
        # GWAS of y0
        gwas_y0 = gwas(mt_sampled.y[2], mt_sampled.GT.n_alt_alleles(), [1.0])
        fn = output_bucket + 'gwas/gwas_y0_sample' + str(i) + '.tsv'
        export_gwas(gwas_y0, fn)

        # GWAS of y0 adjusted for sex
        gwas_y0_adj = gwas(mt_sampled.y[2], mt_sampled.GT.n_alt_alleles(), [1.0, mt_sampled.sex])
        fn = output_bucket + 'gwas/gwas_y0_sample' + str(i) + '_adj.tsv'
        export_gwas(gwas_y0_adj, fn)

        # GWAS of y0 for females
        mt_sampled_f = mt_sampled.filter_cols(mt_sampled.sex == 1)
        gwas_y0_f = gwas(mt_sampled_f.y[2], mt_sampled_f.GT.n_alt_alleles(), [1.0])
        fn = output_bucket + 'gwas/gwas_y0_f_sample' + str(i) + '.tsv'
        export_gwas(gwas_y0_f, fn)

        # GWAS of y0 for males
        mt_sampled_m = mt_sampled.filter_cols(mt_sampled.sex == 0)
        gwas_y0_m = gwas(mt_sampled_m.y[2], mt_sampled_m.GT.n_alt_alleles(), [1.0])
        fn = output_bucket + 'gwas/gwas_y0_m_sample' + str(i) + '.tsv'
        export_gwas(gwas_y0_m, fn)

        # ------------ y1 ---------------
        # GWAS of y1
        gwas_y1 = gwas(mt_sampled.y[3], mt_sampled.GT.n_alt_alleles(), [1.0])
        fn = output_bucket + 'gwas/gwas_y1_sample' + str(i) + '.tsv'
        export_gwas(gwas_y1, fn)

        # GWAS of y1 adjusted for sex
        gwas_y1_adj = gwas(mt_sampled.y[3], mt_sampled.GT.n_alt_alleles(), [1.0, mt_sampled.sex])
        fn = output_bucket + 'gwas/gwas_y1_sample' + str(i) + '_adj.tsv'
        export_gwas(gwas_y1_adj, fn)

        # GWAS of y1 for females
        mt_sampled_f = mt_sampled.filter_cols(mt_sampled.sex == 1)
        gwas_y1_f = gwas(mt_sampled_f.y[3], mt_sampled_f.GT.n_alt_alleles(), [1.0])
        fn = output_bucket + 'gwas/gwas_y1_f_sample' + str(i) + '.tsv'
        export_gwas(gwas_y1_f, fn)

        # GWAS of y1 for males
        mt_sampled_m = mt_sampled.filter_cols(mt_sampled.sex == 0)
        gwas_y1_m = gwas(mt_sampled_m.y[3], mt_sampled_m.GT.n_alt_alleles(), [1.0])
        fn = output_bucket + 'gwas/gwas_y1_m_sample' + str(i) + '.tsv'
        export_gwas(gwas_y1_m, fn)
