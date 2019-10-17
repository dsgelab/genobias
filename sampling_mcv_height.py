#!/usr/bin/env python
# coding: utf-8

import hail as hl
import hail.expr.aggregators as agg
import hail.methods

import requests
from math import exp, log
import pandas as pd
import numpy as np
from scipy.stats import chi2

from typing import *

import random

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
    gann.select('A1', 'A2', 'rsid', 'n', 'beta', 'standard_error', 'p_value').export(fname)
    return


mt = hl.read_matrix_table('gs://mattia-simulations/simEUR350_mcv_height.mt')

output_bucket = 'gs://mattia-simulations/EUR350_sampling_mcv_height/'

# Subsampling with different probs for M and F and running GWAS
# Save mt cols to table and convert to pandas df for sampling
mt = mt.annotate_cols(y0=mt.y[0])
mt = mt.annotate_cols(y1=mt.y[1])
df = mt.cols().select('s', 'y0', 'y1', 'sex').key_by().to_pandas()

OR_y0 = [1.2, 1.5, 1.8, 2, 4]
OR_y1 = [1.2, 1.5, 1.8, 2, 4]

i = 0
for ory0 in OR_y0:
    for ory1 in OR_y1:
        i += 1

        df['z'] = df['y0'] * log(ory0) + df['y1'] * log(ory1)

        df['prob'] = 0

        df.loc[(df['sex'] == 0), 'prob'] = [1 / (1 + exp(-z)) for z in df.loc[(df['sex'] == 0), 'z']]
        df.loc[(df['sex'] == 1), 'prob'] = [1 / (1 + exp(z)) for z in df.loc[(df['sex'] == 1), 'z']]

        sampled = df.sample(frac=0.75, weights='prob', random_state=123)

        samples_to_keep = set(sampled['s'])
        set_to_keep = hl.literal(samples_to_keep)
        mt_sampled = mt.filter_cols(set_to_keep.contains(mt['s']), keep=True)

        # Export phenotypes
        mt_sampled.cols().select('s',
                                 'sex',
                                 'y0',
                                 'y1').key_by().export(
                                 output_bucket + 'phenotypes/pheno_sample' + str(i) + '.tsv')

        # ------------ sex -------------
        # GWAS of sex
        gwas_s = gwas(mt_sampled.sex, mt_sampled.GT.n_alt_alleles(), [1.0])
        fn = output_bucket + 'gwas/gwas_sex_sample' + str(i) + '.tsv'
        export_gwas(gwas_s, fn)

        # ------------ mcv --------------
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

        # ------------ height ---------------
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
