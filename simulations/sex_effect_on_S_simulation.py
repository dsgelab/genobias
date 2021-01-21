#!/usr/bin/env python
# coding: utf-8

import hail as hl
import numpy as np
from math import exp, log
import random

hl.init()

# # # Import MatrixTable with simulated data
mt = hl.read_matrix_table('gs://mattia/mattia-simulations/simEUR350_2.mt')

# # # Select phenotype columns and output bucket
# mt.y:
# y0:y3    h2 = 0.1
# y4:y7    h2 = 0.3

# h2: 0.1
# mt = mt.annotate_cols(y0=mt.y[0])
# mt = mt.annotate_cols(y1=mt.y[1])
# out_bucket = 'gs://.../simulations_0.1/'

# # # # # # #
# h2: 0.1 # #
# # # # # # #
mt = mt.annotate_cols(y0=mt.y[0]+mt.sex*1+hl.rand_norm(0, 0.1))
out_bucket = 'gs://mattia/mattia-simulations/simulation_sex_on_X_0.1_1/'

# Export phenotypes
mt.cols().select('s', 'sex', 'y0').key_by().export(out_bucket + 'phenotypes/pheno_0.tsv')

# Sampling parameters
OR = [1.2, 1.5, 2, 3, 5]
k = 1   # no sex diff effect

df = mt.cols().select('s', 'y0', 'sex').key_by().to_pandas()

for o in OR:
    # Participation bias
    df['z'] = df['y0'] * log(o)
    df['prob'] = [1 / (1 + exp(-z)) for z in df['z']]
    df['sel'] = np.random.binomial(n=1, p=df['prob'], size=len(df.index))

    # Filter MatrixTable and get sample
    samples_to_keep = set(df.loc[(df['sel'] == 1), 's'])
    set_to_keep = hl.literal(samples_to_keep)
    mt_sampled = mt.filter_cols(set_to_keep.contains(mt['s']), keep=True)

    i = '_'+str(k)+'_'+str(o)

    # Export phenotypes
    mt_sampled.cols().select('s', 'sex', 'y0').key_by().export(out_bucket + 'phenotypes/pheno' + i + '.tsv')

    # Unadjusted GWASs
    result_ht = hl.linear_regression_rows(y=[mt_sampled.sex,
                                             mt_sampled.y0],
                                          x=mt_sampled.GT.n_alt_alleles(),
                                          covariates=[1],
                                          pass_through=['rsid'])

    result_ht = result_ht.annotate(A1=result_ht.alleles[0],
                                   A2=result_ht.alleles[1]).key_by()

    file_names = [out_bucket + 'gwas/gwas_sex' + i + '.tsv',
                  out_bucket + 'gwas/gwas_y0' + i + '.tsv']

    for j, file_name in enumerate(file_names):
        result_ht.select(result_ht.locus,
                         result_ht.A1,
                         result_ht.A2,
                         result_ht.rsid,
                         N=result_ht.n,
                         beta=result_ht.beta[j],
                         se=result_ht.standard_error[j],
                         p_value=result_ht.p_value[j]).export(file_name)

# # # # # # #
# h2: 0.3 # #
# # # # # # #
mt = mt.annotate_cols(y0=mt.y[4] + mt.sex * 0.5 + hl.rand_norm(0, 0.1))
out_bucket = 'gs://mattia/mattia-simulations/simulation_sex_on_X_0.3_1/'

# Export phenotypes
mt.cols().select('s', 'sex', 'y0').key_by().export(out_bucket + 'phenotypes/pheno_0.tsv')

# Sampling parameters
OR = [1.2, 1.5, 2, 3, 5]
k = 1  # no sex diff effect

df = mt.cols().select('s', 'y0', 'sex').key_by().to_pandas()

for o in OR:
    # Participation bias
    df['z'] = df['y0'] * log(o)
    df['prob'] = [1 / (1 + exp(-z)) for z in df['z']]
    df['sel'] = np.random.binomial(n=1, p=df['prob'], size=len(df.index))

    # Filter MatrixTable and get sample
    samples_to_keep = set(df.loc[(df['sel'] == 1), 's'])
    set_to_keep = hl.literal(samples_to_keep)
    mt_sampled = mt.filter_cols(set_to_keep.contains(mt['s']), keep=True)

    i = '_' + str(k) + '_' + str(o)

    # Export phenotypes
    mt_sampled.cols().select('s', 'sex', 'y0').key_by().export(out_bucket + 'phenotypes/pheno' + i + '.tsv')

    # Unadjusted GWASs
    result_ht = hl.linear_regression_rows(y=[mt_sampled.sex,
                                             mt_sampled.y0],
                                          x=mt_sampled.GT.n_alt_alleles(),
                                          covariates=[1],
                                          pass_through=['rsid'])

    result_ht = result_ht.annotate(A1=result_ht.alleles[0],
                                   A2=result_ht.alleles[1]).key_by()

    file_names = [out_bucket + 'gwas/gwas_sex' + i + '.tsv',
                  out_bucket + 'gwas/gwas_y0' + i + '.tsv']

    for j, file_name in enumerate(file_names):
        result_ht.select(result_ht.locus,
                         result_ht.A1,
                         result_ht.A2,
                         result_ht.rsid,
                         N=result_ht.n,
                         beta=result_ht.beta[j],
                         se=result_ht.standard_error[j],
                         p_value=result_ht.p_value[j]).export(file_name)
