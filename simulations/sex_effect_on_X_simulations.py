#!/usr/bin/env python
# coding: utf-8

import hail as hl
import numpy as np
from math import exp, log
import random

hl.init()


# # # Methods to run GWAS and export results
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


# # # Import MatrixTable with simulated data
mt0 = hl.read_matrix_table('gs://mattia/mattia-simulations/simEUR350_2.mt')


# # # Select phenotype columns and output bucket
# mt.y:
# y0:y3    h2 = 0.1
# y4:y7    h2 = 0.3

# h2: 0.1
# mt = mt.annotate_cols(y0=mt.y[0])
# mt = mt.annotate_cols(y1=mt.y[1])
# out_bucket = 'gs://.../simulations_0.1/'

# h2: 0.3
mt0 = mt0.annotate_cols(y0=mt0.y[4])
mt0 = mt0.annotate_cols(y1=mt0.y[5])

# # # Randomly sample original pop to 50,100,150k (which will give sampled pop of 25,50,75)
N = [50000, 100000, 150000]

for n in N:
    out_bucket = 'gs://mattia/mattia-simulations/simulation_sample_size/' + str(n/1000) + '/'

    indices = random.sample(range(mt0.count_cols()), n)
    mt = mt0.choose_cols(list(indices))

    # Export phenotypes
    mt.cols().select('s', 'sex', 'y0', 'y1').key_by().export(out_bucket + 'phenotypes/pheno_0.tsv')

    # # # GWAS in full population
    # Unadjusted GWASs
    result_ht = hl.linear_regression_rows(y=[mt.sex,
                                             mt.y0,
                                             mt.y1],
                                          x=mt.GT.n_alt_alleles(),
                                          covariates=[1],
                                          pass_through=['rsid'])

    result_ht = result_ht.annotate(A1=result_ht.alleles[0],
                                   A2=result_ht.alleles[1]).key_by()

    file_names = [out_bucket + 'gwas/gwas_sex_0.tsv',
                  out_bucket + 'gwas/gwas_y0_0.tsv',
                  out_bucket + 'gwas/gwas_y1_0.tsv']

    for j, file_name in enumerate(file_names):
        result_ht.select(result_ht.locus,
                         result_ht.A1,
                         result_ht.A2,
                         result_ht.rsid,
                         N=result_ht.n,
                         beta=result_ht.beta[j],
                         se=result_ht.standard_error[j],
                         p_value=result_ht.p_value[j]).export(file_name)

    # Adjusted GWASs
    result_ht = hl.linear_regression_rows(y=[mt.y0,
                                             mt.y1],
                                          x=mt.GT.n_alt_alleles(),
                                          covariates=[1, mt.sex],
                                          pass_through=['rsid'])

    result_ht = result_ht.annotate(A1=result_ht.alleles[0],
                                   A2=result_ht.alleles[1]).key_by()

    file_names = [out_bucket + 'gwas/gwas_y0_0_adj.tsv',
                  out_bucket + 'gwas/gwas_y1_0_adj.tsv']

    for j, file_name in enumerate(file_names):
        result_ht.select(result_ht.locus,
                         result_ht.A1,
                         result_ht.A2,
                         result_ht.rsid,
                         N=result_ht.n,
                         beta=result_ht.beta[j],
                         se=result_ht.standard_error[j],
                         p_value=result_ht.p_value[j]).export(file_name)

    # # # Sampling and GWAS in sampled population
    # Save columns to Pandas.df for sampling
    df = mt.cols().select('s', 'y0', 'y1', 'sex').key_by().to_pandas()

    # Sampling parameters
    OR = [1.2, 1.5, 1.8, 2, 3]
    K = [-0.5, -0.3, 0, 0.3, 0.7, 1, 1.5]

    for k in K:
        for o in OR:

            # Participation bias
            df['z'] = df['y0'] * log(o) + df['y1'] * log(o) + np.random.normal(0, 0.1, len(df.index))

            # Sex-differential effect Zm = Z*K
            df.loc[(df['sex'] == 0), 'z'] = df.loc[(df['sex'] == 0), 'z'] * k
            df['prob'] = [1 / (1 + exp(-z)) for z in df['z']]
            df['sel'] = np.random.binomial(n=1, p=df['prob'], size=len(df.index))

            # Filter MatrixTable and get sample
            samples_to_keep = set(df.loc[(df['sel'] == 1), 's'])
            set_to_keep = hl.literal(samples_to_keep)
            mt_sampled = mt.filter_cols(set_to_keep.contains(mt['s']), keep=True)

            i = '_'+str(k)+'_'+str(o)

            # Export phenotypes
            mt_sampled.cols().select('s', 'sex', 'y0', 'y1').key_by().export(out_bucket + 'phenotypes/pheno' + i + '.tsv')

            # Unadjusted GWASs
            result_ht = hl.linear_regression_rows(y=[mt_sampled.sex,
                                                     mt_sampled.y0,
                                                     mt_sampled.y1],
                                                  x=mt_sampled.GT.n_alt_alleles(),
                                                  covariates=[1],
                                                  pass_through=['rsid'])

            result_ht = result_ht.annotate(A1=result_ht.alleles[0],
                                           A2=result_ht.alleles[1]).key_by()

            file_names = [out_bucket + 'gwas/gwas_sex' + i + '.tsv',
                          out_bucket + 'gwas/gwas_y0' + i + '.tsv',
                          out_bucket + 'gwas/gwas_y1' + i + '.tsv']

            for j, file_name in enumerate(file_names):
                result_ht.select(result_ht.locus,
                                 result_ht.A1,
                                 result_ht.A2,
                                 result_ht.rsid,
                                 N=result_ht.n,
                                 beta=result_ht.beta[j],
                                 se=result_ht.standard_error[j],
                                 p_value=result_ht.p_value[j]).export(file_name)

            # Adjusted GWASs
            result_ht = hl.linear_regression_rows(y=[mt_sampled.y0,
                                                     mt_sampled.y1],
                                                  x=mt_sampled.GT.n_alt_alleles(),
                                                  covariates=[1, mt_sampled.sex],
                                                  pass_through=['rsid'])

            result_ht = result_ht.annotate(A1=result_ht.alleles[0],
                                           A2=result_ht.alleles[1]).key_by()

            file_names = [out_bucket + 'gwas/gwas_y0' + i + '_adj.tsv',
                          out_bucket + 'gwas/gwas_y1' + i + '_adj.tsv']

            for j, file_name in enumerate(file_names):
                result_ht.select(result_ht.locus,
                                 result_ht.A1,
                                 result_ht.A2,
                                 result_ht.rsid,
                                 N=result_ht.n,
                                 beta=result_ht.beta[j],
                                 se=result_ht.standard_error[j],
                                 p_value=result_ht.p_value[j]).export(file_name)
