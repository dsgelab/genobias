#!/usr/bin/env python
# coding: utf-8

import hail as hl
import numpy as np
from math import exp, log
import subprocess

hl.init()

# Install statsmodels (for glm)
subprocess.call(['/opt/conda/miniconda3/bin/pip', 'install', 'statsmodels'])
import statsmodels.api as sm
import statsmodels.formula.api as smf


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


# input_matrix = 'gs://mattia/mattia-simulations/simEUR350_2.mt'
# output_bucket = 'gs://mattia/mattia-simulations/simulations_heckman_rg_0/'
input_matrix = 'gs://mattia/mattia-simulations/simEUR350_correlated.mt'
output_bucket = 'gs://mattia/mattia-simulations/simulations_heckman/'

# Phenotypes:
# simEUR350_2 (uncorrelated)
# y0:y3     h2 = 0.1
# y04:y7    h2 = 0.3
# x = 4
# y = 5
# u = 6

# simEUR350_correlated:
# y0        h2 = 0.3
# y1,y8     h2 = 0.3, rg(y0,yn) = [-0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5, 0]
x = 0
u = 8

rgs = [-0.3, -0.1, 0, 0.1, 0.3]

# Import matrix and annotate cols with phenotypes
mt = hl.read_matrix_table(input_matrix)

mt = mt.annotate_cols(U=mt.y[u])
mt = mt.annotate_cols(X=mt.y[x])

for i in range(5):
    mt = mt.annotate_cols(Y=mt.y[i+2])

    # GWAS of X, Y in all
    result_ht = hl.linear_regression_rows(y=[mt.X,
                                             mt.Y],
                                          x=mt.GT.n_alt_alleles(),
                                          covariates=[1],
                                          pass_through=['rsid'])

    result_ht = result_ht.annotate(A1=result_ht.alleles[0],
                                   A2=result_ht.alleles[1]).key_by()

    file_names = [output_bucket + 'gwas/gwas_X_' + str(rgs[i]) + '.tsv',
                  output_bucket + 'gwas/gwas_Y_' + str(rgs[i]) + '.tsv']

    for j, file_name in enumerate(file_names):
        result_ht.select(result_ht.locus,
                         result_ht.A1,
                         result_ht.A2,
                         result_ht.rsid,
                         N=result_ht.n,
                         beta=result_ht.beta[j],
                         se=result_ht.standard_error[j],
                         p_value=result_ht.p_value[j]).export(file_name)

    # Save mt cols to table and convert to pandas df for sampling
    df = mt.cols().select('s', 'X', 'Y', 'U', 'sex').key_by().to_pandas()

    # Selection based on X,Y,U no sex-diff bias
    OR = [1.2, 1.5, 1.8, 2, 3]

    for o in OR:
        df['z'] = df['X'] * log(o) + df['Y'] * log(o) + df['U'] * log(2)
        df['prob'] = [1 / (1 + exp(-z)) for z in df['z']]
        df['sel'] = np.random.binomial(n=1, p=df['prob'], size=len(df.index))

        name = '_' + str(rgs[i]) + '_' + str(o)

        # Prediction
        formula = 'sel ~ X + U'
        mod1 = smf.glm(formula=formula, data=df, family=sm.families.Binomial(sm.families.links.probit)).fit()
        print(mod1.summary())
        pred = mod1.predict(df)
        df['pred'] = pred

        # Save phenotypes table
        t = hl.Table.from_pandas(df).key_by('s')
        t.export(output_bucket + 'phenotypes/pheno'+name+'.tsv')

        # Join 'phenotypes' table and matrixtable
        mt = mt.annotate_cols(pheno=t[mt.s])

        # GWAS of sel
        gwas_sel = gwas(mt.pheno.sel, mt.GT.n_alt_alleles(), [1.0])
        fn = output_bucket + 'gwas/gwas_sel'+name+'.tsv'
        export_gwas(gwas_sel, fn)

        # GWAS of X*, Y*, pred*
        mt_sel = mt.filter_cols(mt.pheno.sel == 1)

        result_ht = hl.linear_regression_rows(y=[mt_sel.pheno.pred,
                                                 mt_sel.pheno.X,
                                                 mt_sel.pheno.Y],
                                              x=mt_sel.GT.n_alt_alleles(),
                                              covariates=[1],
                                              pass_through=['rsid'])

        result_ht = result_ht.annotate(A1=result_ht.alleles[0],
                                       A2=result_ht.alleles[1]).key_by()

        file_names = [output_bucket + 'gwas/gwas_pred_star' + name + '.tsv',
                      output_bucket + 'gwas/gwas_X_star' + name + '.tsv',
                      output_bucket + 'gwas/gwas_Y_star' + name + '.tsv']

        for j, file_name in enumerate(file_names):
            result_ht.select(result_ht.locus,
                             result_ht.A1,
                             result_ht.A2,
                             result_ht.rsid,
                             N=result_ht.n,
                             beta=result_ht.beta[j],
                             se=result_ht.standard_error[j],
                             p_value=result_ht.p_value[j]).export(file_name)
