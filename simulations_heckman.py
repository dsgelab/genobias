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
output_bucket = 'gs://mattia-simulations/simulation_heckman/'

# Save mt cols to table and convert to pandas df for sampling
# Phenotypes:
# y0:y3     h2 = 0.1
# y04:y7    h2 = 0.3
# Gonna use y04,y05,y06 here as x,y,Z
mt = mt.annotate_cols(X=mt.y[4])
mt = mt.annotate_cols(Y=mt.y[5])
mt = mt.annotate_cols(Z=mt.y[6])

df = mt.cols().select('s', 'X', 'Y', 'Z', 'sex').key_by().to_pandas()

# Selection based on X,Y,Z
df['z'] = df['X'] * log(3) + df['Y'] * log(3) + df['Z'] * log(1.5)
df['prob'] = [1 / (1 + exp(z)) for z in df['z']]
df['select'] = np.random.binomial(n=1, p=df['prob'], size=len(df.index))

samples_to_keep = set(df.loc[(df['select'] == 1), 's'])
set_to_keep = hl.literal(samples_to_keep)

# Add 'sel' column to mt
mt = mt.annotate_cols(sel=hl.cond(set_to_keep.contains(mt['s']), 1, 0))

# Export phenotypes
mt.cols().select('s',
                 'sex',
                 'X',
                 'Y',
                 'Z',
                 'sel').key_by().export(output_bucket + 'phenotypes/pheno.tsv')

# GWAS of sel
gwas_sel = gwas(mt.sel, mt.GT.n_alt_alleles(), [1.0])
fn = output_bucket + 'gwas_sel.tsv'
export_gwas(gwas_sel, fn)

# GWAS of X,Y,Z in all samples
gwas_x = gwas(mt.X, mt.GT.n_alt_alleles(), [1.0])
fn = output_bucket + 'gwas_x_all.tsv'
export_gwas(gwas_x, fn)

gwas_y = gwas(mt.Y, mt.GT.n_alt_alleles(), [1.0])
fn = output_bucket + 'gwas_y_all.tsv'
export_gwas(gwas_y, fn)

gwas_z = gwas(mt.Z, mt.GT.n_alt_alleles(), [1.0])
fn = output_bucket + 'gwas_z_all.tsv'
export_gwas(gwas_z, fn)

# GWAS of X,Y,Z in selected
mt_sel = mt.filter_cols(mt.sel == 1)

gwas_x = gwas(mt_sel.X, mt_sel.GT.n_alt_alleles(), [1.0])
fn = output_bucket + 'gwas_x_sel.tsv'
export_gwas(gwas_x, fn)

gwas_y = gwas(mt_sel.Y, mt_sel.GT.n_alt_alleles(), [1.0])
fn = output_bucket + 'gwas_y_sel.tsv'
export_gwas(gwas_y, fn)

gwas_z = gwas(mt_sel.Z, mt_sel.GT.n_alt_alleles(), [1.0])
fn = output_bucket + 'gwas_z_sel.tsv'
export_gwas(gwas_z, fn)
