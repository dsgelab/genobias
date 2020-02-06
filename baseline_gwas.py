#!/usr/bin/env python
# coding: utf-8

import hail as hl
import hail.expr.aggregators as agg
import hail.methods

import requests

import pandas as pd

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


mt = hl.read_matrix_table('gs://mattia-simulations/simEUR350.mt')
output_bucket = 'gs://mattia-simulations/EUR350_base/'

# Export phenotypes
mt = mt.annotate_cols(y0=mt.y[0])
mt = mt.annotate_cols(y1=mt.y[1])
mt = mt.annotate_cols(y2=mt.y[2])
mt = mt.annotate_cols(y3=mt.y[3])

mt.cols().select('s',
                 'sex',
                 'y0',
                 'y1',
                 'y2',
                 'y3').key_by().export(
                 output_bucket + 'phenotypes/pheno_0.tsv')


output_bucket = 'gs://mattia-simulations/EUR350_base/10/'

# ------------ sex -------------
# GWAS of sex
gwas_s = gwas(mt.sex, mt.GT.n_alt_alleles(), [1.0])
fn = output_bucket + 'gwas/gwas_sex_0.tsv'
export_gwas(gwas_s, fn)

# ------------ y0 --------------
# GWAS of y0
gwas_y0 = gwas(mt.y[0], mt.GT.n_alt_alleles(), [1.0])
fn = output_bucket + 'gwas/gwas_y0_0.tsv'
export_gwas(gwas_y0, fn)

# GWAS of y0 adjusted for sex
gwas_y0_adj = gwas(mt.y[0], mt.GT.n_alt_alleles(), [1.0, mt.sex])
fn = output_bucket + 'gwas/gwas_y0_0_adj.tsv'
export_gwas(gwas_y0_adj, fn)

# GWAS of y0 for females
mt_f = mt.filter_cols(mt.sex == 1)
gwas_y0_f = gwas(mt_f.y[0], mt_f.GT.n_alt_alleles(), [1.0])
fn = output_bucket + 'gwas/gwas_y0_f_0.tsv'
export_gwas(gwas_y0_f, fn)

# GWAS of y0 for males
mt_m = mt.filter_cols(mt.sex == 0)
gwas_y0_m = gwas(mt_m.y[0], mt_m.GT.n_alt_alleles(), [1.0])
fn = output_bucket + 'gwas/gwas_y0_m_0.tsv'
export_gwas(gwas_y0_m, fn)

# ------------ y1 ---------------
# GWAS of y1
gwas_y1 = gwas(mt.y[1], mt.GT.n_alt_alleles(), [1.0])
fn = output_bucket + 'gwas/gwas_y1_0.tsv'
export_gwas(gwas_y1, fn)

# GWAS of y1 adjusted for sex
gwas_y1_adj = gwas(mt.y[1], mt.GT.n_alt_alleles(), [1.0, mt.sex])
fn = output_bucket + 'gwas/gwas_y1_0_adj.tsv'
export_gwas(gwas_y1_adj, fn)

# GWAS of y1 for females
mt_f = mt.filter_cols(mt.sex == 1)
gwas_y1_f = gwas(mt_f.y[1], mt_f.GT.n_alt_alleles(), [1.0])
fn = output_bucket + 'gwas/gwas_y1_f_0.tsv'
export_gwas(gwas_y1_f, fn)

# GWAS of y1 for males
mt_m = mt.filter_cols(mt.sex == 0)
gwas_y1_m = gwas(mt_m.y[1], mt_m.GT.n_alt_alleles(), [1.0])
fn = output_bucket + 'gwas/gwas_y1_m_0.tsv'
export_gwas(gwas_y1_m, fn)

output_bucket = 'gs://mattia-simulations/EUR350_base/30/'

# ------------ sex -------------
# GWAS of sex
gwas_s = gwas(mt.sex, mt.GT.n_alt_alleles(), [1.0])
fn = output_bucket + 'gwas/gwas_sex_0.tsv'
export_gwas(gwas_s, fn)

# ------------ y0 --------------
# GWAS of y0
gwas_y0 = gwas(mt.y[2], mt.GT.n_alt_alleles(), [1.0])
fn = output_bucket + 'gwas/gwas_y0_0.tsv'
export_gwas(gwas_y0, fn)

# GWAS of y0 adjusted for sex
gwas_y0_adj = gwas(mt.y[2], mt.GT.n_alt_alleles(), [1.0, mt.sex])
fn = output_bucket + 'gwas/gwas_y0_0_adj.tsv'
export_gwas(gwas_y0_adj, fn)

# GWAS of y0 for females
mt_f = mt.filter_cols(mt.sex == 1)
gwas_y0_f = gwas(mt_f.y[2], mt_f.GT.n_alt_alleles(), [1.0])
fn = output_bucket + 'gwas/gwas_y0_f_0.tsv'
export_gwas(gwas_y0_f, fn)

# GWAS of y0 for males
mt_m = mt.filter_cols(mt.sex == 0)
gwas_y0_m = gwas(mt_m.y[2], mt_m.GT.n_alt_alleles(), [1.0])
fn = output_bucket + 'gwas/gwas_y0_m_0.tsv'
export_gwas(gwas_y0_m, fn)

# ------------ y1 ---------------
# GWAS of y1
gwas_y1 = gwas(mt.y[3], mt.GT.n_alt_alleles(), [1.0])
fn = output_bucket + 'gwas/gwas_y1_0.tsv'
export_gwas(gwas_y1, fn)

# GWAS of y1 adjusted for sex
gwas_y1_adj = gwas(mt.y[3], mt.GT.n_alt_alleles(), [1.0, mt.sex])
fn = output_bucket + 'gwas/gwas_y1_0_adj.tsv'
export_gwas(gwas_y1_adj, fn)

# GWAS of y1 for females
mt_f = mt.filter_cols(mt.sex == 1)
gwas_y1_f = gwas(mt_f.y[3], mt_f.GT.n_alt_alleles(), [1.0])
fn = output_bucket + 'gwas/gwas_y1_f_0.tsv'
export_gwas(gwas_y1_f, fn)

# GWAS of y1 for males
mt_m = mt.filter_cols(mt.sex == 0)
gwas_y1_m = gwas(mt_m.y[3], mt_m.GT.n_alt_alleles(), [1.0])
fn = output_bucket + 'gwas/gwas_y1_m_0.tsv'
export_gwas(gwas_y1_m, fn)
