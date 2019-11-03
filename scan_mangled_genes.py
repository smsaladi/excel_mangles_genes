"""
This file scans for any mangled gene identifiers 

using regexes from 
https://sourceforge.net/projects/genenameerrorsscreen/

corresponding paper:
Mark Ziemann, Yotam Eren & Assam El-Osta 
Genome Biology volume 17, Article number: 177 (2016)
https://doi.org/10.1186/s13059-016-1044-7

"""

import re
import datetime

import click

import numpy as np
import pandas as pd
import xlrd

mangled_re = [
    # XX/XX/XXxx
    '[0-9]{1,2}\/[0-9]{1,2}\/[0-9]{1,4}',
    # XX-XX-XXxx
    '[0-9]{1,2}-[0-9]{1,2}-[0-9]{1,4}',
    # D-MMM or DD-MMM
    '[0-9]{1,2}-(?:JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC|Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec|jan|feb|mar|apr|may|jun|jul|aug|sep|oct|nov|dec)',

    # Scientific notation
    '[0-9]\.[0-9]{2}E\+[0-9]{2}'
]


# combine into a single large regex
mangled_re = ['(?:' + x + ')' for x in mangled_re]
re_all = '|'.join(mangled_re)

def has_mangled(ser):
    detect_re = ser.str.contains(re_all)
    detect_dt = ser.apply(lambda x: isinstance(x, datetime.datetime))

    # if datetime, then regex balks
    detect_re.fillna(value=False, inplace=True)

    return detect_re | detect_dt

def select_gene_cols(df, genes, thresh=0.2):
    df = df.select_dtypes('object')
    
    def count_names(ser):
        if ser.size == 0:
            return 0
        return ser.isin(genes).sum() / ser.size

    # find columns with known gene names
    gene_cols = df.apply(count_names, axis=0)
    gene_cols = gene_cols[gene_cols > thresh]

    # check if genes are in rows
    if len(gene_cols) == 0:
        gene_rows = df.apply(count_names, axis=1)
        gene_rows = gene_rows[gene_rows > thresh]
        return df.loc[gene_rows.index].T

    return df[gene_cols.index]

def check_df(df):
    df.fillna('', inplace=True)
    df_detect = df.apply(has_mangled, axis=0)

    # all non-matches will be set to `nan`
    found = df[df_detect]
    found = found.values[~pd.isnull(found)]
    found = found.tolist()

    return df_detect.values.sum(), found


def read_file(fn):
    if fn.endswith('csv'):
        df = pd.read_csv(fn, low_memory=False)
        all_df = {'csv': df}
    else:
        try:
            all_df = pd.read_excel(fn, sheet_name=None)
        except xlrd.biffh.XLRDError:
            print("xlrd unable to read file: " + fn, file=sys.stderr)
            all_df = {}

    return all_df

@click.command()
@click.argument('fns', nargs=-1)
@click.option('--refnames', default='all_gene_ids.ziemann.txt')
def check_all_files(fns, refnames):
    genes = pd.read_csv(refnames, usecols=[0], header=None)
    genes = genes[0]

    for f in fns:
        all_df = read_file(f)
        for name, df in all_df.items():
            df = select_gene_cols(df, genes)
            count, vals = check_df(df)
            print(f, name, count, vals)

    return

if __name__ == '__main__':
    check_all_files()

