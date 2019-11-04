"""
This file scans for any mangled gene identifiers

using regexes from
https://sourceforge.net/projects/genenameerrorsscreen/

corresponding paper:
Mark Ziemann, Yotam Eren & Assam El-Osta
Genome Biology volume 17, Article number: 177 (2016)
https://doi.org/10.1186/s13059-016-1044-7

"""

import sys
import re
import datetime

import click

import numpy as np
import pandas as pd

mangled_re = [
    # XX/XX/XXxx
    '[0-9]{1,2}\/[0-9]{1,2}\/[0-9]{1,4}',
    # mm-dd-YYyy
    '[0,1]{1,2}-[0-9]{1,2}-\d{2}(?:\d{2})?',
    # dd-mm-YYyy

    # XX-XX-XXxx
    '[0-9]{1,2}-[0-9]{1,2}-\d{2}(?:\d{2})?',
    # D-MMM or DD-MMM
    '[0-9]{1,2}-(?:JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC|[Jj]an|[Ff]eb|[Mm]ar|[Aa]pr|[Mm]ay|[Jj]un|[Jj]ul|[Aa]ug|[Ss]ep|[Oo]ct|[Nn]ov|[Dd]ec)$',

    # Scientific notation
    '[0-9]\.[0-9]{2}E\+[0-9]{2}'
]


# combine into a single large regex
mangled_re = ['(?:' + x + ')' for x in mangled_re]
# allow spaces on either side but nothing else
re_all = '|'.join(mangled_re)
re_all = re.compile(re_all)

def not_manged(ser):
    return ser.str.contains('Date')

def has_mangled(ser, maxcell=20):
    """Detects if a datetime/date-like string is found
    """
    # easy when excel stores cell as datetime
    detect_dt = ser.apply(lambda x: isinstance(x, datetime.datetime))

    # pd.Series.str balks if non-string types are mixed in
    ser = ser.astype(str)
    # if the cell is long, probably not a mangled genename
    ser = ser.mask(ser.str.len() > maxcell, "")

    # detect dates as strings
    detect_re = ser.str.contains(re_all)

    return detect_re | detect_dt

def count_names(ser, genes):
    if ser.size == 0:
        return 0
    return ser.isin(genes).sum() / ser.size

def stringlike(x):
    try:
        x.astype(str)
        return True
    except:
        return False

def select_gene_cols(df, genes, mask=True, thresh=0.2):
    """Selects columns with gene names,
       * can also mask real gene names at the same time

       * Assumes genes are unique and lowercase
    """

    # multiple levels to select string/datetime columns
    df = df.select_dtypes('object')

    # find columns with known gene names
    df_str = df.apply(lambda x: x.astype(str).str.upper())
    df_is_gene = df_str.isin(genes)
    gene_cols = df_str.columns[df_is_gene.sum(axis=0) > thresh]
    
    df = df[gene_cols].copy()
    # so values are not detected by regex
    if mask:
        df.mask(df_is_gene, "", inplace=True)

    return df

def is_longer(df):
    return df.shape[0] > df.shape[1]

def check_df(df, genes):
    """Checks a dataframe for manged gene names
    """
    df_sub = select_gene_cols(df, genes)
    if df_sub.size == 0:
        df_sub = select_gene_cols(df.T, genes)
        df_sub = df_sub.T

    df_sub.fillna('', inplace=True)

    # apply over rows/columns (depending on which is fewer)
    df_detect = df_sub.apply(has_mangled, axis=int(is_longer(df)))

    # all non-matches will be set to `nan`
    found = df[df_detect]
    found = found.values[(~pd.isnull(found)).values]
    found = found.tolist()

    return df_detect.values.sum(), found


def read_file(fn):
    if fn.endswith('csv'):
        try:
           df = pd.read_csv(fn, low_memory=False)
           all_df = {'csv': df}
        except Exception as e:
            print(fn + ": Unable to parse." + str(e), file=sys.stderr)
            all_df = {}
    else:
        try:
            all_df = pd.read_excel(fn, sheet_name=None)
        except Exception as e:
            print(fn + ": Unable to read file. " + str(e), file=sys.stderr)
            all_df = {}

    all_df = {k:v for k,v in all_df.items() if v.size > 0}

    return all_df

@click.command()
@click.argument('fns', nargs=-1)
@click.option('--refnames', default='all_gene_ids.ziemann.txt')
def check_all_files(fns, refnames):
    genes = pd.read_csv(refnames, usecols=[0], header=None)
    genes = genes[0].str.upper().unique()

    for f in fns:
        print("# ", f, file=sys.stderr)
        all_df = read_file(f)
        for name, df in all_df.items():
            count, vals = check_df(df, genes)
            print(f, name, count, vals)

    return

if __name__ == '__main__':
    check_all_files()

