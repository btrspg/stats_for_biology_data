#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2018/6/27 下午5:25
# @Author  : Chen Yuelong
# @Mail    : yuelong_chen@yahoo.com
# @File    : pipeline.py
# @Software: PyCharm

from __future__ import absolute_import, unicode_literals
import pandas as pd

__author__ = 'Chen Yuelong'
import os, sys
sys.path.append(os.path.abspath(os.path.dirname(os.path.dirname(__file__))))
from plots.imageshow import echart3D_tsne
from commons.utils import *
from formats.format_utils import get_gc_length_count,get_bed_gc
from formats.formats import Fasta

def cal_bed_gc(args):
    '''

    :param args:
    :return:
    '''
    fasta = Fasta(args.fasta)
    get_bed_gc(args.bed,fasta,args.output)

def pyecharts_3d_scatter_pipeline(args):
    '''

    :param args:
    :return:
    '''
    data = pd.read_csv(args.file)
    echart3D_tsne(data,'label','c1','c2','c3',args.output)

def divide_pipeline(args):
    '''

    :param args:
    :return:
    '''
    if args.method == 'tsne':
        divid_method(trans_to_tsne,args.file,args.outdir)
    else:
        raise ValueError('method can only be tsne')

def stat_pipeline(args):
    '''

    :param args:
    :return:
    '''
    samfile = args.align
    target =args.bed
    threshold = args.threshold
    outdir = args.outdir
    check_dirs(outdir)
    cache_list=[]
    name_list=[]
    with open(target,'r') as f:
        for sam in samfile:
            name_list.append(os.path.basename(sam))
            with pysam.AlignmentFile(sam,'rb') as align:
                f.seek(0,0)
                cache_list.append(stat_amplicon(align,f,threshold))

    pd_data = pd.DataFrame(cache_list).transpose()
    pd_data.columns=name_list
    pd_data.to_csv('{}.stat.csv'.format(outdir))
    pd_data.plot(logy=True, kind='bar', figsize=(max(50, pd_data.shape[0]/10), 10))
    plt.savefig('{}.stat.png'.format(outdir), dpi=300)
    print('Finish!!')

def template_gc_length_pipeline(args):
    '''

    :param args:
    :return:
    '''
    fasta = Fasta(args.fasta)
    check_files(args.fasta,args.align)
    check_dirs(args.outdir)
    with pysam.AlignmentFile(args.align, 'rb') as align:
        get_gc_length_count(align, fasta, args.outdir,args.sample_n)

def main():
    pass


if __name__ == '__main__':
    main()