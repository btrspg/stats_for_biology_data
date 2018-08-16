#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2018/6/6 下午5:28
# @Author  : Chen Yuelong
# @Mail    : yuelong_chen@yahoo.com
# @File    : format_utils.py
# @Software: PyCharm

from __future__ import absolute_import, unicode_literals
import os, sys
import pysam
import pandas as pd
import time
import random

sys.path.append(os.path.abspath(os.path.dirname(os.path.dirname(__file__))))
__author__ = 'Chen Yuelong'
from formats.formats import Fasta
from plots.imageshow import plot3D, echart3D


def get_gc_length_count(align, fasta, outdir, sample_n=1000):
    '''

    :param align:
    :param fasta:
    :param outdir:
    :return:
    '''
    if not isinstance(align, pysam.AlignmentFile) or not isinstance(fasta, Fasta):
        raise TypeError('fasta or alignment!!')

    for i in fasta.fai_dict.keys():
        # print(i, 1, int(fasta.fai_dict[i][0]))
        start = time.time()
        print('contig:' + i)
        # segment = align.fetch(i,1,int(fasta.fai_dict[i][0]))
        segment = align.fetch(reference=i)
        end = time.time()
        print('get Segment use %.2f ms' % ((end - start) * 1000))
        start = end
        loc_list = []
        n = 0
        for seg in segment:
            if seg.is_paired and seg.template_length > 0 and seg.is_proper_pair:
                loc_list.append([i, seg.reference_start, seg.reference_start + seg.template_length])
                n = n + 1
                if n % 10000 == 0:
                    end = time.time()
                    print('{} W segments cost %.2f ms'.format(str(n / 10000)) % ((end - start) * 1000))
                    start = end
        if len(loc_list) == 0:
            continue
        sample_list = random.sample(loc_list, sample_n) if len(loc_list) > sample_n else loc_list.copy()
        loc_list.clear()
        # print(len(sample_list))
        end = time.time()
        print('{} W segments cost %.2f ms'.format(str(n / 10000)) % ((end - start) * 1000))
        # print('finished filter {}'.format((time.time() - start) * 1000))
        start = time.time()
        df = pd.DataFrame(fasta.get_gc_length(sample_list))
        print('finished calculate gc %.2f ms\n' % ((time.time() - start) * 1000))
        start = time.time()
        df.columns = ['gc', 'length', 'count']
        plot3D(df, 'gc', 'length', 'count', '{}/{}.png'.format(outdir, i))
        echart3D(df, 'gc', 'length', 'count', '{}/{}.html'.format(outdir, i))
        df.to_csv('{}/{}.csv'.format(outdir, i))
        print('finished! %.2f ms' % ((time.time() - start) * 1000))


def get_bed_gc(bed, fasta, output):
    '''

    :param bed:
    :param fasta:
    :param output:
    :return:
    '''
    if not isinstance(fasta, Fasta):
        raise TypeError('fasta or alignment!!')
    loc_list = []

    with open(bed, 'r') as fbed:

        line = fbed.readline()
        while line:
            cells = line.strip('\n').split('\t')
            if len(cells) > 0:
                loc_list.append([cells[0], int(cells[1]), int(cells[2])])
            line = fbed.readline()

    df = pd.DataFrame(fasta.get_gc(loc_list))
    df.columns = ['loc', 'gc', 'length']
    df.to_csv(output)


def main():
    pass


if __name__ == '__main__':
    main()
