#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2018/6/2 上午10:18
# @Author  : Chen Yuelong
# @Mail    : yuelong_chen@yahoo.com
# @File    : utils.py
# @Software: PyCharm

from __future__ import absolute_import, unicode_literals

__author__ = 'Chen Yuelong'
import os,sys
import pysam
import pandas as pd
import matplotlib
from sklearn.manifold import TSNE
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.append(os.path.abspath(os.path.dirname(os.path.dirname(__file__))))

from plots.imageshow import echart3D_tsne

def trans_to_tsne(data,n_components=3):
    '''

    :param data:
    :param n_components:
    :return:
    '''
    label = data.pop('label')
    x_tsne = TSNE(n_components=n_components).fit_transform(data)
    x_tsne=pd.DataFrame(x_tsne,columns=['c1','c2','c3'])
    # x_tsne.columns=['c1','c2','c3']
    return pd.concat([label,x_tsne],axis=1)


def cal_overlap(dna_start,dna_end,target_start,target_end):
    '''

    :param dna_start:
    :param dna_end:
    :param target_start:
    :param target_end:
    :return:
    '''

    return (min(dna_end,int(target_end))-max(dna_start,int(target_start)))/(int(target_end)-int(target_start))

def get_overlap_amp(ctg,alignsegment,target_start,target_end,cache,threshold):
    '''

    :param ctg:
    :param alignsegment:
    :param target_start:
    :param target_end:
    :param cache:
    :param threshold:
    :return:
    '''
    cache.setdefault(ctg, 0)
    for seg in alignsegment:
        if seg.template_length > 0 and cal_overlap(seg.reference_start,seg.reference_start+seg.template_length,
                                                   target_start,target_end) > threshold:
            cache[ctg] += 1
    return cache

def get_samfile(file):
    '''

    :param file:
    :return:
    '''
    return pysam.AlignmentFile(file,'rb')

def stat_amplicon(align,buffer,threshold):
    '''

    :param samfile:
    :param target:
    :param threshold:
    :return:
    '''
    cache={}
    for line in buffer.readlines():
        # print(line)
        cells = line.strip('\n').split('\t')
        cache = get_overlap_amp(line.strip('\n').replace('\t','-'),align.fetch(cells[0],int(cells[1]),int(cells[2])),
                                int(cells[1]),int(cells[2]),cache,threshold)
    # print(cache)
    return cache






def divid_method(function,input,outdir):
    '''

    :param function:
    :param input:
    :param outdir:
    :return:
    '''
    pd_data = pd.read_csv(input)
    # pd_data.to_csv()
    trans_data = function(pd_data)
    trans_data.to_csv('{outdir}/TSNE-result.csv'.format(outdir=outdir))
    echart3D_tsne(trans_data,'label','c1','c2','c3','{outdir}/TSNE.html'.format(outdir=outdir))











def check_dirs(*dirs):
    '''

    :param dirs:
    :return:
    '''
    for di in dirs:
        if not os.path.exists(di):
            os.makedirs(di,exist_ok=True)

def check_files(*files):
    '''

    :param files:
    :return:
    '''
    not_exists = []
    for file in files:
        if not os.path.exists(file):
            not_exists.append(file)
    if len(not_exists) == 0:
        print(','.join(files)+' all found!')
        return True
    else:
        err = '{} not Found!'.format(','.join(not_exists))
        raise FileNotFoundError(err)


def main():
    pass



if __name__ == '__main__':
    main()