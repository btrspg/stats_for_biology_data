#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2018/6/4 下午5:36
# @Author  : Chen Yuelong
# @Mail    : yuelong_chen@yahoo.com
# @File    : formats.py
# @Software: PyCharm

from __future__ import absolute_import, unicode_literals

__author__ = 'Chen Yuelong'
import os, sys

sys.path.append(os.path.abspath(os.path.dirname(os.path.dirname(__file__))))
from collections import Counter
import time
from math import atan, pi


class Fasta():
    '''
    fasta file
    '''

    def __init__(self, fasta, fai=''):
        self._fasta = fasta
        self._fai = fai if fai != '' else '{}.fai'.format(self._fasta)
        self._fai_dict = self._get_fai_()

    def _get_fai_(self):
        '''

        :return:
        '''
        tmp_dict = {}
        with open(self._fai, 'r') as f:
            line = f.readline()
            while line:
                cells = line.strip('\n').split('\t')
                contig = cells.pop(0)
                tmp_dict[contig] = cells
                line = f.readline()
        return tmp_dict

    @property
    def fai_dict(self):
        return self._fai_dict

    def _return_seq(self, fbuffer, contig, start, end):
        '''

        :param fbuffer:
        :param contig:
        :param start:
        :param end:
        :return:
        '''
        if contig in self._fai_dict:
            [length, offset, line, size] = list(map(lambda x: int(x), self._fai_dict[contig]))
            if length < end - start:
                raise ValueError('Require length "{}" is larger than contig "{}" length "{}"'.format(
                    end - start, contig, length
                ))
            return self._get_seq(fbuffer, start, end, offset, line, size)

        else:
            raise TypeError('{} not in fasta'.format(contig))

    def get_depth_correction(self, require_list, tmp):
        '''
        未完成
        :param require_list:
        :param tmp:
        :return:
        '''

    def _correct(self, gc, length):
        '''
        这个校正方法只是临时使用，并不是真是的校正方法，未完成
        :param gc:
        :param length:
        :return:
        '''
        return 1 / ((abs(gc - 50) / 20) * abs(length - 300) / 50)

    def get_gc(self, require_list):
        '''

        :param require_list:
        :return:
        '''
        with open(self._fasta, 'r') as fbuffer:
            gc_length = list(
                map(lambda x: self._gc_length(self._return_seq(fbuffer, x[0], x[1], x[2])), require_list))
        return [
            ('{rl0}:{rl1}-{rl2}'.format(
                rl0=rl[0],
                rl1=rl[1],
                rl2=rl[2]),gl[0],gl[1])
            for rl, gl in zip(require_list, gc_length)]

    def get_gc_length(self, require_list):
        '''

        :param require_list:
        :return:
        '''
        start = time.time()
        # print(len(require_list))
        with open(self._fasta, 'r') as fbuffer:
            gc_length = list(
                map(lambda x: self._gc_length(self._return_seq(fbuffer, x[0], x[1], x[2])), require_list))

        end = time.time()
        print('cal gc %.2f ms!' % ((end - start) * 1000))
        start = end
        counter = Counter(gc_length)
        print('count gc %.2f ms!' % ((time.time() - start) * 1000))
        return list(map(lambda x: [x[0], x[1], counter[x]], counter.keys()))

    def _gc_length(self, seq):
        '''

        :param seq:
        :return:
        '''
        at = seq.count('A') + seq.count('a') + seq.count('T') + seq.count('t')
        gc = seq.count('G') + seq.count('g') + seq.count('C') + seq.count('c')
        return int((gc / (gc + at)) * 100), len(seq)

    def _get_seq(self, fbuffer, start, end, offset, line, size):
        '''

        :param fbuffer:
        :param start:
        :param end:
        :param offset:
        :param line:
        :param size:
        :return:
        '''
        location = offset + int(start / line) + start - 1
        length = (int(end / line) - int(start / line)) * (size - line) + end - start + 1
        fbuffer.seek(location, 0)
        sequence = fbuffer.read(length).replace('\n', '')
        # print(sequence)
        return sequence

    def print(self):
        for i in self.__dict__:
            print(i)
            print(self.__dict__[i])


def main():
    pass


if __name__ == '__main__':
    main()

'''
20 63025520 2763883380 60 61
21 48129895 2827959352 60 61

测试例子
>hg19_knownGene_uc002ypa.3_4 range=chr21:33040784-33040891 5'pad=0 3'pad=0 strand=+ repeatMasking=none
GTCCATGAAAAAGCAGATGACTTGGGCAAAGGTGGAAATGAAGAAAGTACAAAGACAGGAAACGCTGGAAGTCGTTTGGCTTGTGGTGTAATTGGGATCGCCCAATAA
samtools提取结果
GTCCATGAAAAAGCAGATGACTTGGGCAAAGGTGGAAATGAAGAAAGTACAAAGACAGGAAACGCTGGAAGTCGTTTGGCTTGTGGTGTAATTGGGATCGCCCAATAA
我提取出的结果
GTCCATGAAAAAGCAGATGACTTGGGCAAAGGTGGAAATGAAGAAAGTACAAAGACAGGAAACGCTGGAAGTCGTTTGGCTTGTGGTGTAATTGGGATCGCCCAATAA
'''
