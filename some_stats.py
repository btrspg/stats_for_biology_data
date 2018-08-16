#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2018/6/2 下午1:30
# @Author  : Chen Yuelong
# @Mail    : yuelong_chen@yahoo.com
# @File    : some_stats.py
# @Software: PyCharm

from __future__ import absolute_import, unicode_literals

__author__ = 'Chen Yuelong'
import argparse
from commons.pipeline import stat_pipeline, template_gc_length_pipeline, divide_pipeline, \
    pyecharts_3d_scatter_pipeline, cal_bed_gc


def test(args):
    print(args)


def main():
    parser = argparse.ArgumentParser(description='Some stats for some biology Format file')
    parser_sub = parser.add_subparsers(dest='subparser_name', help='Sub-commands (use with -h for more info)')
    # parser_sub.required=True
    amplicon_stats_parser = parser_sub.add_parser('Amplicon_stats', help='统计扩增子效率')
    amplicon_stats_parser.add_argument('--alignFile', dest='align', required=True, action='store', nargs='+',
                                       help='Bam/BAMS,<required>')
    amplicon_stats_parser.add_argument('--threshold', dest='threshold', required=True, type=float, action='store',
                                       help='threshold,<required>')
    amplicon_stats_parser.add_argument('--target', dest='bed', required=True, action='store',
                                       help='target region,<required>')
    amplicon_stats_parser.add_argument('--prefix', dest='outdir', required=True, action='store',
                                       help='result outdir prefix,<required>')
    amplicon_stats_parser.set_defaults(func=stat_pipeline)

    bam_template_gc_parser = parser_sub.add_parser('Bam_GC', help='统计DNA模版GC含量')
    '''align,fasta,outdir,sample_n=10000'''
    bam_template_gc_parser.add_argument('--align', dest='align', action='store', required=True,
                                        help='bamfile,<required>')
    bam_template_gc_parser.add_argument('--fasta', dest='fasta', action='store', required=True,
                                        help='fasta,need fasta.fai in same direction,<required>')
    bam_template_gc_parser.add_argument('--outdir', dest='outdir', action='store', required=True,
                                        help='outdir,<required>')
    bam_template_gc_parser.add_argument('--sample', dest='sample_n', action='store', type=int,
                                        help='每个contig抽取的模版数（default=100000）', default=100000)
    bam_template_gc_parser.set_defaults(func=template_gc_length_pipeline)

    divide_parser = parser_sub.add_parser('Check', help='根据数据获得样本的可分性')

    divide_parser.add_argument('--file', dest='file', action='store', required=True,
                               help='csv file,<required>')
    divide_parser.add_argument('--method', dest='method', choices=['tsne'], action='store', required=True,
                               help='区分的方法，目前只有T-SNE(pca等后续会加入),<required>')
    divide_parser.add_argument('--outdir', dest='outdir', action='store', required=True,
                               help='outdir,<required>')
    divide_parser.set_defaults(func=divide_pipeline)

    plot_parser = parser_sub.add_parser('plot3Dscatter', help='画3D scatter图')

    plot_parser.add_argument('--file', dest='file', action='store', required=True,
                             help='csv file（label,c1,c2,c3）,<required>')
    plot_parser.add_argument('--output', dest='output', action='store', required=True,
                             help='output,结尾需要是.html,<required>')
    plot_parser.set_defaults(func=pyecharts_3d_scatter_pipeline)

    calgc_parser = parser_sub.add_parser('Bed_GC', help='统计bed区域GC含量')

    calgc_parser.add_argument('--bed', dest='bed', action='store', required=True,
                              help='bed file,<required>')
    calgc_parser.add_argument('--fasta', dest='fasta', action='store', required=True,
                              help='参考基因组，需要有samtools构建的index在同级目录,<required>')
    calgc_parser.add_argument('--output', dest='output', action='store', required=True,
                              help='output,结尾需要.csv,<required>')
    calgc_parser.set_defaults(func=cal_bed_gc)

    args = parser.parse_args()
    try:
        args.func(args)
    except Exception as e:
        print(e.args)
        parser.print_help()


if __name__ == '__main__':
    main()
