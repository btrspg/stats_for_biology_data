#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2018/6/8 下午11:01
# @Author  : Chen Yuelong
# @Mail    : yuelong_chen@yahoo.com
# @File    : imageshow.py
# @Software: PyCharm

from __future__ import absolute_import, unicode_literals

__author__ = 'Chen Yuelong'
import os, sys
sys.path.append(os.path.abspath(os.path.dirname(os.path.dirname(__file__))))
from mpl_toolkits.mplot3d import Axes3D
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pyecharts import Scatter3D
import json


def echart3D(pddata, x_title, y_title, z_title, output):
    pddata[pddata[y_title] > 1000] = 0
    data = json.loads(pddata[[x_title, y_title, z_title]].to_json(orient='values'))
    # print(type(data))
    # print(data)
    # data =random.sample(data,100)

    range_color = [
        '#313695', '#4575b4', '#74add1', '#abd9e9', '#e0f3f8', '#ffffbf',
        '#fee090', '#fdae61', '#f46d43', '#d73027', '#a50026']
    scatter3D = Scatter3D("3D GC+LENGTH+COUNT", width=1200, height=600)
    scatter3D.add("gc,length,count",
                  data, is_visualmap=True,
                  visual_range_color=range_color,
                  xaxis3d_name='GC(%)',
                  yaxis3d_name='Template Length',
                  zaxis3d_name='Count', )
    scatter3D.render(output)


def echart3D_tsne(pddata, label, x_title, y_title, z_title, output):
    scatter3D = Scatter3D("T-SNE", width=1200, height=600)
    for lab in set(pddata[label]):
        data = pddata[pddata[label] == lab]
        jsondata = json.loads(data[[x_title, y_title, z_title]].to_json(orient='values'))
        scatter3D.add(str(lab),
                      jsondata,
                      xaxis3d_name='Component 1',
                      yaxis3d_name='Component 2',
                      zaxis3d_name='Component 3', )
    scatter3D.render(output)


def plot3D(data, x_title, y_title, z_title, output):
    '''

    :param data:
    :param x_title:
    :param y_title:
    :param z_title:
    :param output:
    :return:
    '''
    fig = plt.figure(figsize=(50, 50))
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(data[x_title], data[y_title], data[z_title])

    ax.set_xlim(0, 100)
    ax.set_ylim(0, 1000)
    ax.set_zlim(0, max(data[z_title]))
    ax.set_xlabel(x_title)
    ax.set_ylabel(y_title)
    ax.set_zlabel(z_title)
    fig.savefig(output, dpi=300)


def main():
    import pandas as pd
    pddata = pd.read_csv('/wes/chenyl/projects/multiplex_PCR_assess/demo_data/14.csv', header=0)
    # print(pddata[['gc','length']])
    # print(pddata.columns)
    echart3D(pddata, 'gc', 'length', 'count', '/wes/chenyl/projects/multiplex_PCR_assess/demo_data/14.html')


if __name__ == '__main__':
    main()
