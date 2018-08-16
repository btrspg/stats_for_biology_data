## 例子

[https://btrspg.github.io/stats_for_biology_data/](https://btrspg.github.io/stats_for_biology_data/)



## 目前需要改进的地方

- [x] 计算bam中每个contig的gc+length count分布的时候，耗费的时间特别长。没有办法对基因组上的一个contig进行统计。
- [ ] 能不能找到一个对列表中的元素同时进行处理的方法（类似与numpy的那种）
- [X] gc+length+count图（一般3D图）
- [X] 如果能加上echart图会更好一些，因为需要展示一些效果
- [X] 在[https://btrspg.github.io/stats_for_biology_data/](https://btrspg.github.io/stats_for_biology_data/)中可以看到使用matplotlib和pyecharts画的3D图（同一套数据），发现3D的图还是pyecharts更加适合一些。
- [ ] 需要加入深度计算的方法
- [ ] 校正方法及参数还没有确定，但是计算或者统计方法可以先写好，后续直接补入算法即可

## 改进方法

- [X] 采取抽样:这个也写了，目前计算的可接受范围内，我认为是每个contig大约100W个模版的计算。
- [X] 改进算法:这个是因为之前的程序有bug所以导致计算时间非常长而且有错误，现阶段已经在可接受的范围内了。
- [ ] 按照一个contig来进行统计，这样存储的时候内存压力不会那么大
- [ ] 每一个contig每一个位置的深度存储为临时文件，最终整合后在输出。
- [ ] 是否能够区分pcr-free和一般的wgs的数据，可以试试使用PCA或者T-SNE来进行一下图形化看看，刚好可以学习一下T-SNE及PCA的使用方法。

## 问题

- [ ] 现在有一个问题，让我一直想不明白，假设有一个均匀分布的样本，由于某些条件的加入，导致他服从了正态分布，我要怎么做才能再将它变换回均匀分布。
- [ ] GC含量的峰值是不是由于基因组本身GC含量所造成的。
