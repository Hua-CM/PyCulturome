# PyCulturome
![License](https://img.shields.io/badge/GPL-3.0-brightgreen) 
![Python](https://img.shields.io/badge/Python-3.8.0-brightgreen)

一个用于微生物组中特殊细菌的高通量培养和鉴定的一键式Python工作流（Culturome）

## 安装
~~~shell
git clone 
cd PyCulturome
pip install -r requirements.txt
python pipe.py
~~~

## 用法

运行示例
~~~shell
python /path/to/PyCulturome/pipe.py \
 -1 /path/to/your/fastq1 \
 -2 /path/to/your/fastq2 \
 -f example/input/fwd_bar.tsv \
 -r example/input/rev_bar.tsv \
 -d /home/database/micro/SILVA138_RESCRIPt.fasta \
 -o example/output
 ~~~
 
**-1**：Read1 输入文件路径
**-2**：Read2 输入文件路径
**-f**：正向条形码信息，以制表符分隔格式，包含两列：孔 id \t 序列
**-r**：反向条形码信息，以制表符分隔格式，包含两列：板 id \t 序列
**-d**：用于 vsearch 序列注释的分类数据库
**-o**：输出目录

> 目前仅支持799F和1193R，正向10bp，反向6bp的测序结果拆分。

## 输出
### 序列
**\<名称>.rep.fa**: 包含ASV代表序列的fasta文件。

### 表格
**\<名称>.abun.tsv**: 宽表格格式的OTU丰度表。由usearch直接生成。

**readscount_well.tsv**: 长表格格式的OTU丰度表。

**purified_well.tsv**: 纯化的孔位信息。

**recommend_well4asv.tsv**: 建议分离每个ASV的孔位数量。

**taxonomy_8.tsv**: 每个ASV的分类信息。

### 图像
**Positive_Negative_Control.pdf**: 负控制和正控制孔中读取数目的箱线图。

**Reads_frequency.pdf**: 每个孔读取数目的直方图。

**Purity_frequency.pdf**: 孔位纯度（从0到100%）的直方图。

**Rarefication.pdf**: 随着孔位数量变化的ASV稀释曲线。

**Tree.pdf**: 分离菌株的分类树。

### 依赖
- [usearch](https://github.com/rcedgar/usearch_old_binaries) > v11
- [vsearch](https://github.com/torognes/vsearch) >= v2.22
- python >= 3.8
- Biopython == 1.80
    - pandas > 2.0.0
    - numpy > 1.23.0
    - ete3 > 3.0.0
    - seaborn > 0.13.0
    - matplotlib > 3.6.0

## TODO
1. 比较和/或合并多个菌库的结果
2. 增加参数支持其他引物和barcode

## 参考
[Culturome-YongxinLiu](https://github.com/YongxinLiu/Culturome/tree/master)