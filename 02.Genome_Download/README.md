# 指定物种基因组批量下载脚本

## 运行环境

```
Author:wangwe&liuzihao
Date:2021-4-9
CentOS Linux release 7.8.2003 (Core)	
Architecture:          x86_64
CPU op-mode(s):        32-bit, 64-bit
Byte Order:            Little Endian
CPU(s):                144
Model name:            Intel(R) Xeon(R) Gold 5220 CPU @ 2.20GHz
```

## 说明

实例为下载大肠杆菌属的基因组信息，可以通过修改脚本中以下的信息下载其他物种的基因组，默认存放目录为用户家目录

```shell
targetTaxa=Escherichia_coli
ftp://ftp.ncbi.nih.gov/genomes/genbank/bacteria/Escherichia_coli/assembly_summary.txt
```

## 运行

```shell
bash 02.Genome_Download.sh
```

