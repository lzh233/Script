# 应用场景

按不同的分类层级对out表进行注释并计算相对丰度和counts

# 方法

## 输入文件

otutable

```R
#OTU_ID	DD10_2	DD20_1	DD20_6	DD20_3	DD20_4	DD20_5	DD10_4	DD10_3	DD10_1	DD20_2	DD10_5	DD10_6
OTU4342881	1	1	0	1	0	0	1	2	3	0	1	2
OTU215575	3	1	4	2	4	3	1	1	5	1	1	2
OTU222196	12	4	1	8	9	3	6	6	9	6	9	7
```

注释文件

```R
OTUOTUID	Kingdom	Phylum	Class	Order	Family	Genus	Species
OTU187144	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__	 g__	 s__
OTU836974	k__Bacteria	 p__Cyanobacteria	 c__Chloroplast	 o__Cercozoa	 f__	 g__	 s__
OTU310669	k__Bacteria	 p__Firmicutes	 c__Clostridia	 o__Clostridiales	 f__	 g__	 s__

```

## 使用

```r
source(".//script/NoteAndExtract.R")
#all: 按照所有分类水平分别提取
#extract.by： 可指定分类水平（只能指定一个） 如，Phylum
#caclulate.ra是否计算丰度
#verbose是否显示输出信息
#otutable otu表 和注释表的out名一致
#note.table 注释表 和otu表的out名一致
NoteAndExtract(otutable = otu,
               note.table = note,
               extract.by = "all",
               caclulate.ra = T,
               verbose = T)
#指定分类水平提取
NoteAndExtractSelect(otutable = otu,
                     note.table = note,
                     extract.select = c("Phylum","Family","Class"),
                     calculate.ra = T,
                     verbose = T)
#批量提取extract.by可指定分类水平（只能指定一个） 如，Phylum
NoteAndExtractAll(otu.dir = "./data2/",
                  note = note,
                  extract.by = "all",
                  calculate.ra = T,
                  sep.otu = ",",
                  verbose = T)

```

