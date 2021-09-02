from collections import defaultdict
import pysam
import numpy as np

#构建count_dict
def genDict(dim=3, valType=int):
    if dim == 1:
        return defaultdict(valType)
    else:
        return defaultdict(lambda: genDict(dim - 1, valType=valType))

count_dict = genDict(dim=3, valType=int)

samfile = pysam.AlignmentFile("sample1_Aligned.sortedByCoord.out.bam.featureCounts.bam","rb")

#构建一个gene_list和match_list用于模拟

gene_list = []
match_barcode = []

for readr in samfile:
    try:
        gene_name = readr.get_tag('GN')
        gene_list.append(gene_name)
    except KeyError:
        continue
    try:
        barcode = readr.get_tag('CB')
        UMI = readr.get_tag('UB')
        match_barcode.append(barcode)
    except KeyError:
        continue
    count_dict[barcode][gene_name][UMI] += 1

gene_list = gene_list
match_barcode = match_barcode
"""
count_dict结构字典套字典 ({barcode1:{gene_name1:{umi:1,umi2:2,umi3:3....}}},.....)
count_dict[barcode]
defaultdict(<function genDict.<locals>.<lambda> at 0x7f80bcae8f28>, {'TFAM': defaultdict(<class 'int'>, {'TGATGAAT': 2}),.....}) 
"""
total_reads = 0
enriched_reads = 0
enriched_reads_in_cells = 0
enriched_reads_per_cell_list = []

#读取barcode
for barcode in count_dict:
    #每个cell富集的umi数量
    cell_enriched_reads = 0
    #取出每个barcode中包含的基因名
    for gene_name in count_dict[barcode]:
        #每个基因中对应的umi个数 所以用len()
        #gene_UMI = len(count_dict[barcode][gene_name])
        #计算每个reads被富集的数量
        gene_reads = sum(count_dict[barcode][gene_name].values())
        #每次总的umi 累加
        total_reads += gene_reads
        #判断基因是否在需要的基因列表中，如果该基因在，则认为该基因为要的基因 对应的umi个数(len())为富集的umi数量
        if gene_name in gene_list:
            enriched_reads += gene_reads
            #进一步判断这个barcode是不是在barcode列表中，如果在 那么这个为富集在细胞内的umi
            if barcode in match_barcode:
                enriched_reads_in_cells += gene_reads
                cell_enriched_reads += gene_reads

    if barcode in match_barcode:
        enriched_reads_per_cell_list.append(cell_enriched_reads)
