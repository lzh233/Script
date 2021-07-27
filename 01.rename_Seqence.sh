#!/bin/bash
#Date:20201020
#Author:lzh
#定义行号变量，line_num 
#当前序列替换条数，seq_count
#当前序列替换条数，seq_count
#总条数，seq
#序列名，seq_name
#输入文件：seq.fasta、name.txt
#输出文件：seq_name.fasta

line_num=1
seq_count=0
cp -a ./seq.fasta ./seq_name.fasta
cp -a ./name.txt ./name.txt.tmp
seq=$(wc -l seq_name.fasta | awk '{print $1}')

#名称开头加大于号
function format(){
	sed -i 's/^/\>/g' name.txt.tmp
}

#批量命名，按照行号1 3 5 7 ....插入
function rename(){
	for seq_name in $(cat name.txt.tmp)
		do
			sed -i "${line_num} i ${seq_name}" seq_name.fasta
			line_num=$(($line_num+2))
			seq_count=$(($seq_count+1))
			echo -ne "共${seq}条序列，已经命名${seq_count}条\r"
		done
}
#若名称提前加入“>”,在语句format前加#---> #format
#format
rename
rm -rf name.txt.tmp