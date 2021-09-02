#!/bin/bash
#Author:wangwe&liuzihao
#Date:20210221
#切换到工作目录

targetTaxa=Escherichia_coli
mkdir -p ~/testgene2/references/${targetTaxa}
cd ~/testgene2/references/${targetTaxa}
count=0
wget ftp://ftp.ncbi.nih.gov/genomes/genbank/bacteria/Escherichia_coli/assembly_summary.txt
#将第12列为complete genome的第20列FTP地址打印出来，存在assembly_summary_complete_genomes.txt

awk -F '\t' '{if($12=="Complete Genome") print $20}' assembly_summary.txt > assembly_summary_genomes.txt

awk '{print $0 "/*"}' assembly_summary_genomes.txt  > assembly_summary_complete_genomes.txt

function create_dir(){
	#网址中提取菌株信息
		bac_name=$(sed -n ${count}p assembly_summary_complete_genomes.txt | cut -d "/" -f 10)
		echo "${bac_name}"
		#创建存放数据目录	
		mkdir -p ${bac_name}
		cd ${bac_name}
}

function get_data(){
	for next in $(cat assembly_summary_complete_genomes.txt); 
		do 
			count=$((count+1))
			create_dir
			wget -P ~/testgene2/references/${targetTaxa}/${bac_name} ${next} >> log.txt
			cd ~/testgene2/references/${targetTaxa}
			echo "下载${count}个"

		done
}
#运行程序
get_data