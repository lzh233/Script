# get_cover_condition.py

## environment

```python
python 3.6 or later
require: pysam samtools
```

```apl
Script: /SGRNJ03/randd/user/liuzihao/statistics_bam/
Test data: /SGRNJ03/randd/user/liuzihao/statistics_bam/data
```

## options

`--data_dir`: The directory which saved the data, make sure this directory only save the result which you need, if not, please use  `--manual_sample_list` and `--sample_list_dir` to provide a sample name list.

 `--results_dir` : Output directory.

`--manual_sample_list`: If provide this option, the script will not create samlple name list by auto.

`--sample_list_dir`: When use `--manual_sample_list` you must use  `--sample_list_dir` to provide a sample name list file.

`--bp_min`: Provide a min length to filter data (default: 60)

## data_dir  sample list and COSMIC_DATABASE.conf

### data_dir

**The bam file may not have index file,  the script will run `samtools index <bam file>` by auto to make a bai file**

```shell
$ ls data/
test1  test2  test3
$ ls ./data/test1/
B0831_R1_V2FJ_filtered_sorted.bam  B0831_R1_V2FJ_filtered_sorted.bam.bai  B0831_R1_V2ZL_tsne_coord.tsv  blood_Aligned.sortedByCoord.out.bed
```

### sample  list

**if you provide a manual sample list ,it should not have duplicate name, if have, the script will drop it by auto**

```shell
test1
test2
test3
```

### COSMIC_DATABASE.conf

The script read the  `Database.conf` to get the location of  `Cosmic Database`

```shell
[Database]
db = /SGRNJ/Database/script/database/annovar/humandb/hg38_cosmic70.txt
```

if you want to change the database,  you should change the `Database.conf` file. 

For example

```shell
[Database]
db = /Personal/....
```

## usage

```shell
#Run test data
bash -x run_test.sh
```

```shell
#Quick start
python get_cover_condition.py --data_dir ./data/ --results_dir ./results
```

```shell
# If you want set a min length of reads use option of --bp_min n
python  get_cover_condition.py --data_dir ./data/ --results_dir ./results --bp_min 100
```

```shell
# If you want to manual make sample name list
python get_cover_condition.py --data_dir ./snp_data/ --results_dir ./results --manual_sample_list --sample_list_dir ./sample_list.txt
```

## results

**One sample will out put 3 files**

```python
"""
Format
{sample name}_position_barcode_cluster_name.tsv 
"""
$ head test1_position_barcode_cluster_name.tsv 
name	position	barcode	cluster
WT1_PCR3	32396281	TGGTGGTAAGAGTCAACCGACAAC	6
WT1_PCR3	32396281	AGTACAAGGTGTTCTACACCTTAC	2
WT1_PCR3	32396281	CAGCGTTAAACGCTTACTGAGCCA	2
WT1_PCR3	32396281	ATTGGCTCAACGTGATAGCAGGAA	1
WT1_PCR3	32396281	CATACCAAACAGCAGAGTCGTAGA	1
WT1_PCR3	32396281	ACAGATTCCTGGCATAATCATTCC	2
WT1_PCR3	32396281	AATCCGTCAATGTTGCACGTATCA	1
WT1_PCR3	32396281	AAACATCGCAAGGAGCACAGCAGA	1
WT1_PCR3	32396281	CACTTCGAAGTCACTACCAGTTCA	1
......
```

```python
"""
Format
{sample name}_umi_reads_message.tsv
"""
$ head test1_umi_reads_message.tsv 
position	name	cell	reads_number	umi_number	mean_reads	mean_umi	medium_umi	medium_reads
7675105	TP53-2轮扩增-4	1964	7293	2970	3.7133401221995928	1.5122199592668024	1.0	2.0
7673765	TP-53-扩增-2	3082	20369	5977	6.609020116807268	1.9393251135626217	2.0	4.0
37171397	Blood-PIM1-2轮扩增-3	3618	74403	13326	20.564676616915424	3.683250414593698	3.0	12.0
25227401	Blood-KRAS-2轮扩增-2	5	6	5	1.2	1.0	1.0	1.0
......
```

```python
"""
Format
{sample name}_intersection_with_cosmic.tsv
"""
$ cat test1_intersection_with_cosmic.tsv 
1	69345	69345.1	C	A	ID=COSM911918;OCCURENCE=1(endometrium)
2	25245272	25245281	CCTCCAACGA	-	ID=COSM1583065;OCCURENCE=1(haematopoietic_and_lymphoid_tissue)
7	55181306	55181306	T	C	ID=COSM28943;OCCURENCE=1(lung)
7	55181312	55181312	G	A	ID=COSM12989;OCCURENCE=1(lung)
......
```

## log.txt

```python
#You can check log.txt to konw is it successful
$ cat log.txt 
test1	 Sucess
test2	 Sucess
test3	 Sucess
test1	 Find intersection sucess
test2	 Find intersection sucess
test3	 Find intersection sucess
```
