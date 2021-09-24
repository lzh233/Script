# get_cover_condition.py

## environment

```python
python 3.6 or later
require: pysam samtools
```

## options

`--data_dir`: The directory which saved the data, make sure this directory only save the result which you need, if not, please use  `--manual_sample_list` and `--sample_list_dir` to provide a sample name list.

 `--results_dir` : Output directory.

`--manual_sample_list`: If provide this option, the script will not create samlple name list by auto.

`--sample_list_dir`: When use `--manual_sample_list` you must use  `--sample_list_dir` to provide a sample name list file.

`--bp_min`: Provide a min length to filter data (default: 50)

## data_dir  sample list and Database.conf

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
```

## usage

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
