import pysam
import numpy as np
import pandas as pd
import subprocess
import glob
import os
import time
from collections import defaultdict
from collections import Counter


class Static_Bam:
    def __init__(self,data_dir,results_dir,manual_sample_list,sample_list_dir,bp_min):
        #make sample list
        if manual_sample_list == False:
            cmd = f'ls -l {data_dir} | grep "^d"  | rev | cut -d " " -f 1 | rev &> sample_list.txt'
            subprocess.check_call(cmd,shell=True)
            with open("sample_list.txt") as fd:
                self.sample_names = [sample_names.strip() for sample_names in fd.readlines()]
        else:
            with open(f"{sample_list_dir}") as fd:
                self.sample_names = [sample_names.strip() for sample_names in fd.readlines()]
        
        #check if sample list have duplicate name
        if len(set(self.sample_names)) != len(self.sample_names):
            dup = dict(Counter(self.sample_names))
            dup_lst = [dup_name for dup_name,value in dup.items() if value > 1]
            dup_lst = " ".join(dup_lst)
            print(f"\n WARNING: \n Your sample name list have duplicate name, Please Check!\n This script will drop the duplicate name auto!\n Duplicate name: {dup_lst}\n ")
            self.sample_names = set(self.sample_names)
        
        self.result_dir = results_dir
        if os.path.exists(self.result_dir):
            pass
        else:
            cmd = f'mkdir {self.result_dir}'
            subprocess.check_call(cmd,shell=True)
        
        self.data_dir = data_dir
        self.bp_min = bp_min
        
    def get_effective_position(self,sample_name):
        df_bed_file = pd.read_table(glob.glob(f"{self.data_dir}/{sample_name}/*Aligned.sortedByCoord.out.bed")[0],header=None)

        df_bed_file.loc[:,"len_bp"] = df_bed_file.iloc[:,2] - df_bed_file.iloc[:,1]
        
        bed_file_position_filter = df_bed_file[df_bed_file.loc[:,"len_bp"] >= self.bp_min]
        
        #get start end position
        bed_file_pos_s = list(bed_file_position_filter.iloc[:,1])
        bed_file_pos_e = list(bed_file_position_filter.iloc[:,2])
        bed_file_chr = list(bed_file_position_filter.iloc[:,0])
        
        bed_file_position_message = (bed_file_chr,bed_file_pos_s,bed_file_pos_e)
        df_bed_file_position_name = df_bed_file.iloc[:,[1,3]]
        
        #get position name file
        df_bed_file_position_name.columns = ["position","name"]

        return (bed_file_position_message,df_bed_file_position_name)

    def get_cluster(self,sample_name):
        df_cluster_file = pd.read_table((glob.glob(f"{self.data_dir}/{sample_name}/*tsne_coord.tsv")[0]))
        df_cluster_file.rename(columns={'Unnamed: 0':'barcode'}, inplace = True)
        df_barcode_cluster = df_cluster_file.loc[:,["barcode","cluster"]]
        return df_barcode_cluster
    
    def gte_samfile(self,sample_name):
        bam_file = glob.glob(f"{self.data_dir}/{sample_name}/*filtered_sorted.bam")[0]
        #if index have
        if len(glob.glob(f"{self.data_dir}/{sample_name}/*.bai")) == 0:
            print(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}\t - Make index of {sample_name} bam file...")
            
            cmd = f"samtools index {bam_file}"
            subprocess.check_call(cmd,shell=True)
            
            bam_file_index = glob.glob(f"{self.data_dir}/{sample_name}/*bam.bai")[0]
            samfile = pysam.AlignmentFile(bam_file,"rb",index_filename = bam_file_index)

        else:
            bam_file_index = glob.glob(f"{self.data_dir}/{sample_name}/*.bai")[0]
            samfile = pysam.AlignmentFile(bam_file,"rb",index_filename = bam_file_index)
        return samfile
    
    def get_cover_condition(self):
        for sample in self.sample_names:
            print(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}\t - Calculate the {sample}....")
            samfile_sample = self.gte_samfile(sample_name = sample)
            position_lst = self.get_effective_position(sample_name = sample)
            position_name = position_lst[1]
            cluster = self.get_cluster(sample_name = sample)

            #make a  defaultdict to store data
            position_barcode = defaultdict(list)
            position_star_end = position_lst[0]

            
            #pysam fetch
            for chr,start,end in zip(position_star_end[0],position_star_end[1],position_star_end[2]):
                for readr in samfile_sample.fetch(str(chr),int(start),int(end)):
                    length_overlap = readr.get_overlap(int(start),int(end))
                    #print(readr.get_refernece_name())
                    #length_overlap = readr.reference_length
                    #length_overlap  = readr.query_alignment_length
                    #length_overlap = readr.infer_query_length()
                    #length_overlap = readr.infer_read_length()
                    """
                    https://pysam.readthedocs.io/_/downloads/en/latest/pdf/ 
                    page 14
                    get_overlap(self, uint32_t start, uint32_t end)
                    return number of aligned bases of read overlapping the interval start and end on the reference sequence.
                    Return None if cigar alignment is not available.
                    """
                    if length_overlap >= self.bp_min:
                        #add position message
                        position_barcode["chr"].append(str(chr))
                        position_barcode["start"].append(int(start))
                        position_barcode["end"].append(int(end))
                        position_barcode["len_overlap"].append(int(length_overlap))
                        try:
                            barcode = readr.get_tag("CB")
                            position_barcode["barcode"].append(barcode)
                        except:
                            position_barcode["barcode"].append("no_CB")
                        try:
                            position_barcode["UMI"].append(readr.get_tag("UB"))
                        except:
                            position_barcode["UMI"].append("no_UB")
                        try:
                            r_name = readr.query_name.split('_')[:2]
                            r_name = "_".join(r_name)
                            position_barcode["reads_name"].append(r_name)
                        except:
                            position_barcode["reads_name"].append("no_name")
            
            #position_barcode barcode umi dataframe
            df_pos = pd.DataFrame(position_barcode)
            #filter
            df_pos = df_pos[df_pos.loc[:,"barcode"] != "no_CB"]
            df_pos = df_pos[df_pos.loc[:,"UMI"] != "no_UB"]
            df_pos = df_pos[df_pos.loc[:,"reads_name"] != "no_name"]
            
            #save message
            umi_barcode_mesage = defaultdict(list)
            #add cell
            pos_start = position_star_end[1]
            pos_end = position_star_end[2]
            chr_s = position_star_end[0]
            #print(df_pos)
            
            for chr,start,end in zip(chr_s,pos_start,pos_end):
                df_pos_tmp = df_pos[(df_pos.loc[:,"start"] == start) & (df_pos.loc[:,"end"] == end) & (df_pos.loc[:,"chr"] == str(chr))]
                #add position message
                umi_barcode_mesage["chr"].append(str(chr))
                umi_barcode_mesage["start"].append(start)
                umi_barcode_mesage["end"].append(end)
                
                #add cell
                cell = Counter(list(df_pos_tmp.loc[:,"barcode"]))
                umi_barcode_mesage["cell"].append(len(cell.keys()))
                
                #add reads
                reads = Counter(list(df_pos_tmp.loc[:,"reads_name"]))
                umi_barcode_mesage["reads"].append(sum(reads.values()))
                umi_barcode_mesage["medium_reads"].append(np.median(list(cell.values())))
                
                #add umi
                umi = Counter(list(df_pos_tmp.loc[:,"UMI"]))
                umi_barcode_mesage["UMI"].append(len(umi.keys()))
                umi_barcode_mesage["medium_UMI"].append(np.median(list(umi.values())))
          
            df_umi_barcode_mesage = pd.DataFrame(umi_barcode_mesage)
            df_umi_barcode_mesage.loc[:,"mean_reads"] = df_umi_barcode_mesage.loc[:,"reads"] / df_umi_barcode_mesage.loc[:,"cell"]
            df_umi_barcode_mesage.loc[:,"mean_UMI"] = df_umi_barcode_mesage.loc[:,"UMI"] / df_umi_barcode_mesage.loc[:,"cell"]

            df_umi_barcode_mesage = pd.merge(left=df_umi_barcode_mesage,
                                             right = position_name,
                                             left_on = "start",
                                             right_on ="position",
                                             how = "left").drop("position",axis = 1).fillna(0)
            cluster_message = pd.merge(left = df_pos,
                                       right = cluster,
                                       on = "barcode",
                                        how = "left")
            cluster_message = pd.merge(left = cluster_message,
                                       right = position_name,
                                       left_on = "start",
                                       right_on="position",
                                       how = "left").drop(["UMI","reads_name","position","len_overlap"],axis = 1)
            cluster_message.to_csv(f"{self.result_dir}/{sample}_cluster_message.tsv",index=None,sep = "\t")
            df_umi_barcode_mesage.to_csv(f"{self.result_dir}/{sample}_umi_barcode_mesage.tsv",index=None,sep = "\t")





            

                




