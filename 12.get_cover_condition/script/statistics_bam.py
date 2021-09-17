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
    def __init__(self,data_dir,results_dir,manual_sample_list,sample_list_dir,bp_min,cosmic_database):
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
        
        self.cosmic_database = cosmic_database
        self.data_dir = data_dir
        self.bp_min = bp_min
        self.log = open("log.txt","w")
        
    def get_effective_position(self,sample_name):
        """
        two return values
        1.
        [pos1,pos2,pos3....]
        default: bp>=50
        2.
        pos     name
        pos1    name1
        pos2    name2
        pos3    name3
        ......
        """
        df_bed_file = pd.read_table(glob.glob(f"{self.data_dir}/{sample_name}/*Aligned.sortedByCoord.out.bed")[0],header=None)
        df_bed_file.loc[:,"len_bp"] = df_bed_file.iloc[:,2] - df_bed_file.iloc[:,1]
        bed_file_position_filter = list(df_bed_file[df_bed_file.loc[:,"len_bp"] >= self.bp_min].iloc[:,1])
        df_bed_file_position_name = df_bed_file.iloc[:,[1,3]]
        df_bed_file_position_name.columns = ["position","name"]

        return (bed_file_position_filter,df_bed_file_position_name)
    
    def get_cosmic_data(self):
        #read the database and get the position of database
        print(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}\t - Loading Database... (location: {self.cosmic_database})")
        database = pd.read_table(self.cosmic_database,low_memory=False,header=None)
        pos_database = set(database.iloc[:,1])
        
        try:
            for sample in self.sample_names:
                print(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}\t - Find the intersection with bed file.... . Sample: {sample}")
                pos_effective = set(self.get_effective_position(sample_name = sample)[0])
                pos_intersection = pos_effective.intersection(pos_database)
                if len(pos_intersection) == 0:
                    print(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}\t - No detect intersection with bed file.... . Sample: {sample}")
                    self.log.write(f"{sample}\t No detect intersection\n")
                else:
                #get data 
                    df_pos_intersection = database[database.iloc[:,1].isin(pos_intersection)]
                    df_pos_intersection.to_csv(f"{self.result_dir}/{sample}_intersection_with_cosmic.tsv",sep = "\t",index=None)
                    self.log.write(f"{sample}\t Find intersection sucess\n")
        except:
            self.log.write(f"{sample}\t Find intersection sucess Failture\n")

    def get_cluster(self,sample_name):
        """
        barcode cluster
        bar1    1
        bar2    3
        bar3    5
        ......
        """
        df_cluster_file = pd.read_table((glob.glob(f"{self.data_dir}/{sample_name}/*tsne_coord.tsv")[0]))
        df_cluster_file.rename(columns={'Unnamed: 0':'barcode'}, inplace = True)
        df_barcode_cluster = df_cluster_file.loc[:,["barcode","cluster"]]
        return df_barcode_cluster
    
    def gte_samfile(self,sample_name):
        """
        get the iter of samfile 
        """
        bam_file = glob.glob(f"{self.data_dir}/{sample_name}/*filtered_sorted.bam")[0]
        #if index have
        if len(glob.glob(f"{self.data_dir}/{sample_name}/*.bai")) == 0:
            print(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}\t - No detect the index of {sample_name} bam file")
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
        """
        position barcode cluster    name
        pos1     bar1    1          name1
        pos2     bar2    2          name2
        ......
        """
        try:
            for sample in self.sample_names:
                print(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}\t - Calculate the {sample}....")
                samfile_sample = self.gte_samfile(sample_name = sample)
                position_lst = self.get_effective_position(sample_name = sample)
                cluster = self.get_cluster(sample_name = sample)
                #make a  defaultdict to store data
                position_barcode = defaultdict(list)
                for reader in samfile_sample:
                    try:
                        barcode = reader.get_tag('CB')
                        position_barcode["barcode"].append(barcode)
                    except:
                        position_barcode["UMI"].append("no_barcode")
                    try:             
                        umi = reader.get_tag('UB')
                        position_barcode["UMI"].append(umi)
                    except:
                        position_barcode["UMI"].append("no_umi")
                
                    try:
                        gene_name = reader.get_tag('GN')
                        position_barcode["gene"].append(gene_name)
                    except:
                        position_barcode["gene"].append("no_gene_name")  
                
                    position_barcode["gene_len"].append(reader.query_length)
                    position_barcode["position"].append(reader.reference_start)
                #to dataframe  
                position_barcode = pd.DataFrame(position_barcode)
            
                 #filter by gene_len
                position_barcode = position_barcode[position_barcode.loc[:,"gene_len"] >= self.bp_min]

                #filter by position
                position_barcode = position_barcode[position_barcode.loc[:,"position"].isin(position_lst[0])]
                #merge cluster and name
                position_barcode_cluster = pd.merge(left=position_barcode,
                                                    right=cluster,
                                                    on = "barcode",
                                                    how = "left").fillna("no_cluster")
                position_barcode_cluster_name = pd.merge(left=position_barcode_cluster,
                                                    right=position_lst[1],
                                                    on = "position",
                                                    how = "left").fillna("no_name")
                #save data of position_barcode_cluster_name
                data_save_colum = ["name","position","barcode","cluster"]
                position_barcode_cluster_name.loc[:,data_save_colum].to_csv(f"{self.result_dir}/{sample}_position_barcode_cluster_name.tsv",sep = "\t",index=None)
            
           
                #create a collection object to save the umi reads message
                read_umi_message = defaultdict(list)
            
                 #get position list
                pos_lst = set(position_barcode.loc[:,"position"])

                for position in pos_lst:
                    df_pos =  position_barcode[position_barcode.loc[:,"position"] == position]
                    #add position
                    read_umi_message["position"].append(position)

                     #add cell
                    cell = Counter(list(df_pos.loc[:,"barcode"]))
                    read_umi_message["cell"].append(len(cell.keys()))
                
                    #add medium_reads
                    medium_reads = np.median(list(cell.values()))
                    read_umi_message["medium_reads"].append(medium_reads)

                    #add UMI message
                    umi = Counter(list(df_pos.loc[:,"UMI"]))
                    read_umi_message["umi_number"].append(len(umi.keys()))

                    #add reads_number
                    read_umi_message["reads_number"].append(sum(umi.values()))

                    #add medium_umi
                    medium_umi = np.median(list(umi.values()))
                    #medium_umi = np.median(list(df_pos.groupby("barcode")["UMI"].value_counts()))
                
                    read_umi_message["medium_umi"].append(medium_umi)
            
                df2_umi_reads_message = pd.DataFrame(read_umi_message)
                #add mean_reads
                df2_umi_reads_message.loc[:,"mean_reads"] = df2_umi_reads_message.loc[:,"reads_number"] / df2_umi_reads_message.loc[:,"cell"]
            
                #add mean_umis
                df2_umi_reads_message.loc[:,"mean_umi"] = df2_umi_reads_message.loc[:,"umi_number"] / df2_umi_reads_message.loc[:,"cell"]
            
                #add name
                df2_umi_reads_message = pd.merge(left = df2_umi_reads_message,
                                                    right = position_lst[1],
                                                    on = "position",
                                                    how = "left")

                #save data
                sort_index = ["position","name","cell","reads_number","umi_number","mean_reads","mean_umi","medium_reads","medium_umi"]
                df2_umi_reads_message.loc[:,sort_index].to_csv(f"{self.result_dir}/{sample}_umi_reads_message.tsv",sep = "\t",index=None)
                self.log.write(f"{sample}\t Sucess\n")
        except:
            self.log.write(f"{sample}\t Failture\n")
