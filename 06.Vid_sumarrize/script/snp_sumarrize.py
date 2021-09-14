import glob
import os
import subprocess
import pandas as pd
from collections import Counter


class Sumarrize_Variant:
    """
    Create summarize report of vid from ref or variant.
    Input: result dir of snp 
    Output: ./results/snp_summarize_capture_vid.tsv
    """
    def __init__(self,data_dir,results_dir,detect_share,manual_sample_list,sample_list_dir):
        if manual_sample_list == False:
            cmd = f'ls -l {data_dir} | grep "^d"  | rev | cut -d " " -f 1 | rev &> sample_list.txt'
            subprocess.check_call(cmd,shell=True)
            with open("sample_list.txt") as fd:
                self.sample_names = [sample_names.strip() for sample_names in fd.readlines()]
        else:
            with open(f"{sample_list_dir}") as fd:
                self.sample_names = [sample_names.strip() for sample_names in fd.readlines()]
        if len(set(self.sample_names)) != len(self.sample_names):
            dup = dict(Counter(self.sample_names))
            dup_lst = [dup_name for dup_name,value in dup.items() if value > 1]
            #dup_lst = " ".join(dup_lst)
            print(f"\n WARNING: Your sample name list have duplicate name, Please Check!\n Duplicate name is {dup_lst}\n")

        self.data_dir = data_dir
        self.results_dir = results_dir
        self.detect_share = detect_share
        self.log = open("./log.txt","w")
    
    def get_filter_variant_data(self,sample_name):
        
        df_filter,*other = glob.glob(f"{self.data_dir}/{sample_name}/07.variant_calling/*_filter_variant_count.tsv")
        return df_filter
    
    def get_protein_value(self,sample_name):
        """
        VID Protein
        1   aaa,bbb,cc
        2   eee,ddd
        3   no_detect
        ......
        """
        protein_value,*other = glob.glob(f"{self.data_dir}/{sample_name}/08.analysis_snp/*_variant_table.tsv")
        df_protein_value = pd.read_table(protein_value).loc[:,["VID","Protein"]].fillna("no_detect")

        return df_protein_value

    def get_vid_summarize(self):
        if os.path.exists(self.results_dir):
            pass
        else:
            cmd = f'mkdir {self.results_dir}'
            subprocess.check_call(cmd,shell=True)

        for sample in self.sample_names:
            try:
                filter_variant_count = pd.read_table(self.get_filter_variant_data(sample_name=sample))
                protein_value = self.get_protein_value(sample_name=sample)

                filter_variant_count.loc[:,"vid_judge"] = filter_variant_count.loc[:,"ref_count"] + filter_variant_count.loc[:,"alt_count"]
                filter_variant_count = filter_variant_count[filter_variant_count.loc[:,"vid_judge"] > 0]

                #summarize
                vid_summarize = {}
                vid_summarize["VID"] = list(set(filter_variant_count.loc[:,"VID"]))
                vid_summarize["nCell_with_read_count"] = list(filter_variant_count.groupby("VID")["vid_judge"].count())

                variant_count =  (filter_variant_count.loc[:,"alt_count"] != 0).astype(int)
                ref_count =  (filter_variant_count.loc[:,"ref_count"] != 0).astype(int)
            
                variant_count["VID"] = filter_variant_count.loc[:,"VID"]
                ref_count["VID"] = filter_variant_count.loc[:,"VID"]

                vid_summarize["with_ref_read"] = list(ref_count.groupby("VID").sum())
                vid_summarize["with_variant_read"] = list(variant_count.groupby("VID").sum())

                vid_summarize = pd.DataFrame(vid_summarize)
                
                

                if (self.detect_share == True):
                    vid_summarize.loc[:,"both_ref_and_variant"] =  (vid_summarize.loc[:,"with_ref_read"] + vid_summarize.loc[:,"with_variant_read"]) - vid_summarize.loc[:,"nCell_with_read_count"] 
                    vid_summarize.loc[:,"with_ref_read"] = vid_summarize.loc[:,"with_ref_read"] - vid_summarize.loc[:,"both_ref_and_variant"]
                    vid_summarize.loc[:,"with_variant_read"] = vid_summarize.loc[:,"with_variant_read"] - vid_summarize.loc[:,"both_ref_and_variant"]
                    vid_summarize = pd.merge(left=vid_summarize,right=protein_value,on = "VID")

                    vid_summarize.to_csv(f'{self.results_dir}/{sample}_summarize_capture_vid.tsv', sep='\t', index=False)
                    self.log.write(f"{sample}\tSuccess\n")
                else:
                    vid_summarize = pd.merge(left=vid_summarize,right=protein_value,on = "VID")
                    vid_summarize.to_csv(f'{self.results_dir}/{sample}_summarize_capture_vid.tsv', sep='\t', index=False)
                    self.log.write(f"{sample}\tSuccess\n")
            except:
                    self.log.write(f"{sample}\tFailure\n")


