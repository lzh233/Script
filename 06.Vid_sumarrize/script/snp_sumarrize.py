import glob
import os
import subprocess
import pandas as pd
from collections import Counter
from collections import defaultdict 


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
        
        #check if sample list have duplicate name
        if len(set(self.sample_names)) != len(self.sample_names):
            dup = dict(Counter(self.sample_names))
            dup_lst = [dup_name for dup_name,value in dup.items() if value > 1]
            dup_lst = " ".join(dup_lst)
            print(f"\nWARNING: Your sample name list have duplicate name, Please Check!\nDuplicate name is {dup_lst}\n")

        self.data_dir = data_dir
        self.results_dir = results_dir
        self.detect_share = detect_share
        self.log = open("./log.txt","w")
    
    def get_cid_cluster(self,sample_name):
        """
        CID cluster
        1   1
        2   5
        3   4
        ......
        This cluster message is from the celescop rna results !!
        """
        barcode_cid_file = pd.read_table(glob.glob(f"{self.data_dir}/{sample_name}/*_tsne_coord.tsv")[0])
        barcode_cid_file.columns = ["barcode","tSNE_1","tSNE_2" ,"cluster" ,"Gene_Counts"]
        barcode_cid_file = barcode_cid_file.loc[:,["barcode","cluster"]]

        #read CID file
        cid_barcode_file = pd.read_table(glob.glob(f"{self.data_dir}/{sample_name}/07.variant_calling/*_CID.tsv")[0]).loc[:,["CID","barcode"]]
        
        cid2cluster = pd.merge(left = barcode_cid_file,
                               right = cid_barcode_file,
                               on = "barcode",
                               how = "left")
        cid2cluster = cid2cluster.loc[:,["CID","cluster"]] 
        return cid2cluster
    
    def get_filter_variant_data(self,sample_name):
        df_filter = glob.glob(f"{self.data_dir}/{sample_name}/07.variant_calling/*_filter_variant_count.tsv")[0]
        return df_filter
    
    def get_protein_value(self,sample_name):
        """
        VID Protein
        1   aaa,bbb,cc
        2   eee,ddd
        3   no_detect
        ......
        """
        protein_value = glob.glob(f"{self.data_dir}/{sample_name}/08.analysis_snp/*_variant_table.tsv")[0]
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
                cid_cluster = self.get_cid_cluster(sample_name=sample)
                vid_summarize = {}
                
                #drop the vid_judge = 0
                filter_variant_count.loc[:,"vid_judge"] = filter_variant_count.loc[:,"ref_count"] + filter_variant_count.loc[:,"alt_count"]
                filter_variant_count = filter_variant_count[filter_variant_count.loc[:,"vid_judge"] > 0]

                #add data to vid_summarize
                vid_summarize["VID"] = list(set(filter_variant_count.loc[:,"VID"]))
                vid_summarize["nCell_with_read_count"] = list(filter_variant_count.groupby("VID")["vid_judge"].count())
                #add the cell number of alt_count and ref_count
                variant_count =  (filter_variant_count.loc[:,"alt_count"] != 0).astype(int)
                ref_count =  (filter_variant_count.loc[:,"ref_count"] != 0).astype(int)
                
                #add VID colums to merge data
                variant_count["VID"] = filter_variant_count.loc[:,"VID"]
                ref_count["VID"] = filter_variant_count.loc[:,"VID"]
                
                #add number of alt_count and ref_count to vid_summarize
                vid_summarize["with_ref_read"] = list(ref_count.groupby("VID").sum())
                vid_summarize["with_variant_read"] = list(variant_count.groupby("VID").sum())
                vid_summarize = pd.DataFrame(vid_summarize)

                #add cluster column to filter_variant_count
                variant_with_cluster = pd.merge(left = filter_variant_count.loc[:,["VID","CID"]],
                                                right=cid_cluster, 
                                                on = "CID",
                                                how = "left").fillna("no_cluster")
                #get vid and cluster 
                vid_set = set(variant_with_cluster.loc[:,"VID"])
                cluster_set = set(variant_with_cluster.loc[:,"cluster"])
                
                #set function class the cluster
                class_cluster = lambda cluster_count_list,cluster :cluster_count_list[cluster]
                
                #set a defaultdict to store results
                cluster_result = defaultdict(list)

                #get result
                for vid in vid_set:
                    cluster_count = Counter(list(variant_with_cluster[variant_with_cluster.loc[:,"VID"] == vid].loc[:,"cluster"]))
                    for cluster in cluster_set:
                        if cluster == "no_cluster":
                            cluster_result[f"no_cluster"].append(class_cluster(cluster_count_list = cluster_count,cluster = cluster))
                        else:
                            cluster_result[f"cluster_{int(cluster)}"].append(class_cluster(cluster_count_list = cluster_count,cluster = cluster))
                    
                cluster_result = pd.DataFrame(cluster_result)
                
            
                if (self.detect_share == True):
                    #get the both_ref_and_variant
                    vid_summarize.loc[:,"both_ref_and_variant"] =  (vid_summarize.loc[:,"with_ref_read"] + vid_summarize.loc[:,"with_variant_read"]) - vid_summarize.loc[:,"nCell_with_read_count"] 
                    vid_summarize.loc[:,"with_ref_read"] = vid_summarize.loc[:,"with_ref_read"] - vid_summarize.loc[:,"both_ref_and_variant"]
                    vid_summarize.loc[:,"with_variant_read"] = vid_summarize.loc[:,"with_variant_read"] - vid_summarize.loc[:,"both_ref_and_variant"]
                    #add protrie data
                    vid_summarize = pd.merge(left=vid_summarize,right=protein_value,on = "VID")
                    #add cluster data
                    vid_summarize = pd.concat([vid_summarize,cluster_result],axis = 1)
                    #save vid_summarize
                    vid_summarize.to_csv(f'{self.results_dir}/{sample}_summarize_capture_vid.tsv', sep='\t', index=False)
                    #add log
                    self.log.write(f"{sample}\tSuccess\n")
                else:
                    vid_summarize = pd.merge(left=vid_summarize,right=protein_value,on = "VID")
                    vid_summarize = pd.concat([vid_summarize,cluster_result],axis = 1)

                    vid_summarize.to_csv(f'{self.results_dir}/{sample}_summarize_capture_vid.tsv', sep='\t', index=False)
                    self.log.write(f"{sample}\tSuccess\n")
            except:
                self.log.write(f"{sample}\tFailure\n")


