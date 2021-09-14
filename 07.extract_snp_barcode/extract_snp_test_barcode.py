import pandas as pd
import argparse

class Extract_Data:
    def __init__(self,count_file,variant_table):
        #tsne
        count = pd.read_table(count_file)
        #variant
        variant = pd.read_table(variant_table).fillna("no_pro")
        #barcode和cid的对应关系
        self.cid2barcode = count.set_index(['CID'])['barcode'].to_dict()
        #每个barcode的突变数量
        self.barcode_vcount = count.set_index(['barcode'])['alt_count'].to_dict()
        
        #筛选包含T790M和L858R的突变的细胞
        gene_T790M = variant[variant['Protein'].str.contains("T790M")].loc[:,"CID"]
        gene_L858R = variant[variant['Protein'].str.contains("L858R")].loc[:,"CID"]
        
        #得到cid列表
        gene_T790M_cid = set(list(gene_T790M)[0].strip("(").strip(")").split(","))
        gene_L858R_cid = set(list(gene_L858R)[0].strip("(").strip(")").split(","))
        
        #共有T790M和L858R的barcode(cid)
        self.cid = gene_T790M_cid.intersection(gene_L858R_cid)
        
        
        self.barcode = {}
        for i in self.cid:
            #cid转换为barcode
            key = self.cid2barcode[eval(eval(i))]
            #得到barcao的突变数量
            value = self.barcode_vcount[key]
            self.barcode[key] = value

    def get_brcode_test(self):
        #将count数排序
        count_lst = sorted(list(self.barcode.values()),reverse = True)
        #挑选计数位点最多的barcode前五名
        barcode_test = [bar for bar,value in self.barcode.items() if value in count_lst[0:5]]
        if len(barcode_test) > 5:
            barcode_test = barcode_test[0:5]
        
        with open("./barcode_file.txt","w") as fd:
            for line in barcode_test:
                fd.write(f"{line}\n")
        #return barcode_test

test = Extract_Data(count_file="./data/D0817_PHA6567_SR_11LGR_TS_count_tsne.tsv",
                    variant_table="./data/D0817_PHA6567_SR_11LGR_TS_variant_table.tsv")
test.get_brcode_test()

"""
def main():
    parser = argparse.ArgumentParser(description='extract barcode to make test data')
    parser.add_argument("--count_file", help="count_data", required=True)
    parser.add_argument("--variant_table", help="variant_table", required=True)
    #parser.add_argument("--cell_number",help = "number of cells wants to get")

    args = parser.parse_args()

    test = Extract_Data(count_file=str(args.count_file),variant_table=str(args.variant_table))
    
    test.get_brcode_test()

if __name__ == "__mian__":
    main()

"""




