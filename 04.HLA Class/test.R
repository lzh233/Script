setwd("D:\\Desktop\\GitHub\\Data_private\\Demo_work\\liuzihao")
library(tidyverse)

HLA_data <- read_tsv(".//data/test_input.tsv") 
class <- read_tsv(".//data/hla_known_new.txt")

source(".//RScript/Main.R")
Cell_Class(data = HLA_data,
           col_remove = c("_baf","_tpm"),
           HLA_classification = class,
           prefix_classification = c("P1","P2"),
           qc_by = "_count",
           qc_threshold = 0)

# -- ----------------------------------------------------------------------

Cell_Class <- function(data,
                       col_remove,
                       HLA_classification,
                       qc = F,
                       max_different = 2,
                       qc_by = NULL,
                       qc_threshold = NULL){
  #创建结果存放目录
  if(!dir.exists(".//results")){
    dir.create(".//results")
  }
  
  # - -----------------------------------------------------------------------
  
  #质控数据
   #print("Data QC")
  # if(qc == T){
  #   source("RScript/Table2CleanData.R")
  #   cell.table <- get_cleandata(data = data,
  #                               col_remove = col_remove,
  #                               qc_by = qc_by,
  #                               qc_threshold = qc_threshold)
  #   write_csv(cell.table,"cell_table.csv")
  # }

  # - -----------------------------------------------------------------------
  
  #数据筛选, 得到来自于P1 P2和others的cell
  #获取每个cell的位点信息的cleandata
  cell_table <- 
    data %>%  
    select(-contains(col_remove)) %>% 
    gather(key = "HLA",value = "Type", -sample) %>% 
    select(sample,Type) %>% 
    drop_na() 
  #计算每个位点在cell中的出现次数
  cell_table_count <- 
    cell_table %>% 
    group_by(sample,Type) %>% 
    count() %>%
    rename("Counts" = "n")
  #提取位点名称
  type_name <- cell_table_count$Type
  #注释位点来源，P1 P2 P1_P2 others
  type_source <- HLA_classification[type_name,] %>% 
    as.data.frame() %>% 
    rename("source"=".")
  #将NA即没有注释到的位点转换为other
  type_source$source[is.na(type_source$source)] <- "other"
  #将位点来源写入cell_table_count并计算每个cell中各个来源的位点总数
  cell_table_count <- cbind(cell_table_count,source = type_source)
  #以cell和来源为factor计算counts
  cell_table_count <- 
    cell_table_count %>% 
    ungroup(Type) %>% 
    select(-Type) %>% 
    group_by(sample,source) %>% 
    summarise_each(sum) %>% 
    spread(key = source,value = Counts,fill = 0)
  write_csv(cell_table_count,"Cell_HLA_count.csv")

# -- ----------------------------------------------------------------------

  
  #进一步筛选数据
  source(".//RScript//Find_cell_Source.R")
  Find_Cell_Source(Class_data = cell_table_count,
                   max_different = max_different)
}



# -- ----------------------------------------------------------------------

library(tidyverse)

HLA_data <- read_tsv(".//data/test_input.tsv") 
class <- read_tsv(".//data/hla_known_new.txt") %>% column_to_rownames("type")
Cell_Class(data = HLA_data,col_remove = c("_baf","_tpm","_count"),HLA_classification = class)




co <- read_csv("Cell_Source.CSV")
co %>% mutate(class = case_when(
  source_num >
))



