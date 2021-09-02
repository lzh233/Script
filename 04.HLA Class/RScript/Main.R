Cell_Class <- function(data,
                       col_remove,
                       HLA_classification,
                       max_different = 2,
                       qc = F,
                       qc_by = NULL,
                       qc_threshold = NULL){

# -- ----------------------------------------------------------------------

  
  #创建结果存放目录
  if(!dir.exists(".//results")){
    dir.create(".//results")
  }
  
  # - -----------------------------------------------------------------------
  
  #质控数据
  #print("Data QC")
   if(qc == T){
     source("RScript/Table2CleanData.R")
     cell_table <- get_qc(data = data,
                          col_remove = col_remove,
                          qc_by = qc_by,
                          qc_threshold = qc_threshold)
   }
  
  # - -----------------------------------------------------------------------
  
  #数据筛选, 得到来自于P1 P2和others的cell
  if(qc == F){
    cell_table <- 
      data %>%  
      select(-contains(col_remove)) %>% 
      gather(key = "HLA",value = "Type", -colnames(data[,1])) %>% 
      select(sample,Type) %>% 
      drop_na()
  }

  #计算每个位点在cell中的出现次数
  cell_table_count <- 
    cell_table %>% 
    group_by(sample,Type) %>% 
    count() %>%
    rename("Counts" = "n")
  write_tsv(cell_table_count,"results/Cell_count_message.tsv")
  
  #提取位点名称
  type_name <- cell_table_count$Type
  
  #注释位点来源，P1 P2 P1_P2 others
  type_source <- 
    HLA_classification[type_name,] %>% 
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
  
  # -- ----------------------------------------------------------------------
  
  
  #根据细胞各个来源的counts数注释细胞
  source(".//RScript//Find_cell_Source.R")
  Annote_Cell_Source(Class_data = cell_table_count,
                   max_different = max_different)
}