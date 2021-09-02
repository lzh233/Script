#对数据进行质控
get_qc <- function(data = data,col_remove = col_remove,qc_by = qc_by,qc_threshold = qc_threshold){
  cell_table_gene <- 
    data %>% 
    select(-contains(col_remove),-contains(qc_by)) %>% 
    gather(key = "HLA",value = "Type",-colnames(data[,1])) %>% 
    arrange(sample)
  
  #得到count数的信息
  cell_table_gene_count <- 
    data %>% 
    select(sample,contains(qc_by)) %>% 
    gather(key = "HLA",value = "count",-colnames(data[,1]))%>% 
    arrange(sample)
  #write_csv(cell.table.counts,"count.csv")
  #得到cleandata
  cell_table <- cell_table_gene %>% 
    mutate(count = cell_table_gene_count$count) %>% 
    drop_na() %>% 
    filter(count > qc_threshold) %>% 
    select(-count) 
  write_csv(cell_table,".//results//qc_Data.csv")
  return(cell_table)
}