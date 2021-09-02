Annote_Cell_Source <- function(Class_data,
                 max_different){
  Cell_Source <- 
    Class_data %>% 
    mutate(source_num = P1-P2) %>% 
    mutate(class = case_when(
      source_num >= max_different ~ "P1",
      source_num <= 0-max_different ~ "P2",
      TRUE ~ "Others"
    )) %>% 
    select(-source_num)
  write_csv(Cell_Source,".//results//Cell_Source.CSV")
}