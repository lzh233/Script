########################################
#Author:liuzihao                       #
#Date: 20210510                        #
#fuction: Extracting data from otutable#
#Version: 1.0                          #
########################################
#source("MainSet.R")
#输出信息补齐长度与符号
print.len <- 62
print.pad <- "-"

# - -----------------------------------------------------------------------


#创建目录
if(dir.exists(".\\extractData\\") == FALSE) {
  dir.create(".\\extractData\\")
  print("result directory creat")
} else {
  print("result directory exits!")
}

#输出目录设定
data.save <- ".\\extractData\\"


# -- ----------------------------------------------------------------------



#定义提取和计算丰度的功能
calculate.fetures <- function(otutable,
                             table,
                             extract.by,
                             calculate.ra){
  #将分组名称转换为内部变量名extract.factor
  names(table) <- "extract.factor"
  #data.save <- ".\\extractData\\"
  #以extract.factor为因子进行分组求和
  otutable_extract_f <- 
    data.frame(table,otutable) %>% 
    group_by(extract.factor) %>% 
    summarise_each(sum) %>% 
    column_to_rownames("extract.factor")
  
  #计算丰度
  if (calculate.ra == TRUE){
    otutable_extract_f_RA <- apply(otutable_extract_f, 2, function(x){x/sum(x)})
    write.csv(otutable_extract_f_RA,str_c(data.save,extract.by,"_RA.csv"))
  }
  
  
  #保存丰度表与计数表
  
  write.csv(otutable_extract_f,str_c(data.save,extract.by,".csv"))
}


# -- ----------------------------------------------------------------------


NoteAndExtract <- function(otutable,
                           note.table,
                           extract.by = "all",
                           calculate.ra = TRUE,
                           verbose = TRUE){
  if(verbose == TRUE){
    print(str_pad(" Extracting Features",width = print.len,side = "both",pad = print.pad))
  }  
    #输入行名为features，列名为处理
    otutable_ex <- 
      otutable %>% 
      as.data.frame()
    if(verbose == TRUE){
      print(str_c("Nmuber of features: ",nrow(otutable_ex) ))
    }
    
    
    #根据注释文件进行注释
    otutable_note <- note.table[rownames(otutable),]
    write.csv(otutable_note,str_c(data.save,"Features_Note.csv"))
  
  #根据extract.by指定的信息进行提取
  if (extract.by != "all") {
    if(verbose == TRUE ){
      print(str_pad(str_c("Features Extract by ",extract.by),width = print.len/2,side = "right",pad = " "))
    }
    #将对应分组进行提取并转换为数据框
    otutable_extract <- otutable_note[,extract.by] %>% as.data.frame()
    #调用calculate.fetures进行features提取与计算丰度
    calculate.fetures(otutable = otutable,table = otutable_extract,extract.by = extract.by,calculate.ra = calculate.ra)
  }
    
    
  if (extract.by == "all") {
    facs <- names(note.table)
    
    for (ex.factor in facs){
      if(verbose == TRUE){
        print(str_pad(str_c("Features Extract by ",ex.factor),width = print.len/2,side = "right",pad = " "))
      }
      #将对应分组进行提取并转换为数据框
      otutable_extract <- otutable_note[,ex.factor] %>% as.data.frame()
      
      #调用calculate.fetures进行features提取与计算丰度
      calculate.fetures(otutable = otutable,table = otutable_extract,extract.by = ex.factor,calculate.ra = calculate.ra)
    } 
  }
}


# -- ----------------------------------------------------------------------


NoteAndExtractAll <- function(otu.dir,
                              note,
                              extract.by = "all",
                              calculate.ra = TRUE,
                              sep.otu = "\t",
                              verbose = TRUE){
  
  otu.files <- dir(otu.dir)
  
  for(i in otu.files){
    dir.name <- str_split(i,pattern = "\\.",simplify = T,n = 2)[1]
    if(dir.exists(str_glue(".\\extractData\\{dir.name}")) == FALSE){
      dir.create(str_glue(".\\extractData\\{dir.name}"))
    }else{
      if(verbose == TRUE){print("result directory exits!")}
    }
    
    data.save <<- str_glue(".\\extractData\\{dir.name}\\")
    otu <- read.delim(str_glue("{otu.dir}\\{i}"),header = T,row.names = 1,sep = sep.otu,encoding = "UTF-8")
    NoteAndExtract(otutable = otu,note.table = note,extract.by = extract.by,calculate.ra = calculate.ra,verbose = verbose)
  }
  data.save <<- ".\\extractData\\"
}


# -- ----------------------------------------------------------------------

NoteAndExtractSelect <- function(otutable,
                                 note.table,
                                 extract.select = c("Phylum","Family"),
                                 calculate.ra = TRUE,
                                 verbose = TRUE){
  data.save <<- ".\\extractData\\"
  for (i in extract.select){
    NoteAndExtract(otutable = otutable,note.table=note.table,extract.by = i,calculate.ra = TRUE,verbose = verbose)
  }
}








