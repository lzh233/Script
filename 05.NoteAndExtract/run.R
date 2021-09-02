library(tidyverse)
rm(list = ls())
setwd("G:\\Desktop\\Github\\NoteAndExtract")
source(".//script/NoteAndExtract.R")

otu <- read.delim(".//data1/otu.csv",header = T,row.names = 1,encoding = "UTF-8",sep = ",")
note <- read.delim(".//data1/note.txt",header = T,row.names = 1,encoding = "UTF-8")

NoteAndExtract(otutable = otu,
               note.table = note,
               extract.by = "all",
               caclulate.ra = T,
               verbose = T)

NoteAndExtractSelect(otutable = otu,
                     note.table = note,
                     extract.select = c("Phylum","Family","Class"),
                     calculate.ra = T,
                     verbose = T)

NoteAndExtractAll(otu.dir = "./data2/",
                  note = note,
                  extract.by = "all",
                  calculate.ra = T,
                  sep.otu = ",",
                  verbose = T)
