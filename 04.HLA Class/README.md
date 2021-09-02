# 应用场景

根据单细胞HLA分型的结果，来判断每一个细胞是来自于供体还是受体。

# 输入文件

- `reads.tsv`: 每个cell barcode的HLA在每个位点的分型结果，**注意第一列不能该转换为行名**，使用默认的读入数据的函数时，不能指定`row.names`

- `HLA_classification.txt`: 已知的两个人（供体和受体）的HLA在每个位点的分型，列名不能重复，需要通过一些字段来区分供体和受体

# 输入文件格式

- 需要两个输入文件，并将其置于`data`目录中

  `reads.tsv`: 每个cell barcode的HLA在每个位点的分型结果，该表格的第一列的为`cellbarcode`，第一列的列名需要命名为`sample`

  ```r
  # A tibble: 5 x 7
    sample                   A_allele1 A_allele1_count A_allele1_tpm A_allele2 A_allele2_count A_allele2_tpm
    <chr>                    <chr>               <dbl>         <dbl> <chr>               <dbl>         <dbl>
  1 CATCAAGTAGCACCTCAATGTTGC NA                     NA            NA NA                     NA            NA
  2 ACGTATCACAGATCTGCATACCAA NA                     NA            NA NA                     NA            NA
  3 CGAACTTATGAAGAGACTGAGCCA A*01:01                 0             0 A*11:01                 4         13157
  4 ATTGAGGAAACAACCAAACAACCA A*02:19                 1           715 A*02:804                0            NA
  5 CAAGGAGCACAGCAGAACCACTGT A*30:01                 1          1184 A*30:01                 1          1184
  ```
  
   `HLA_classification.txt`: 已知的两个人（供体和受体）的HLA在每个位点的分型，具体的输入格式如下表，该文件输入允许行数不等
  
  ```R
  # A tibble: 12 x 4
     P1_0702_CA1_gene1 P1_0702_CA1_gene2 P2_0702_CA1_gene1 P2_0702_CA1_gene2
     <chr>             <chr>             <chr>             <chr>            
   1 A*02:06           A*02:10           A*24:418          A*30:01          
   2 B*13:02           B*40:06           B*54:01           B*13:02          
   3 C*08:01           C*03:04           C*06:02           C*03:04          
   4 DQA1*05:05        DQA1*06:01        DQA1*02:01        DQA1*03:03       
   5 DQB1*03:01        DQB1*03:01        DQB1*02:02        DQB1*04:02  
  ```

# 输出文件格式

- 输出文件存放于脚本运行目录的`results`目录内，共有四个输出文件，

  `cell_HLA_Class.csv`: 每个细胞依据`HLA_classification.txt`的分型结果汇总

  `cell_P1.csv`: 来自于`P1`的细胞列表

  `cell_P2.csv`: 来自于`P2`的细胞列表

  `cell_Unclass.csv` : 无法判断来自于`P1`或是`P2`的细胞列表
  
  `HLA_CleanData.csv`: 质控与转换后的数据

# 实现方法

在`R_Script`中共有三个脚本，分别为，

- `Main.R`: 主函数，用于实现来判断每一个细胞是来自于供体还是受体，提供的函数为`Cell_Class()`，并生成输出文件`cell_HLA_Class.csv`
- `Table2CleanData.R`: 用于输入数据的质控和转换，即根据`counts`删除`counts = 0`的位点，提供函数为`get_cleandata `, 输出文件`HLA_CleanData.csv`
- `Find_Unclass_Cell.R`: 用于对最后结果的拆分，即确定哪些细胞来自于`P1`，哪些细胞来自于`P2`，哪些无法判断, 并生成输出文件`cell_P1.csv`、`cell_P2.csv`和`cell_Unclass.csv` 

## 具体实现算法

**Step1 （质控）**：将数据输入主函数后，`Cell_Class()`函数首先调用`get_cleandata()`对数据进行一个质控与转换，第一步，将数据首先转换为`tidydata`，进一步筛选`counts`数量大于0的位点，第二步将所有包含的`NA`的行删除，得到质控后的`cleandata`，并输出`HLA_CleanData.csv`

```r
# A tibble: 17,306 x 3
   sample                   HLA          type      
   <chr>                    <chr>        <chr>     
 1 AAACATCGAACAACCACCGAAGTA A_allele1    A*02:01   
 2 AAACATCGAACAACCACCGAAGTA A_allele2    A*02:01   
 3 AAACATCGAACAACCACCGAAGTA B_allele1    B*13:02   
 4 AAACATCGAACAACCACCGAAGTA B_allele2    B*40:06   
# ... with 17,296 more rows
```

**Step2  注释表`HLA_classification.txt`的处理**：根据`cell_Class()`函数中的`prefix_classification`提供的用于区分样品的`index`对注释表拆分，后根据拆分的结果匹配细胞信息

**Step3 （匹配细胞信息）：**首先取出一对等位基因中的第一个，根据Step1得到的`cleandata`,去匹配`type`列等于第一个基因的细胞，得到`cell_gene_1`，然后使用`cell_gene_1`的细胞名去匹配`cleandata`中的对应的细胞得到`cell_gene_2`, 然后使用等位基因中第二个基因`gene2`, 去匹配`cell_gene_2`的`type` 列，然后使用`inner_join()`函数将`cell_gene_2`和`cell_gene_1`共有的部分合并，即为含有这一对等位基因的细胞, 同时将每个细胞的来源标注，存储于`Source`列，通过循环将每一对等位基因进行匹配，将每次匹配的结果存放于`cell.filter`列表中，最后将列表中的每个元素合并成一个`Dataframe`，格式如下，

```r
# A tibble: 857 x 6
   sample                   HLA.x     type.x  HLA.y     type.y  source
   <chr>                    <chr>     <chr>   <chr>     <chr>   <chr> 
 1 GAATCTGAGTCTGTCAATCATTCC A_allele1 A*02:06 A_allele2 A*02:10 P1    
 2 AAACATCGAACAACCACCGAAGTA B_allele1 B*13:02 B_allele2 B*40:06 P1    
 3 AACAACCAACAGATTCCAACCACA B_allele1 B*13:02 B_allele2 B*40:06 P1    
 4 AAGACGGAAAGGTACAGTGTTCTA B_allele2 B*13:02 B_allele1 B*40:06 P1    
# ... with 847 more rows
```

**Step3-1 (等位基因名称一致)：** 当等位基因名称一致时，使用`inner_join()`会产生冗余数据，因此需要对数据进一步去冗余，去冗余方法为删除`HLA.x`和`HLA.y`相同的行，然后对细胞进行去重处理

**Step4 （对筛选出的细胞进一步筛分）：**调用`Find_Unclass_Cell.R`中的`Find_Unclass_Cell()`函数，根据细胞的来源对细胞区分，同时将无法区分的细胞单独保存，输出文件为，`cell_P1.csv`;  `cell_P2.csv`;  `cell_Unclass.csv` ，并且输出汇总数据`cell_HLA_Class.csv`，具体格式如下，

**cell_P1.csv**

```r
# A tibble: 206 x 2
   sample                   source
   <chr>                    <chr> 
 1 GAATCTGAGTCTGTCAATCATTCC P1    
 2 AAACATCGAACAACCACCGAAGTA P1    
 3 AAGGTACAAGCACCTCGGTGCGAA P1    
 4 ACATTGGCTATCAGCACATCAAGT P1   
# ... with 202 more rows
```

**cell_P2.csv**

```r
> aa
# A tibble: 404 x 2
   sample                   source
   <chr>                    <chr> 
 1 TTCACGCAAGGCTAACACTATGCA P2    
 2 AAACATCGGAACAGGCCACCTTAC P2    
 3 AACAACCAACAGATTCACCTCCAA P2    
 4 AACCGAGACTGAGCCAAAGAGATC P2    
# ... with 400 more rows
```

**cell_Unclass.csv**

```r
# A tibble: 100 x 2
   sample                   source 
   <chr>                    <chr>  
 1 AACAACCAACAGATTCCAACCACA unclass
 2 AAGACGGAAAGGTACAGTGTTCTA unclass
 3 ACTATGCAACATTGGCAACAACCA unclass
 4 ATTGAGGAAACGCTTAACAAGCTA unclass
# ... with 93 more rows
```

**cell_HLA_Class.csv**

```r
# A tibble: 857 x 6
   sample                   HLA.x     type.x  HLA.y     type.y  source
   <chr>                    <chr>     <chr>   <chr>     <chr>   <chr> 
 1 GAATCTGAGTCTGTCAATCATTCC A_allele1 A*02:06 A_allele2 A*02:10 P1    
 2 AAACATCGAACAACCACCGAAGTA B_allele1 B*13:02 B_allele2 B*40:06 P1    
 3 AAGGTACAAGCACCTCGGTGCGAA B_allele1 B*13:02 B_allele2 B*40:06 P1    
 4 ACATTGGCTATCAGCACATCAAGT B_allele1 B*13:02 B_allele2 B*40:06 P1    
# ... with 853 more rows
```

# 使用方法

```r
source(".//RScript/Main.R")
library(tidyverse)
cell_Class(data,
           col_remove,
           HLA_classification,
           prefix_classification,
           qc_by = "_count",
           qc_threshold = 0)
```

- 选项参数

  `data`: 输入数据，需要保证每列可以有明确的区分的字段

  `col_remove`：后跟字符串或字符串向量，需要移除的列的区分字段

  `HLA_classification`：已知的两个人（供体和受体）的HLA在每个位点的分型

  `prefix_classification`：后跟字符串向量，可以区分供体和受体的字段

  `qc_by`：哪一列为qc的指标，目前仅支持使用`counts`列

   `qc_threshold`: 设定质控的阈值，默认为0

# 测试数据

```r
setwd("yourPath")
#install.packages("tidyverse")
library(tidyverse)
HLA_data <- read_tsv("yourPath//data//reads.tsv") 
class <- read_tsv("yourPath//data/HLA_classification.txt")
#载入函数
source(".//RScript/Main.R")
cell_Class(data = HLA_data,
           col_remove = c("_baf","_tpm"),
           HLA_classification = class,
           prefix_classification = c("P1","P2"),
           qc_by = "_count",
           qc_threshold = 0)
###输出信息会简单显示来自于P1/P2/unclass的cell数量
[1] "P1: 206"
[1] "P2: 404"
[1] "Unclass: 97"
```





