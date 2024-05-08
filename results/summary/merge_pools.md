Merge barcode-variant tables for libraries combined into different
experimental pools
================
Tyler Starr
02/20/2024

This notebook takes different pools as defined in the config.yaml file
and makes a merged barcode-variant lookup table that pools different
libraries according to experimental specifications.

``` r
#list of packages to install/load
packages = c("yaml","data.table","tidyverse","gridExtra","plotly","withr","htmlwidgets","knitr")
#install any packages not already installed
installed_packages <- packages %in% rownames(installed.packages())
if(any(installed_packages == F)){
  install.packages(packages[!installed_packages],
                   lib=c(paste("/uufs/chpc.utah.edu/common/home/",Sys.getenv("USER"),"/RLibs/",Sys.getenv("R_VERSION"),sep="")),
                   repos=c("http://cran.us.r-project.org"))
}
#load packages
invisible(lapply(packages, library, character.only=T))

knitr::opts_chunk$set(echo = T)
knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))

#read in config file
config <- read_yaml("config.yaml")
```

Session info for reproducing environment:

``` r
sessionInfo()
```

    ## R version 4.1.3 (2022-03-10)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Rocky Linux 8.8 (Green Obsidian)
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /uufs/chpc.utah.edu/sys/spack/linux-rocky8-nehalem/gcc-8.5.0/intel-oneapi-mkl-2021.4.0-h43nkmwzvaltaa6ii5l7n6e7ruvjbmnv/mkl/2021.4.0/lib/intel64/libmkl_rt.so.1
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] knitr_1.37        htmlwidgets_1.5.4 withr_2.5.0       plotly_4.10.1    
    ##  [5] gridExtra_2.3     forcats_0.5.1     stringr_1.4.0     dplyr_1.0.8      
    ##  [9] purrr_0.3.4       readr_2.1.2       tidyr_1.2.0       tibble_3.1.6     
    ## [13] ggplot2_3.4.1     tidyverse_1.3.1   data.table_1.14.2 yaml_2.3.5       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.1.2  xfun_0.30         haven_2.4.3       colorspace_2.0-3 
    ##  [5] vctrs_0.5.2       generics_0.1.2    viridisLite_0.4.0 htmltools_0.5.2  
    ##  [9] utf8_1.2.2        rlang_1.0.6       pillar_1.7.0      glue_1.6.2       
    ## [13] DBI_1.1.2         dbplyr_2.1.1      modelr_0.1.8      readxl_1.3.1     
    ## [17] lifecycle_1.0.3   munsell_0.5.0     gtable_0.3.0      cellranger_1.1.0 
    ## [21] rvest_1.0.2       evaluate_0.15     tzdb_0.2.0        fastmap_1.1.0    
    ## [25] fansi_1.0.2       broom_0.7.12      Rcpp_1.0.11       backports_1.4.1  
    ## [29] scales_1.2.1      jsonlite_1.8.7    fs_1.5.2          hms_1.1.1        
    ## [33] digest_0.6.29     stringi_1.7.6     grid_4.1.3        cli_3.6.0        
    ## [37] tools_4.1.3       magrittr_2.0.2    lazyeval_0.2.2    crayon_1.5.0     
    ## [41] pkgconfig_2.0.3   ellipsis_0.3.2    xml2_1.3.3        reprex_2.0.1     
    ## [45] lubridate_1.8.0   rstudioapi_0.13   assertthat_0.2.1  rmarkdown_2.13   
    ## [49] httr_1.4.7        R6_2.5.1          compiler_4.1.3

Read in tables of per-library barcode-variant lookups, and store the
underlying “library” in a new column called “sublibrary”

``` r
#read in barcode-variant lookup tables
dt_MERS_rpk <- data.table(read.csv(file=config$codon_variant_table_file_MERS_rpk,stringsAsFactors=F)); dt_MERS_rpk[,sublibrary:=library]
dt_PDF2180 <- data.table(read.csv(file=config$codon_variant_table_file_PDF2180,stringsAsFactors=F)); dt_PDF2180[,sublibrary:=library]
dt_panmerbeco <- data.table(read.csv(file=config$nt_variant_table_file_panmerbeco,stringsAsFactors=F)); dt_panmerbeco[,sublibrary:=library]
#remove panmerbeco wts with muts
dt_panmerbeco <- dt_panmerbeco[has_substitutions=="False",]

#for panmerbeco, need to add in the manually associated barcodes, and update column heads to match the codon variant tables
dt_manual <- data.table(read.csv(file=config$spike_in_barcodes,stringsAsFactors=F)); dt_manual[,sublibrary:=library]
dt_panmerbeco <- rbind(dt_panmerbeco,dt_manual)

dt_panmerbeco[,codon_substitutions:=""]
dt_panmerbeco[,aa_substitutions:=""]
dt_panmerbeco[,n_codon_substitutions:=0]
dt_panmerbeco[,n_aa_substitutions:=0]

dt <- rbindlist(list(dt_MERS_rpk,
                     dt_PDF2180,
                     dt_panmerbeco[,.(target, library, barcode, variant_call_support, codon_substitutions, aa_substitutions, n_codon_substitutions, n_aa_substitutions, sublibrary)]),
                use.names=T)

head(dt)
```

    ##      target library          barcode variant_call_support codon_substitutions
    ## 1: MERS_rpk   lib84 AAAAAAAAAAAACCTG                    1             GAT8AAA
    ## 2: MERS_rpk   lib84 AAAAAAAAAAGAATTA                    1 GCC106CCC GTT179GGT
    ## 3: MERS_rpk   lib84 AAAAAAAAAATGTGAA                    5            AGC43TAT
    ## 4: MERS_rpk   lib84 AAAAAAAAACCCCTGA                   12            TTA41CCA
    ## 5: MERS_rpk   lib84 AAAAAAAAACTCTAAA                    1            TGC49CAT
    ## 6: MERS_rpk   lib84 AAAAAAAAATACATAG                    2           CAA192ATG
    ##    aa_substitutions n_codon_substitutions n_aa_substitutions sublibrary
    ## 1:              D8K                     1                  1      lib84
    ## 2:      A106P V179G                     2                  2      lib84
    ## 3:             S43Y                     1                  1      lib84
    ## 4:             L41P                     1                  1      lib84
    ## 5:             C49H                     1                  1      lib84
    ## 6:            Q192M                     1                  1      lib84

Make new tables for each pool (better to do ‘by pool’ than ‘by library’,
since each pool is unique but any library can in theory be part of
multiple pools). Note, this code must be manually updated when new pools
are designated in the config.yaml file.

``` r
#pool1
dt_pool1 <- dt[sublibrary %in% config$pool1]
dt_pool1[,library:="pool1"]
head(dt_pool1)
```

    ##      target library          barcode variant_call_support codon_substitutions
    ## 1: MERS_rpk   pool1 AAAAAAAAAAAACCTG                    1             GAT8AAA
    ## 2: MERS_rpk   pool1 AAAAAAAAAAGAATTA                    1 GCC106CCC GTT179GGT
    ## 3: MERS_rpk   pool1 AAAAAAAAAATGTGAA                    5            AGC43TAT
    ## 4: MERS_rpk   pool1 AAAAAAAAACCCCTGA                   12            TTA41CCA
    ## 5: MERS_rpk   pool1 AAAAAAAAACTCTAAA                    1            TGC49CAT
    ## 6: MERS_rpk   pool1 AAAAAAAAATACATAG                    2           CAA192ATG
    ##    aa_substitutions n_codon_substitutions n_aa_substitutions sublibrary
    ## 1:              D8K                     1                  1      lib84
    ## 2:      A106P V179G                     2                  2      lib84
    ## 3:             S43Y                     1                  1      lib84
    ## 4:             L41P                     1                  1      lib84
    ## 5:             C49H                     1                  1      lib84
    ## 6:            Q192M                     1                  1      lib84

``` r
#pool2
dt_pool2 <- dt[sublibrary %in% config$pool2]
dt_pool2[,library:="pool2"]
head(dt_pool2)
```

    ##      target library          barcode variant_call_support codon_substitutions
    ## 1: MERS_rpk   pool2 AAAAAAAAAAAACGAA                    4           ACC113GCT
    ## 2: MERS_rpk   pool2 AAAAAAAAAAACATAT                    5            ATT66ACT
    ## 3: MERS_rpk   pool2 AAAAAAAAAAATGCAA                    1           GTT185GAT
    ## 4: MERS_rpk   pool2 AAAAAAAAAACCACAA                    4           ACA184GAT
    ## 5: MERS_rpk   pool2 AAAAAAAAACAAATGG                    1                    
    ## 6: MERS_rpk   pool2 AAAAAAAAACACACAA                    2             GCG2ATG
    ##    aa_substitutions n_codon_substitutions n_aa_substitutions sublibrary
    ## 1:            T113A                     1                  1      lib85
    ## 2:             I66T                     1                  1      lib85
    ## 3:            V185D                     1                  1      lib85
    ## 4:            T184D                     1                  1      lib85
    ## 5:                                      0                  0      lib85
    ## 6:              A2M                     1                  1      lib85

``` r
#pool3
dt_pool3 <- dt[sublibrary %in% config$pool3]
dt_pool3[,library:="pool3"]
head(dt_pool3)
```

    ##     target library          barcode variant_call_support codon_substitutions
    ## 1: PDF2180   pool3 AAAAAAAAAAGCCTCC                    3           GTT137TTT
    ## 2: PDF2180   pool3 AAAAAAAAAGTCAAAC                    2           AGC117AAA
    ## 3: PDF2180   pool3 AAAAAAAAAGTCACCG                    1            CAG95TGT
    ## 4: PDF2180   pool3 AAAAAAAACGTCCTAG                    1            CTG35GTT
    ## 5: PDF2180   pool3 AAAAAAAACTACATGT                    7           GCC143GAT
    ## 6: PDF2180   pool3 AAAAAAAACTAGATTA                    7  GCT58--- AGC111AGT
    ##    aa_substitutions n_codon_substitutions n_aa_substitutions sublibrary
    ## 1:            V137F                     1                  1      lib53
    ## 2:            S117K                     1                  1      lib53
    ## 3:             Q95C                     1                  1      lib53
    ## 4:             L35V                     1                  1      lib53
    ## 5:            A143D                     1                  1      lib53
    ## 6:             A58-                     2                  1      lib53

``` r
#pool4
dt_pool4 <- dt[sublibrary %in% config$pool4]
dt_pool4[,library:="pool4"]
head(dt_pool4)
```

    ##     target library          barcode variant_call_support codon_substitutions
    ## 1: PDF2180   pool4 AAAAAAAAAATACCGG                    2            GAT50ACT
    ## 2: PDF2180   pool4 AAAAAAAAAGATAGCC                    2           GAA200AAA
    ## 3: PDF2180   pool4 AAAAAAAAAGCCCATA                    3            ACT10AGA
    ## 4: PDF2180   pool4 AAAAAAAAATGTACTA                    2            ACA66TGT
    ## 5: PDF2180   pool4 AAAAAAAAATTTGAAT                    4           ACC176TGG
    ## 6: PDF2180   pool4 AAAAAAAACATCGAAC                    5           GCT199CAT
    ##    aa_substitutions n_codon_substitutions n_aa_substitutions sublibrary
    ## 1:             D50T                     1                  1      lib54
    ## 2:            E200K                     1                  1      lib54
    ## 3:             T10R                     1                  1      lib54
    ## 4:             T66C                     1                  1      lib54
    ## 5:            T176W                     1                  1      lib54
    ## 6:            A199H                     1                  1      lib54

``` r
#pool5
dt_pool5 <- dt[sublibrary %in% config$pool5]
dt_pool5[,library:="pool5"]
head(dt_pool5)
```

    ##      target library          barcode variant_call_support codon_substitutions
    ## 1: MERS_rpk   pool5 AAAAAAAAAAAACCTG                    1             GAT8AAA
    ## 2: MERS_rpk   pool5 AAAAAAAAAAGAATTA                    1 GCC106CCC GTT179GGT
    ## 3: MERS_rpk   pool5 AAAAAAAAAATGTGAA                    5            AGC43TAT
    ## 4: MERS_rpk   pool5 AAAAAAAAACCCCTGA                   12            TTA41CCA
    ## 5: MERS_rpk   pool5 AAAAAAAAACTCTAAA                    1            TGC49CAT
    ## 6: MERS_rpk   pool5 AAAAAAAAATACATAG                    2           CAA192ATG
    ##    aa_substitutions n_codon_substitutions n_aa_substitutions sublibrary
    ## 1:              D8K                     1                  1      lib84
    ## 2:      A106P V179G                     2                  2      lib84
    ## 3:             S43Y                     1                  1      lib84
    ## 4:             L41P                     1                  1      lib84
    ## 5:             C49H                     1                  1      lib84
    ## 6:            Q192M                     1                  1      lib84

``` r
#pool6
dt_pool6 <- dt[sublibrary %in% config$pool6]
dt_pool6[,library:="pool6"]
head(dt_pool6)
```

    ##      target library          barcode variant_call_support codon_substitutions
    ## 1: MERS_rpk   pool6 AAAAAAAAAAAACGAA                    4           ACC113GCT
    ## 2: MERS_rpk   pool6 AAAAAAAAAAACATAT                    5            ATT66ACT
    ## 3: MERS_rpk   pool6 AAAAAAAAAAATGCAA                    1           GTT185GAT
    ## 4: MERS_rpk   pool6 AAAAAAAAAACCACAA                    4           ACA184GAT
    ## 5: MERS_rpk   pool6 AAAAAAAAACAAATGG                    1                    
    ## 6: MERS_rpk   pool6 AAAAAAAAACACACAA                    2             GCG2ATG
    ##    aa_substitutions n_codon_substitutions n_aa_substitutions sublibrary
    ## 1:            T113A                     1                  1      lib85
    ## 2:             I66T                     1                  1      lib85
    ## 3:            V185D                     1                  1      lib85
    ## 4:            T184D                     1                  1      lib85
    ## 5:                                      0                  0      lib85
    ## 6:              A2M                     1                  1      lib85

Eliminate barcodes that are repeated between different variants within a
single pool.

``` r
#pool1
duplicates_pool1 <- dt_pool1[duplicated(dt_pool1,by=c("barcode","library")),.(library,barcode)] #the data.table duplciates function annoyingly only flags the first of each duplicate so doesn't intrinsically allow removal of both of the entries of the duplicate. So, flag what are duplicates, and then remove
dt_pool1[,duplicate:=FALSE]
if(nrow(duplicates_pool1) > 0){
  for(i in 1:nrow(duplicates_pool1)){
    dt_pool1[library==duplicates_pool1[i,library] & barcode==duplicates_pool1[i,barcode],duplicate:=TRUE]
  }
}
dt_pool1 <- dt_pool1[duplicate==FALSE,]; dt_pool1[,duplicate:=NULL]
print(paste("Removed", nrow(duplicates_pool1), "repeated barcodes from pool1"))
```

    ## [1] "Removed 0 repeated barcodes from pool1"

``` r
#pool2
duplicates_pool2 <- dt_pool2[duplicated(dt_pool2,by=c("barcode","library")),.(library,barcode)] #the data.table duplciates function annoyingly only flags the first of each duplicate so doesn't intrinsically allow removal of both of the entries of the duplicate. So, flag what are duplicates, and then remove
dt_pool2[,duplicate:=FALSE]
if(nrow(duplicates_pool2) > 0){
  for(i in 1:nrow(duplicates_pool2)){
    dt_pool2[library==duplicates_pool2[i,library] & barcode==duplicates_pool2[i,barcode],duplicate:=TRUE]
  }
}
dt_pool2 <- dt_pool2[duplicate==FALSE,]; dt_pool2[,duplicate:=NULL]
print(paste("Removed", nrow(duplicates_pool2), "repeated barcodes from pool2"))
```

    ## [1] "Removed 1 repeated barcodes from pool2"

``` r
#pool3
duplicates_pool3 <- dt_pool3[duplicated(dt_pool3,by=c("barcode","library")),.(library,barcode)] #the data.table duplciates function annoyingly only flags the first of each duplicate so doesn't intrinsically allow removal of both of the entries of the duplicate. So, flag what are duplicates, and then remove
dt_pool3[,duplicate:=FALSE]
if(nrow(duplicates_pool3) > 0){
  for(i in 1:nrow(duplicates_pool3)){
    dt_pool3[library==duplicates_pool3[i,library] & barcode==duplicates_pool3[i,barcode],duplicate:=TRUE]
  }
}
dt_pool3 <- dt_pool3[duplicate==FALSE,]; dt_pool3[,duplicate:=NULL]
print(paste("Removed", nrow(duplicates_pool3), "repeated barcodes from pool3"))
```

    ## [1] "Removed 1 repeated barcodes from pool3"

``` r
#pool4
duplicates_pool4 <- dt_pool4[duplicated(dt_pool4,by=c("barcode","library")),.(library,barcode)] #the data.table duplciates function annoyingly only flags the first of each duplicate so doesn't intrinsically allow removal of both of the entries of the duplicate. So, flag what are duplicates, and then remove
dt_pool4[,duplicate:=FALSE]
if(nrow(duplicates_pool4) > 0){
  for(i in 1:nrow(duplicates_pool4)){
    dt_pool4[library==duplicates_pool4[i,library] & barcode==duplicates_pool4[i,barcode],duplicate:=TRUE]
  }
}
dt_pool4 <- dt_pool4[duplicate==FALSE,]; dt_pool4[,duplicate:=NULL]
print(paste("Removed", nrow(duplicates_pool4), "repeated barcodes from pool4"))
```

    ## [1] "Removed 0 repeated barcodes from pool4"

``` r
#pool5
duplicates_pool5 <- dt_pool5[duplicated(dt_pool5,by=c("barcode","library")),.(library,barcode)] #the data.table duplciates function annoyingly only flags the first of each duplicate so doesn't intrinsically allow removal of both of the entries of the duplicate. So, flag what are duplicates, and then remove
dt_pool5[,duplicate:=FALSE]
if(nrow(duplicates_pool5) > 0){
  for(i in 1:nrow(duplicates_pool5)){
    dt_pool5[library==duplicates_pool5[i,library] & barcode==duplicates_pool5[i,barcode],duplicate:=TRUE]
  }
}
dt_pool5 <- dt_pool5[duplicate==FALSE,]; dt_pool5[,duplicate:=NULL]
print(paste("Removed", nrow(duplicates_pool5), "repeated barcodes from pool5"))
```

    ## [1] "Removed 4 repeated barcodes from pool5"

``` r
#pool6
duplicates_pool6 <- dt_pool6[duplicated(dt_pool6,by=c("barcode","library")),.(library,barcode)] #the data.table duplciates function annoyingly only flags the first of each duplicate so doesn't intrinsically allow removal of both of the entries of the duplicate. So, flag what are duplicates, and then remove
dt_pool6[,duplicate:=FALSE]
if(nrow(duplicates_pool6) > 0){
  for(i in 1:nrow(duplicates_pool6)){
    dt_pool6[library==duplicates_pool6[i,library] & barcode==duplicates_pool6[i,barcode],duplicate:=TRUE]
  }
}
dt_pool6 <- dt_pool6[duplicate==FALSE,]; dt_pool6[,duplicate:=NULL]
print(paste("Removed", nrow(duplicates_pool6), "repeated barcodes from pool6"))
```

    ## [1] "Removed 11 repeated barcodes from pool6"

Merge the per-pool tables back into one aggregate data table and save.

``` r
dt_final <- rbind(dt_pool1, dt_pool2, dt_pool3, dt_pool4, dt_pool5, dt_pool6)

dt_final %>%
  write.csv(file=config$codon_variant_table_file_pools, row.names=F)

head(dt_final)
```

    ##      target library          barcode variant_call_support codon_substitutions
    ## 1: MERS_rpk   pool1 AAAAAAAAAAAACCTG                    1             GAT8AAA
    ## 2: MERS_rpk   pool1 AAAAAAAAAAGAATTA                    1 GCC106CCC GTT179GGT
    ## 3: MERS_rpk   pool1 AAAAAAAAAATGTGAA                    5            AGC43TAT
    ## 4: MERS_rpk   pool1 AAAAAAAAACCCCTGA                   12            TTA41CCA
    ## 5: MERS_rpk   pool1 AAAAAAAAACTCTAAA                    1            TGC49CAT
    ## 6: MERS_rpk   pool1 AAAAAAAAATACATAG                    2           CAA192ATG
    ##    aa_substitutions n_codon_substitutions n_aa_substitutions sublibrary
    ## 1:              D8K                     1                  1      lib84
    ## 2:      A106P V179G                     2                  2      lib84
    ## 3:             S43Y                     1                  1      lib84
    ## 4:             L41P                     1                  1      lib84
    ## 5:             C49H                     1                  1      lib84
    ## 6:            Q192M                     1                  1      lib84
