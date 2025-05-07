# mRNA-HSCedit

# install
dnf install ncurses-compat-libs  -y #need root


``` 
conda create --name CUT_RUNTools_2 --file conda_spec_list.txt  -y

conda activate CUT_RUNTools_2

conda install -c bioconda samtools -y
conda install -c bioconda bedtools -y
conda install -c bioconda deeptools -y

conda install -c bioconda genrich  -y
conda install -c bioconda picard  -y
conda install ucsc-gtftogenepred -y
conda install ucsc-genepredtobed -y
conda install ucsc-bedtobigbed -y

pip install pyGenomeTracks

``` 
hg38-blacklist.v2.bed  https://github.com/Boyle-Lab/Blacklist

R
``` 
install.packages('XML')
install.packages("stringi")
install.packages("stringr")
install.packages("broom") 
install.packages("knitr") 
install.packages("numDeriv") 
install.packages("magrittr")  
install.packages("dplyr")     
install.packages("ggplot2")    
install.packages("viridis")    
install.packages("survival") 
install.packages("corrplot") 
install.packages("locfit") 
   
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")


BiocManager::install("GenomicRanges")
BiocManager::install("chromVAR")
BiocManager::install("rtracklayer")
BiocManager::install("BSgenome")
install.packages("https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.4.tar.gz", repos=NULL, type="source")
install.packages("https://cran.r-project.org/src/contrib/Archive/pbkrtest/pbkrtest_0.5.1.tar.gz", repos=NULL, type="source")
install.packages("https://cran.r-project.org/src/contrib/Archive/ggimage/ggimage_0.3.0.tar.gz", repos=NULL, type="source")
BiocManager::install("DESeq2")
install.packages("ggpubr")
BiocManager::install("ChIPseeker")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

```



# run script
```
reads_clean_pe.sh 13_BKDL210033002-1a_1.clean.fq.gz 13_BKDL210033002-1a_2.clean.fq.gz /data3/02-Data/20_clean_data
reads_clean_pe.sh 14_BKDL210033003-1a_1.clean.fq.gz 14_BKDL210033003-1a_2.clean.fq.gz /data3/02-Data/20_clean_data
reads_clean_pe.sh 15_BKDL210033004-1a_1.clean.fq.gz 15_BKDL210033004-1a_2.clean.fq.gz /data3/02-Data/20_clean_data
```


cutrun.sh

ATAC-Seq.sh



