fastqc  -t 16  *


cores=16

ref="/home/miniconda3/softwave/hg38.chromFa_22+xym/bowtie2/hg38.fa"

/home/miniconda3/softwave/SEACR-master/SEACR_1.3.sh


projPath="/data3/1_mageckenv/XSJ/cutrun_408"

mkdir -p ${projPath}/alignment/sam/bowtie2_summary
mkdir -p ${projPath}/alignment/bam
mkdir -p ${projPath}/alignment/bed
mkdir -p ${projPath}/alignment/bedgraph


for histName in "sg408-WJL_BKDL210030535-1a" "IgG_rep1" "IgG_rep2"
do
	echo ${histName}
	
	bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 1000 -p ${cores} -x ${ref} -1 ${projPath}/fastq/${histName}_R1.fastq -2 ${projPath}/fastq/${histName}_R2.fastq -S ${projPath}/alignment/sam/${histName}_bowtie2.sam &> ${projPath}/alignment/sam/bowtie2_summary/${histName}_bowtie2.txt

	
done
 
 
 
 
 
 
 

conda activate CUT_RUNTools_2
##=== R command ===## 
R
library(magrittr)  
library(dplyr)     
library(ggplot2)    
library(viridis)    
library(ggpubr) 
library(corrplot)

projPath="/data3/1_mageckenv/XSJ/cutrun_408"
sampleList = c("sg408-WJL_BKDL210030535-1a", "IgG_rep1", "IgG_rep2")
 
histList = c("sg408-WJL", "IgG", "IgG")
 
## Collect the alignment results from the bowtie2 alignment summary files
alignResult = c()
for(hist in sampleList){
  alignRes = read.table(paste0(projPath, "/alignment/sam/bowtie2_summary/", hist, "_bowtie2.txt"), header = FALSE, fill = TRUE)
  alignRate = substr(alignRes$V1[6], 1, nchar(as.character(alignRes$V1[6]))-1)
  histInfo = strsplit(hist, "_")[[1]]
  alignResult = data.frame(Histone = histInfo[1], Replicate = histInfo[2], 
                           SequencingDepth = alignRes$V1[1] %>% as.character %>% as.numeric, 
                           MappedFragNum_hg38 = alignRes$V1[4] %>% as.character %>% as.numeric + alignRes$V1[5] %>% as.character %>% as.numeric, 
                           AlignmentRate_hg38 = alignRate %>% as.numeric)  %>% rbind(alignResult, .)
}
alignResult$Histone = factor(alignResult$Histone)
df<- alignResult %>% mutate(AlignmentRate_hg38 = paste0(AlignmentRate_hg38, "%"))
write.csv(df,file="alignResult.csv")



##=== R command ===## 
## Generate sequencing depth boxplot
fig3A = alignResult %>% ggplot(aes(x = Histone, y = SequencingDepth/1000000, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Sequencing Depth per Million") +
    xlab("") + 
    ggtitle("A. Sequencing Depth")

fig3B = alignResult %>% ggplot(aes(x = Histone, y = MappedFragNum_hg38/1000000, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Mapped Fragments per Million") +
    xlab("") +
    ggtitle("B. Alignable Fragment (hg38)")

fig3C = alignResult %>% ggplot(aes(x = Histone, y = AlignmentRate_hg38, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("% of Mapped Fragments") +
    xlab("") +
    ggtitle("C. Alignment Rate (hg38)")



cairo_pdf(file="Generate_sequencing_depth_boxplot.pdf",width=10)  
print(ggarrange(fig3A, fig3B, fig3C,  ncol = 2, nrow=2, common.legend = TRUE, legend="bottom"))
dev.off()





##== linux command ==##
mkdir -p $projPath/alignment/sam/fragmentLen


for histName in "sg408-WJL_BKDL210030535-1a" "IgG_rep1" "IgG_rep2"
do
	echo ${histName}



	## Extract the 9th column from the alignment sam file which is the fragment length
	samtools view -F 0x04 $projPath/alignment/sam/${histName}_bowtie2.sam | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' >$projPath/alignment/sam/fragmentLen/${histName}_fragmentLen.txt



done

 


##=== R command ===## 
## Collect the fragment size information
fragLen = c()
for(hist in sampleList){
  
  histInfo = strsplit(hist, "_")[[1]]
  fragLen = read.table(paste0(projPath, "/alignment/sam/fragmentLen/", hist, "_fragmentLen.txt"), header = FALSE) %>% mutate(fragLen = V1 %>% as.numeric, fragCount = V2 %>% as.numeric, Weight = as.numeric(V2)/sum(as.numeric(V2)), Histone = histInfo[1], Replicate = histInfo[2], sampleInfo = hist) %>% rbind(fragLen, .) 
}
fragLen$sampleInfo = factor(fragLen$sampleInfo, levels = sampleList)
fragLen$Histone = factor(fragLen$Histone, levels = histList)
## Generate the fragment size density plot (violin plot)
fig5A = fragLen %>% ggplot(aes(x = sampleInfo, y = fragLen, weight = Weight, fill = Histone)) +
    geom_violin(bw = 5) +
    scale_y_continuous(breaks = seq(0, 800, 50)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 20) +
    ggpubr::rotate_x_text(angle = 20) +
    ylab("Fragment Length") +
    xlab("")

fig5B = fragLen %>% ggplot(aes(x = fragLen, y = fragCount, color = Histone, group = sampleInfo, linetype = Replicate)) +
  geom_line(size = 1) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500))

# ggarrange(fig5A, fig5B, ncol = 2)

pdf(file="fragment_size_information.pdf",width=20)  
print(ggarrange(fig5A, fig5B, ncol = 2))
dev.off()



##== linux command ==##
## Filter and keep the mapped read pairs


# for histName in "IgG_rep1" "IgG_rep2" "IgG_rep3" "IgG_rep4"
for histName in "sg408-WJL_BKDL210030535-1a" "IgG_rep1" "IgG_rep2"
do
	echo ${histName}



	samtools view -q 2 -bS -F 0x04 $projPath/alignment/sam/${histName}_bowtie2.sam >$projPath/alignment/bam/${histName}_bowtie2.mapped.bam

	## Convert into bed file format
	bedtools bamtobed -i $projPath/alignment/bam/${histName}_bowtie2.mapped.bam -bedpe >$projPath/alignment/bed/${histName}_bowtie2.bed

	## Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
	awk '$1==$4 && $6-$2 < 1000 {print $0}' $projPath/alignment/bed/${histName}_bowtie2.bed >$projPath/alignment/bed/${histName}_bowtie2.clean.bed

	## Only extract the fragment related columns
	cut -f 1,2,6 $projPath/alignment/bed/${histName}_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n  >$projPath/alignment/bed/${histName}_bowtie2.fragments.bed


done



##== linux command ==##
## We use the mid point of each fragment to infer which 500bp bins does this fragment belong to.

# for histName in "IgG_rep1" "IgG_rep2" "IgG_rep3" "IgG_rep4"
for histName in "sg408-WJL_BKDL210030535-1a" "IgG_rep1" "IgG_rep2"
# for histName in "K562_IgG_XLink.hg19"
do
	echo ${histName}

	binLen=500
	# awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' $projPath/alignment/bed/${histName}_bowtie2.fragments.bed | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n  >$projPath/alignment/bed/${histName}_bowtie2.fragmentsCount.bin$binLen.bed
	awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' $projPath/alignment/bed/${histName}_bowtie2.fragments.bed | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n  >$projPath/alignment/bed/${histName}_bowtie2.fragmentsCount.bin$binLen.bed

done



# ##== linux command ==##
# ## We use the mid point of each fragment to infer which 500bp bins does this fragment belong to.
# binLen=500


conda activate CUT_RUNTools_2
##== R command ==##
R
library(magrittr)  
library(dplyr)     
library(ggplot2)    
library(viridis)    
library(ggpubr) 
library(corrplot)
library(dplyr)
library(stringr)
library(ggplot2)
library(viridis)
library(GenomicRanges)
library(chromVAR) ## For FRiP analysis and differential analysis
library(DESeq2) ## For differential analysis section
library(ggpubr) ## For customizing figures
library(corrplot) ## For correlation plot







sampleList = c("sg408-WJL_BKDL210030535-1a", "IgG_rep1", "IgG_rep2")
# sampleList = c("IgG_rep1", "IgG_rep2")
projPath="/data3/1_mageckenv/XSJ/cutrun_408"

reprod = c()
fragCount = NULL
for(hist in sampleList){
  
  if(is.null(fragCount)){
    
    fragCount = read.table(paste0(projPath, "/alignment/bed/", hist, "_bowtie2.fragmentsCount.bin500.bed"), header = FALSE) 
    colnames(fragCount) = c("chrom", "bin", hist)
  
  }else{
    
    fragCountTmp = read.table(paste0(projPath, "/alignment/bed/", hist, "_bowtie2.fragmentsCount.bin500.bed"), header = FALSE)
    colnames(fragCountTmp) = c("chrom", "bin", hist)
    fragCount = full_join(fragCount, fragCountTmp, by = c("chrom", "bin"))
    
  }
}

M = cor(fragCount %>% select(-c("chrom", "bin")) %>% log2(), use = "complete.obs") 



pdf(file=paste0(length(sampleList),"fragment_size_information.pdf"),width=20)  
print(corrplot(M, method = "color", outline = T, addgrid.col = "darkgray", addrect = 3, rect.col = "black", rect.lwd = 3,cl.pos = "b", tl.col = "indianred4", tl.cex = 1, cl.cex = 1, addCoef.col = "black", number.digits = 2, number.cex = 1, col = colorRampPalette(c("midnightblue","white","darkred"))(100)))
dev.off()











##== linux command ==##
chromSize="/home/miniconda3/softwave/hg38.chromFa_22+xym/bowtie2/hg38.fa.fai"

mkdir -p $projPath/alignment/bedgraph

for histName in "sg408-WJL_BKDL210030535-1a" "IgG_rep1" "IgG_rep2"
# for histName in "K562_IgG_XLink.hg19"
do
	echo ${histName}

    bedtools genomecov -bg  -i $projPath/alignment/bed/${histName}_bowtie2.fragments.bed -g $chromSize > $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph
done




##== linux command ==##
seacr="/home/miniconda3/softwave/SEACR-master/SEACR_1.3.sh"

mkdir -p $projPath/peakCalling/SEACR

histControl="IgG_rep2"
histControl="IgG_rep1"

# for histName in "sg408-WJL_BKDL210030535-1a" "IgG_rep1" "IgG_rep2"


for histName in "sg408-WJL_BKDL210030535-1a"
# for histName in "K562_IgG_XLink.hg19"
do
	echo ${histName}


	bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph \
		 $projPath/alignment/bedgraph/${histControl}_bowtie2.fragments.normalized.bedgraph \
		 non stringent $projPath/peakCalling/SEACR/${histName}_seacr_control.peaks

	bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph 0.01 non stringent $projPath/peakCalling/SEACR/${histName}_seacr_top0.01.peaks



done




##=== R command ===## 
peakN = c()
peakWidth = c()
peakType = c("control", "top0.01")
for(hist in sampleList){
  histInfo = strsplit(hist, "_")[[1]]
  if(histInfo[1] != "IgG"){
    for(type in peakType){
      peakInfo = read.table(paste0(projPath, "/peakCalling/SEACR/", hist, "_seacr_", type, ".peaks.stringent.bed"), header = FALSE, fill = TRUE)  %>% mutate(width = abs(V3-V2))
      peakN = data.frame(peakN = nrow(peakInfo), peakType = type, Histone = histInfo[1], Replicate = histInfo[2]) %>% rbind(peakN, .)
      peakWidth = data.frame(width = peakInfo$width, peakType = type, Histone = histInfo[1], Replicate = histInfo[2])  %>% rbind(peakWidth, .)
    }
  }
}
peakN %>% select(Histone, Replicate, peakType, peakN)

write.csv(peakN,file="peakN.csv")


##=== R command ===## 
histL = c("sg408-WJL")
repL = paste0("rep", 1:2)
repL = c("BKDL210030535-1a")
peakType = c("control", "top0.01")
peakOverlap = c()
for(type in peakType){
  for(hist in histL){
    overlap.gr = GRanges()
    for(rep in repL){
      peakInfo = read.table(paste0(projPath, "/peakCalling/SEACR/", hist, "_", rep, "_seacr_", type, ".peaks.stringent.bed"), header = FALSE, fill = TRUE)
      peakInfo.gr = GRanges(peakInfo$V1, IRanges(start = peakInfo$V2, end = peakInfo$V3), strand = "*")
      if(length(overlap.gr) >0){
        overlap.gr = overlap.gr[findOverlaps(overlap.gr, peakInfo.gr)@from]
      }else{
        overlap.gr = peakInfo.gr
        
      }
    }
    peakOverlap = data.frame(peakReprod = length(overlap.gr), Histone = hist, peakType = type) %>% rbind(peakOverlap, .)
  }
}

peakReprod = left_join(peakN, peakOverlap, by = c("Histone", "peakType")) %>% mutate(peakReprodRate = peakReprod/peakN * 100)
peakReprod %>% select(Histone, Replicate, peakType, peakN, peakReprodNum = peakReprod, peakReprodRate)

write.csv(peakReprod,file="peakReprod.csv")




##=== R command ===## 
# BiocManager::install("chromVAR")

library(chromVAR)

bamDir = paste0(projPath, "/alignment/bam")
inPeakData = c()
## overlap with bam file to get count
for(hist in histL){
  for(rep in repL){
    peakRes = read.table(paste0(projPath, "/peakCalling/SEACR/", hist, "_", rep, "_seacr_control.peaks.stringent.bed"), header = FALSE, fill = TRUE)
    peak.gr = GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*")
    bamFile = paste0(bamDir, "/", hist, "_", rep, "_bowtie2.mapped.bam")
    fragment_counts <- getCounts(bamFile, peak.gr, paired = TRUE, by_rg = FALSE, format = "bam")
    inPeakN = counts(fragment_counts)[,1] %>% sum
    inPeakData = rbind(inPeakData, data.frame(inPeakN = inPeakN, Histone = hist, Replicate = rep))
  }
}

frip = left_join(inPeakData, alignResult, by = c("Histone", "Replicate")) %>% mutate(frip = inPeakN/MappedFragNum_hg38 * 100)
frip %>% select(Histone, Replicate, SequencingDepth, MappedFragNum_hg38, AlignmentRate_hg38, FragInPeakNum = inPeakN, FRiPs = frip)

write.csv(frip,file="frip.csv")



fig7A = peakN %>% ggplot(aes(x = Histone, y = peakN, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    facet_grid(~peakType) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Number of Peaks") +
    xlab("")

fig7B = peakWidth %>% ggplot(aes(x = Histone, y = width, fill = Histone)) +
    geom_violin() +
    facet_grid(Replicate~peakType) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    scale_y_continuous(trans = "log", breaks = c(400, 3000, 22000)) +
    theme_bw(base_size = 18) +
    ylab("Width of Peaks") +
    xlab("")

fig7C = peakReprod %>% ggplot(aes(x = Histone, y = peakReprodRate, fill = Histone, label = round(peakReprodRate, 2))) +
    geom_bar(stat = "identity") +
    geom_text(vjust = 0.1) +
    facet_grid(Replicate~peakType) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("% of Peaks Reproduced") +
    xlab("")

fig7D = frip %>% ggplot(aes(x = Histone, y = frip, fill = Histone, label = round(frip, 2))) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("% of Fragments in Peaks") +
    xlab("")


cairo_pdf(file=paste0(length(sampleList),"peakN.pdf"),width=10,height=12)  
	print(ggarrange(fig7A, fig7B, fig7C, fig7D, ncol = 2, nrow=2, common.legend = TRUE, legend="bottom"))
dev.off()


##== linux command ==##
mkdir -p $projPath/alignment/bigwig  

# for histName in "sg408-WJL_BKDL210030535-1a"
# # for histName in "K562_IgG_XLink.hg19"
# do


# done

histName="IgG_rep1"

echo ${histName}

																																	 
samtools sort -@ 16 -o $projPath/alignment/bam/${histName}.sorted.bam $projPath/alignment/bam/${histName}_bowtie2.mapped.bam                                                     
samtools index $projPath/alignment/bam/${histName}.sorted.bam    


histName="IgG_rep2"

echo ${histName}

																																	 
samtools sort -@ 16 -o $projPath/alignment/bam/${histName}.sorted.bam $projPath/alignment/bam/${histName}_bowtie2.mapped.bam                                                     
samtools index $projPath/alignment/bam/${histName}.sorted.bam    











histName="sg408-WJL_BKDL210030535-1a"

echo ${histName}

																																	 
samtools sort -@ 16 -o $projPath/alignment/bam/${histName}.sorted.bam $projPath/alignment/bam/${histName}_bowtie2.mapped.bam                                                     
samtools index $projPath/alignment/bam/${histName}.sorted.bam    

bamCoverage -b $projPath/alignment/bam/${histName}.sorted.bam -o $projPath/alignment/bigwig/${histName}_raw.bw                                                             
bamCoverage  --binSize 10  --numberOfProcessors 16 --effectiveGenomeSize 2862010578 --normalizeUsing RPGC --outFileFormat bigwig  -b $projPath/alignment/bam/${histName}.sorted.bam -o $projPath/alignment/bigwig/${histName}_RPGC.bw                                         




##== linux command ==##

mkdir -p  $projPath/data/hg38_gene

# 7.2.1 Heatmap over transcription units

cores=16
computeMatrix scale-regions -S $projPath/alignment/bigwig/sg408-WJL_BKDL210030535-1a_raw.bw \
							  -R /home/miniconda3/softwave/hg38.chromFa_22+xym/knownGene_hg38.bed \
							  --beforeRegionStartLength 3000 \
							  --regionBodyLength 5000 \
							  --afterRegionStartLength 3000 \
							  --skipZeros -o $projPath/data/hg38_gene/matrix_gene.mat.gz -p $cores

plotHeatmap -m $projPath/data/hg38_gene/matrix_gene.mat.gz -out $projPath/data/hg38_gene/Histone_gene.png --sortUsing sum




# 7.2.2. Heatmap on CUT&Tag peaks

##== linux command ==##

repL = c("BKDL210030535-1a")
repName="BKDL210030535-1a"

awk '{split($6, summit, ":"); split(summit[2], region, "-"); print summit[1]"\t"region[1]"\t"region[2]}' $projPath/peakCalling/SEACR/${histName}_seacr_control.peaks.stringent.bed >$projPath/peakCalling/SEACR/${histName}_seacr_control.peaks.summitRegion.bed

 

computeMatrix reference-point -S $projPath/alignment/bigwig/${histName}_raw.bw \
              -R $projPath/peakCalling/SEACR/${histName}_seacr_control.peaks.summitRegion.bed \
              --skipZeros -o $projPath/peakCalling/SEACR/${histName}_SEACR.mat.gz -p $cores -a 3000 -b 3000 --referencePoint center

plotHeatmap -m $projPath/peakCalling/SEACR/${histName}_SEACR.mat.gz -out $projPath/peakCalling/SEACR/${histName}_SEACR_heatmap.png --sortUsing sum --startLabel "Peak Start" -\
-endLabel "Peak End" --xAxisLabel "" --regionsLabel "Peaks" --samplesLabel "${histName} ${repName}"


 






