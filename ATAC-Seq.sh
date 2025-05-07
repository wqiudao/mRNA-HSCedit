
conda activate CUT_RUNTools_2

# conda install -c bioconda genrich  -y
# conda install -c bioconda picard  -y
# conda install ucsc-gtftogenepred -y
# conda install ucsc-genepredtobed -y
# conda install ucsc-bedtobigbed -y
 


conda activate pygenometracks

gtfToGenePred  -geneNameAsName2 -ignoreGroupsWithoutExons -genePredExt gencode.v38.basic.annotation.gtf annotation.genePred
awk '{ print $12"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"}' annotation.genePred > annotation.gene_name.genePred
genePredToBed annotation.gene_name.genePred annotation.bed12

# sort bed12
sort -k1,1 -k2,2n annotation.bed12 > hg38_annotation.sorted.bed

 
 




pyGenomeTracks  --tracks  config.ini -dpi 300 --region chr11:4181821-4220502 --outFileName test.pdf --width 20 --height 20 

pyGenomeTracks --tracks my_tra.ini --region SL3.0ch01:10000-50000 --outFileName test.pdf --height 20 --width 20


cores=16

ref="/home/miniconda3/softwave/hg38.chromFa_22+xym/bowtie2/hg38.fa"
 

projPath="/data3/1_mageckenv/XSJ/ATAC-SEQ"

mkdir -p ${projPath}/alignment/sam/bowtie2_summary
mkdir -p ${projPath}/alignment/bam
mkdir -p ${projPath}/alignment/bed
mkdir -p ${projPath}/alignment/bedgraph
mkdir -p ${projPath}/GenrichDir


for histName in "sg411-cd34_BKDL210030532-1a" "sg411-hudep2_BKDL210030533-1a"
do
	echo ${histName}
	
	bowtie2 --very-sensitive -k 10 -X 2000 -p ${cores} -x ${ref}  -1 ${projPath}/fastq/${histName}_R1.fastq.gz -2 ${projPath}/fastq/${histName}_R2.fastq.gz -S ${projPath}/alignment/sam/${histName}_bowtie2.sam &> ${projPath}/alignment/sam/bowtie2_summary/${histName}_bowtie2.txt
	
	# bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 1000 -p ${cores} -x ${ref} -1 ${projPath}/fastq/${histName}_R1.fastq -2 ${projPath}/fastq/${histName}_R2.fastq -S ${projPath}/alignment/sam/${histName}_bowtie2.sam &> ${projPath}/alignment/sam/bowtie2_summary/${histName}_bowtie2.txt
	
done
 
 
for histName in "sg411-cd34_BKDL210030532-1a" "sg411-hudep2_BKDL210030533-1a"
do
	echo ${histName}
	
	samtools view -@ ${cores} -b -h  ${projPath}/alignment/sam/${histName}_bowtie2.sam | samtools sort  -@ ${cores} -n -O BAM   -o  ${projPath}/alignment/sam/${histName}_bowtie2_sorted.bam
	
	# gatk3 MarkDuplicates -nt ${cores}  -I ${projPath}/alignment/sam/${histName}_bowtie2_sorted.bam  -M ${projPath}/alignment/sam/${histName}_dup_matrics.txt -O ${projPath}/alignment/sam/${histName}_MD.bam
	
	# gatk3 CollectInsertSizeMetrics -H ${projPath}/alignment/sam/${histName}_InsertSize.pdf -I ${projPath}/alignment/sam/${histName}_MD.bam -O ${projPath}/alignment/sam/${histName}_InsertSize.txt
	# bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 1000 -p ${cores} -x ${ref} -1 ${projPath}/fastq/${histName}_R1.fastq -2 ${projPath}/fastq/${histName}_R2.fastq -S ${projPath}/alignment/sam/${histName}_bowtie2.sam &> ${projPath}/alignment/sam/bowtie2_summary/${histName}_bowtie2.txt
	
done
 
 
# for histName in "sg411-cd34_BKDL210030532-1a" "sg411-hudep2_BKDL210030533-1a"
# do
	# echo ${histName}
	
	# # samtools view -@ ${cores} -b -h  ${projPath}/alignment/sam/${histName}_bowtie2.sam | samtools sort  -@ ${cores} -n -O BAM   -o  ${projPath}/alignment/sam/${histName}_bowtie2_sorted.bam
	
	# picard MarkDuplicates  -I ${projPath}/alignment/sam/${histName}_bowtie2_sorted.bam  -M ${projPath}/alignment/sam/${histName}_dup_matrics.txt -O ${projPath}/alignment/sam/${histName}_MD.bam
	
	# picard CollectInsertSizeMetrics    -H ${projPath}/alignment/sam/${histName}_InsertSize.pdf -I ${projPath}/alignment/sam/${histName}_MD.bam -O ${projPath}/alignment/sam/${histName}_InsertSize.txt

# done
 
 
 
 
projPath="/data3/1_mageckenv/XSJ/ATAC-SEQ"
 
BamDir=${projPath}/alignment/sam

for histName in "sg411-cd34_BKDL210030532-1a" "sg411-hudep2_BKDL210030533-1a"
do
	echo ${histName}
	Sample=${histName}
		 
	# bedtools intersect -nonamecheck -v -a ${BamDir}/${Sample}_bowtie2_sorted.bam -b ${Blacklist} > ${BamDir}/${Sample}_RMBL.bam
	samtools view -@ 16 -h -F 1804 -q 30 -O SAM ${BamDir}/${Sample}_bowtie2_sorted.bam | grep -v -e "^chrM" | samtools sort   -@ 16  -O BAM -o ${BamDir}/${Sample}_Filtered.bam 
	samtools index ${BamDir}/${Sample}_Filtered.bam 
	alignmentSieve --numberOfProcessors 8 --ATACshift -b ${BamDir}/${Sample}_Filtered.bam -o ${BamDir}/${Sample}_Shift.bam
	
	
	samtools sort   -@ 16  -O BAM -o ${BamDir}/${Sample}_ShiftSort.bam  ${BamDir}/${Sample}_Shift.bam
	
	samtools index ${BamDir}/${Sample}_ShiftSort.bam
	bamCoverage --binSize 10  --numberOfProcessors 16 --effectiveGenomeSize 2862010578 --normalizeUsing RPGC --outFileFormat bigwig -b ${BamDir}/${Sample}_ShiftSort.bam -o ${BamDir}/${Sample}.bw


	# samtools view -c ${BamDir}/${Sample}_Shift.bam
	# bedtools intersect -u -a ${BamDir}/${Sample}_Shift.bam -b ${PeakDir}/${Sample}_peaks.narrowPeak | samtools view -c > ${PeakDir}/${Sample}_peaks.narrowPeak.log

done
  
 
  
GenrichDir=${projPath}/GenrichDir
  
Genrich -t ${BamDir}/sg411-cd34_BKDL210030532-1a_bowtie2_sorted.bam,${BamDir}/sg411-hudep2_BKDL210030533-1a_bowtie2_sorted.bam -o ${GenrichDir}/sg411.narrowPeak -f ${GenrichDir}/sg411_pq.bed -k ${GenrichDir}/sg411_pileup_p.bed -b ${GenrichDir}/sg411_reads.bed -r -e chrM   -m 30 -j



bedtools sort -i ${GenrichDir}/sg411_reads.bed | bgzip > ${GenrichDir}/sg411_reads.bed.gz
tabix -p bed ${GenrichDir}/sg411_reads.bed.gz

# bedtools sort -i WT_reads.bed | bgzip > WT_reads_sorted.bed.gz
# tabix -p bed WT_reads_sorted.bed.gz





conda activate pygenometracks




projPath="/data3/1_mageckenv/XSJ/ATAC-SEQ"
 
bwDir=${projPath}/alignment/sam
TrackDir=${projPath}/alignment/sam


"sg411-cd34_BKDL210030532-1a" "sg411-hudep2_BKDL210030533-1a"

make_tracks_file --trackFiles ${bwDir}/sg408-WJL_BKDL210030535-1a_RPGC.bw ${bwDir}/sg411-cd34_BKDL210030532-1a.bw ${bwDir}/sg411-hudep2_BKDL210030533-1a.bw  /home/miniconda3/softwave/hg38.chromFa_22+xym/hg38_annotation.sorted.bed -o ${TrackDir}/Track1.ini

 
pyGenomeTracks --tracks ${TrackDir}/Track1.ini   --region chr11:5241591-5261967 -o image.pdf   --trackLabelFraction 0.2 --width 38 --dpi 300
pyGenomeTracks --tracks Track1.ini   --region chr11:5241591-5261967 -o image.png

file = dm3_subset_BDGP5.78.gtf.gz
height = 10
title = gtf from ensembl one entry per gene
merge_transcripts = true
prefered_name = gene_name
fontsize = 12
file_type = bed

 --trackLabelFraction 0.2



 



R

library(ChIPseeker)
library(ggupset)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(clusterProfiler)

peak_data<-"GenrichDir/sg411_reads.bed.gz"

peak <- readPeakFile("GenrichDir/sg411_reads.bed.gz")
peak

covplot(peak, weightCol="V5")
covplot(peak, weightCol="V5", chrs=c("chr17", "chr18"), xlim=c(4.5e7, 5e7))
covplot(peak, chrs=c("chr11"))

peakHeatmap(peak_data, TxDb=txdb, upstream=3000, downstream=3000, color="red")



plotAvgProf2(peak_data, TxDb=txdb, upstream=3000, downstream=3000, xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")



peakAnno <- annotatePeak(peak_data, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")


write.csv(peakAnno,file="peakAnno.csv")

cairo_pdf(file=paste0("sg411_","annotatePeak.pdf"),width=10,height=12)  
	 


plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)

vennpie(peakAnno)

upsetplot(peakAnno)

upsetplot(peakAnno, vennpie=TRUE)




plotDistToTSS(peakAnno, title="Distribution of transcription factor-binding loci\nrelative to TSS")


dev.off()







  
  
  
  
  