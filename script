#1.Map the reads for each sample to the reference genome:
hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188044_chrX_1.fastq.gz -2 chrX_data/samples/ERR188044_chrX_2.fastq.gz -S ERR188044_chrX.sam

hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188104_chrX_1.fastq.gz -2 chrX_data/samples/ERR188104_chrX_2.fastq.gz -S ERR188104_chrX.sam

hisat2 -p 8 --dta -x ./chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188234_chrX_1.fastq.gz -2 chrX_data/samples/ERR188234_chrX_2.fastq.gz -S ERR188234_chrX.sam

hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188245_chrX_1.fastq.gz -2 chrX_data/samples/ERR188245_chrX_2.fastq.gz -S ERR188245_chrX.sam

hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188257_chrX_1.fastq.gz -2 chrX_data/samples/ERR188257_chrX_2.fastq.gz -S ERR188257_chrX.sam

hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188273_chrX_1.fastq.gz -2 chrX_data/samples/ERR188273_chrX_2.fastq.gz -S ERR188273_chrX.sam

hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188337_chrX_1.fastq.gz -2 chrX_data/samples/ERR188337_chrX_2.fastq.gz -S ERR188337_chrX.sam

hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188383_chrX_1.fastq.gz -2 chrX_data/samples/ERR188383_chrX_2.fastq.gz -S ERR188383_chrX.sam

hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188401_chrX_1.fastq.gz -2 chrX_data/samples/ERR188401_chrX_2.fastq.gz -S ERR188401_chrX.sam

hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188428_chrX_1.fastq.gz -2 chrX_data/samples/ERR188428_chrX_2.fastq.gz -S ERR188428_chrX.sam

hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188454_chrX_1.fastq.gz -2 chrX_data/samples/ERR188454_chrX_2.fastq.gz -S ERR188454_chrX.sam

hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR204916_chrX_1.fastq.gz -2 chrX_data/samples/ERR204916_chrX_2.fastq.gz -S ERR204916_chrX.sam

#2.Sort and convert the SAM files to BAM:
samtools sort -@ 8 -o bam/ERR188044_chrX.bam bam/ERR188044_chrX.sam
samtools sort -@ 8 -o bam/ERR188104_chrX.bam bam/ERR188104_chrX.sam
samtools sort -@ 8 -o bam/ERR188234_chrX.bam bam/ERR188234_chrX.sam
samtools sort -@ 8 -o bam/ERR188245_chrX.bam bam/ERR188245_chrX.sam
samtools sort -@ 8 -o bam/ERR188257_chrX.bam bam/ERR188257_chrX.sam
samtools sort -@ 8 -o bam/ERR188273_chrX.bam bam/ERR188273_chrX.sam
samtools sort -@ 8 -o bam/ERR188337_chrX.bam bam/ERR188337_chrX.sam
samtools sort -@ 8 -o bam/ERR188383_chrX.bam bam/ERR188383_chrX.sam
samtools sort -@ 8 -o bam/ERR188401_chrX.bam bam/ERR188401_chrX.sam
samtools sort -@ 8 -o bam/ERR188428_chrX.bam bam/ERR188428_chrX.sam
samtools sort -@ 8 -o bam/ERR188454_chrX.bam bam/ERR188454_chrX.sam
samtools sort -@ 8 -o bam/ERR204916_chrX.bam bam/ERR204916_chrX.sam

#3.Assemble transcripts for each sample
stringtie -p 8 -G chrX_data/genes/chrX.gtf -o \
ERR188044_chrX.gtf -l ERR188044 bam/ERR188044_chrX.bam

stringtie -p 8 -G chrX_data/genes/chrX.gtf -o \
ERR188104_chrX.gtf -l ERR188104 ERR188104_chrX.bam 

stringtie -p 8 -G chrX_data/genes/chrX.gtf -o \
ERR188234_chrX.gtf -l ERR188234 bam/ERR188234_chrX.bam

stringtie -p 8 -G chrX_data/genes/chrX.gtf -o \
ERR188245_chrX.gtf -l ERR188245 bam/ERR188245_chrX.bam

stringtie -p 8 -G chrX_data/genes/chrX.gtf -o \
ERR188257_chrX.gtf -l ERR188257 bam/ERR188257_chrX.bam

stringtie -p 8 -G chrX_data/genes/chrX.gtf -o \
ERR188273_chrX.gtf -l ERR188273 bam/ERR188273_chrX.bam

stringtie -p 8 -G chrX_data/genes/chrX.gtf -o \
ERR188337_chrX.gtf -l ERR188337 bam/ERR188337_chrX.bam

stringtie -p 8 -G chrX_data/genes/chrX.gtf -o \
ERR188383_chrX.gtf -l ERR188383 bam/ERR188383_chrX.bam

stringtie -p 8 -G chrX_data/genes/chrX.gtf -o \
ERR188401_chrX.gtf -l ERR188401 bam/ERR188401_chrX.bam

stringtie -p 8 -G chrX_data/genes/chrX.gtf -o \
ERR188428_chrX.gtf -l ERR188428 bam/ERR188428_chrX.bam

stringtie -p 8 -G chrX_data/genes/chrX.gtf -o \
ERR188454_chrX.gtf -l ERR188454 bam/ERR188454_chrX.bam

stringtie -p 8 -G chrX_data/genes/chrX.gtf -o \
ERR204916_chrX.gtf -l ERR204916 bam/ERR204916_chrX.bam

#4.Merge transcripts from all samples:
stringtie --merge -p 8 -G chrX_data/genes/chrX.gtf -o stringtie_merged.gtf chrX_data/mergelist.txt

#5.Examine how the transcripts compare with the reference annotation (optional):
gffcompare -r chrX_data/genes/chrX.gtf -G -o merged stringtie_merged.gtf

#6.Estimate transcript abundances and create table counts for Ballgown:
stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188044/ERR188044_chrX.gtf bam/ERR188044_chrX.bam

stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188104/ERR188104_chrX.gtf bam/ERR188104_chrX.bam

stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188234/ERR188234_chrX.gtf bam/ERR188234_chrX.bam

stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188245/ERR188245_chrX.gtf bam/ERR188245_chrX.bam

stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188257/ERR188257_chrX.gtf bam/ERR188257_chrX.bam

stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188273/ERR188273_chrX.gtf bam/ERR188273_chrX.bam

stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188337/ERR188337_chrX.gtf bam/ERR188337_chrX.bam

stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188383/ERR188383_chrX.gtf bam/ERR188383_chrX.bam

stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188401/ERR188401_chrX.gtf bam/ERR188401_chrX.bam

stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188428/ERR188428_chrX.gtf bam/ERR188428_chrX.bam

stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188454/ERR188454_chrX.gtf bam/ERR188454_chrX.bam

stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR204916/ERR204916_chrX.gtf bam/ERR204916_chrX.bam

#7.Load relevant R packages. These include the Ballgown package that you will use for performing most of the analyses, as well as a few other packages that help with these analyses, specifically RSkittleBrewer (for setting up colors), genefilter (for fast calculation of means and variances), dplyr (for sorting and arranging results) and devtools (for reproducibility and installing packages):
# R
# 下载RSkittleBrewer包的步骤
# install.packages("devtools")
# library(devtools)
# options(unzip = 'internal')
# install_github("alyssafrazee/RSkittleBrewer")
# >library(ballgown)
# >library(RSkittleBrewer)
# >library(genefilter)
# >library(dplyr)
# >library(devtools)

8.pheno_data = read.csv("chrX_data/geuvadis_phenodata.csv")

9.bg_chrX = ballgown(dataDir = "ballgown", samplePattern = "ERR", pData=pheno_data)

10.bg_chrX_filt = subset(bg_chrX,"rowVars(texpr(bg_chrX)) >1",genomesubset=TRUE)

11.results_transcripts = stattest(bg_chrX_filt,feature="transcript",covariate="sex",adjustvars =c("population"), getFC=TRUE, meas="FPKM")

12,results_genes = stattest(bg_chrX_filt, feature="gene",covariate="sex", adjustvars = c("population"), getFC=TRUE,meas="FPKM")

13. Add gene names and gene IDs to the results_transcripts data frame:
results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_chrX_filt),geneIDs=ballgown::geneIDs(bg_chrX_filt), results_transcripts)

14| Sort the results from the smallest P value to the largest:
results_transcripts = arrange(results_transcripts,pval)
results_genes = arrange(results_genes,pval)

15| Write the results to a csv file that can be shared and distributed:
write.csv(results_transcripts, "chrX_transcript_results.csv",row.names=FALSE)
write.csv(results_genes, "chrX_gene_results.csv",row.names=FALSE)

16| Identify transcripts and genes with a q value <0.05:
subset(results_transcripts,results_transcripts$qval<0.05)
subset(results_genes,results_genes$qval<0.05)

17.Make the plots pretty. This step is optional, but if you do run it you will get the plots in the nice colors that we used to generate our figures:
tropical= c('darkorange', 'dodgerblue','hotpink', 'limegreen', 'yellow')
palette(tropical)

18.fpkm = texpr(bg_chrX,meas="FPKM")
fpkm = log2(fpkm+1)
boxplot(fpkm,col=as.numeric(pheno_data$sex),las=2,ylab='log2(FPKM+1)')

19.ballgown::transcriptNames(bg_chrX)[12]
ballgown::geneNames(bg_chrX)[12]
plot(fpkm[12,] ~ pheno_data$sex, border=c(1,2),main=paste(ballgown::geneNames(bg_chrX)[12],' : ',ballgown::transcriptNames(bg_chrX)[12]),pch=19, xlab="Sex",ylab='log2(FPKM+1)')
points(fpkm[12,] ~ jitter(as.numeric(pheno_data$sex)),col=as.numeric(pheno_data$sex))\

20.plotTranscripts(ballgown::geneIDs(bg_chrX)[1729], bg_chrX, main=c('Gene XIST in sample ERR188234'), sample=c('ERR188234'))

21.plotMeans('MSTRG.56', bg_chrX_filt,groupvar="sex",legend=FALSE)
