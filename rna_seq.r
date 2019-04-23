
#7. 下载RSkittleBrewer包的步骤
install.packages("devtools")
library(devtools)
options(unzip = 'internal')
install_github("alyssafrazee/RSkittleBrewer")
library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)

#8. 读入数据
pheno_data = read.csv("chrX_data/geuvadis_phenodata.csv")

#9.dataDir指的是数据储存的地方，samplePattern则依据样本的名字来，pheno_data则指明了样本数据的关系，这个里面第一列样本名需要和ballgown下面的文件夹的样本名一样，不然会报错
bg_chrX = ballgown(dataDir = "ballgown", samplePattern = "ERR", pData=pheno_data)

#10.滤掉低丰度的基因
bg_chrX_filt = subset(bg_chrX,"rowVars(texpr(bg_chrX)) >1",genomesubset=TRUE)

#11.确认组间有差异的转录本，Ballgown的统计算法基于标准的线性模型，对于组内数据较少(<4)时较为准确。这里也可以使用其它的一些软件例如limma/DESeq/edgeR等。
results_transcripts = stattest(bg_chrX_filt,feature="transcript",covariate="sex",adjustvars =c("population"), getFC=TRUE, meas="FPKM")

#12,确认组间有差异的基因
results_genes = stattest(bg_chrX_filt, feature="gene",covariate="sex", adjustvars = c("population"), getFC=TRUE,meas="FPKM")

#13. 对结果 result_transcripts增加基因名
results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_chrX_filt),geneIDs=ballgown::geneIDs(bg_chrX_filt), results_transcripts)

#14.按照P值排序（从小到大）
results_transcripts = arrange(results_transcripts,pval)
results_genes = arrange(results_genes,pval)

#15.把结果写入csv文件
write.csv(results_transcripts, "chrX_transcript_results.csv",row.names=FALSE)
write.csv(results_genes, "chrX_gene_results.csv",row.names=FALSE)

#16.筛选出q值小于0.05的transcripts和genes，也就是在性别间表达有差异的基因了
subset(results_transcripts,results_transcripts$qval<0.05)
subset(results_genes,results_genes$qval<0.05)

#17.数据可视化之颜色设定
tropical= c('darkorange', 'dodgerblue','hotpink', 'limegreen', 'yellow')
palette(tropical)
# or rainbow()函数
#palette(rainbow(5))

#18.以FPKM为参考值作图，以性别作为区分条件
fpkm = texpr(bg_chrX,meas="FPKM")  # 提取FPKM值
fpkm = log2(fpkm+1)                #方便作图将其log转换，+1是为了避免出现log2(0)的情况
boxplot(fpkm,col=as.numeric(pheno_data$sex),las=2,ylab='log2(FPKM+1)')    #画图

#19.就单个转录本的查看在样品中的分布
ballgown::transcriptNames(bg_chrX)[12]
ballgown::geneNames(bg_chrX)[12]
# 绘制箱线图
plot(fpkm[12,] ~ pheno_data$sex, border=c(1,2),main=paste(ballgown::geneNames(bg_chrX)[12],' : ',ballgown::transcriptNames(bg_chrX)[12]),pch=19, xlab="Sex",ylab='log2(FPKM+1)')
points(fpkm[12,] ~ jitter(as.numeric(pheno_data$sex)),col=as.numeric(pheno_data$sex))\

#20.查看某一基因位置上所有的转录本
# plotTranscripts函数可以根据指定基因的id画出在特定区段的转录本
# 可以通过sample函数指定看在某个样本中的表达情况，这里选用id=1750, sample=ERR188234
plotTranscripts(ballgown::geneIDs(bg_chrX)[1729], bg_chrX, main=c('Gene XIST in sample ERR188234'), sample=c('ERR188234'))

#21.以性别为区分查看基因表达情况
# 这里以id=575的基因为例（对应上一步作图）
plotMeans('MSTRG.56', bg_chrX_filt,groupvar="sex",legend=FALSE)
