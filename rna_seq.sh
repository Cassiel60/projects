#1.将reads比对上genome
# 少量样本
#hisat2 -p 8 --dta -x chrX_data/indexes/chrX_tran -1 chrX_data/samples/ERR188044_chrX_1.fastq.gz -2 chrX_data/samples/ERR188044_chrX_2.fastq.gz -S ERR188044_chrX.sam
# 批量操作
for i in *1.fastq.gz; 
do
i=${i%1.fastq.gz*}; 
nohup hisat2 -p 8 --dta -x ~/RNA-seq/chrx/chrX_data/indexes/chrX_tran -1 ${i}1.fastq.gz -2 ${i}2.fastq.gz -S ~/RNA-seq/chrx/chrX_data/result/${i}align.sam 2>~/RNA-seq/chrx/chrX_data/result/${i}align.log & 
done


#2.将sam文件转化为bam文件
# 少量样本时
#samtools sort -@ 8 -o bam/ERR188044_chrX.bam bam/ERR188044_chrX.sam
# 批量处理
for i in *.sam;
do 
i=${i%_align.sam*}; nohup samtools sort -@ 8 -o ${i}.bam ${i}_align.sam & 
done


#3.组装转录本并定量表达基因
# stringtie -p 8 -G chrX_data/genes/chrX.gtf -o ERR188044_chrX.gtf -l ERR188044 bam/ERR188044_chrX.bam
# 批处理, in result dir/
for i in *.bam; 
do 
i=${i%.bam*}; nohup stringtie -p 8 -G ~/RNA-seq/chrx/chrX_data/genes/chrX.gtf -o ${i}.gtf ${i}.bam & 
done


#4.将所有转录本合并 此处的mergelist.txt是自己创建的，需要包含之前output.gtf文件的路径
#stringtie --merge -p 8 -G chrX_data/genes/chrX.gtf -o stringtie_merged.gtf chrX_data/mergelist.txt
stringtie --merge -p 8 -G  ~/RNA-seq/chrx/chrX_data/genes/chrX.gtf -o stringtie_merged.gtf  ~/RNA-seq/chrx/chrX_data/mergelist.txt

#5.检测相对于注释基因组，转录本的组装情况
#gffcompare -r chrX_data/genes/chrX.gtf -G -o merged stringtie_merged.gtf
gffcompare -r ~/RNA-seq/chrx/chrX_data/genes/chrX.gtf  -G -o merged stringtie_merged.gtf

#6.重新组装转录本并估算基因表达丰度，并为ballgown创建读入文件
#stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/ERR188044/ERR188044_chrX.gtf bam/ERR188044_chrX.bam
# 批处理, in result dir/
for i in *_chrX.bam; 
do 
i=${i%_chrX.bam*}; nohup stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/${i}/${i}_chrX.gtf  ${i}_chrX.bam & 
done
