# projects


the  tutorial version of hisat  scripts uesage:
./rnaseq_pipeline.sh  outputdir  
when you run this script ,you may count the problem 
[2019-04-24 16:27:41] #> START:  ./rnaseq_pipeline.sh outputdir
[2019-04-24 16:27:41] Processing sample: *
[2019-04-24 16:27:41]    * Alignment of reads to genome (HISAT2)

the problem was that the data is can not found 
the right path is BASEDIR="/media/gsadmin/vd2/mysql/chrX_data/chrX_data"
