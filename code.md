##### **Index**

Reference=GCF\_000001405.39



bowtie2-build Reference Reference.fasta



##### Mapping, Alignment, and Coverage





bowtie2 --maxins 800 -t -x Reference -p 8 --fr -1 Seq1.fastq.gz -2 Seq2.fastq.gz | samtools view -Sb -@8 | samtools sort --@8 > Doc.sorted.bam 



samtools coverage Doc.sorted.bam > coverage\_report.txt 



samtools view Doc.sorted.bam | awk '{print $3}' | sort | uniq -c | awk '{print $2 \\t $1}' > gene\_counts.txt 



##### Counting table



genera\_tabla\_conteos<-function(fnPath,fnExt){

&nbsp; fnPattern=paste(fnPath, fnExt,"$",sep="",collapse="")

&nbsp; fnLista\_archivos=dir(fnPath, fnExt)

&nbsp; print(fnLista\_archivos)

&nbsp; fnNomExperimentos = sub(fnExt,'', fnLista\_archivos)

&nbsp; fnTablaConteos = data.frame(matrix(vector(), 0, length(fnLista\_archivos)))

&nbsp; names(fnTablaConteos)<-fnNomExperimentos

&nbsp; for(i in fnNomExperimentos)

&nbsp; {

&nbsp;   print(paste(fnPath, i, fnExt, sep = ""))

&nbsp;   fnColExp<-read.delim(paste(fnPath, i, fnExt, sep = ""), header = T,row.names = 1) 

&nbsp;   print(head(fnColExp))

&nbsp;   fnTablaConteos\[row.names(fnColExp),i]<-fnColExp\[,1]

&nbsp; }

&nbsp; fnTablaConteos\[is.na(fnTablaConteos)]=0

&nbsp; fnTablaConteos<-fnTablaConteos\[rowSums(fnTablaConteos)>0,]

&nbsp; return(fnTablaConteos)

}

myPath<-"D:/AZ/subset/tabla\_de\_conteos/"

myExt<-".txt"

t\_d\_c<-genera\_tabla\_conteos(myPath,myExt)

t\_d\_c

write.table(t\_d\_c,"D:/AZ/subset/tabla\_de\_conteos/BX\_conteos.txt", sep="\\t", quote=F)



names(t\_d\_c)



