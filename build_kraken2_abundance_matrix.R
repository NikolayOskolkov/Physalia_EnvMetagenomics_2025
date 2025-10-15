#This is a script for combining Kraken2 outputs from multiple samples and producing species abundance matrix.
#Run this script as: Rscipt build_kraken2_abundance_matrix.R input_dir output_dir n_unique_kmers n_tax_reads
# For example: Rscipt build_kraken2_abundance_matrix.R kraken_input_dir results 1000 200

args = commandArgs(trailingOnly=TRUE)
input_dir<-as.character(args[1])
output_dir<-as.character(args[2])
n_unique_kmers<-as.integer(args[3])
n_reads<-as.integer(args[4])

#Reading and merging kraken outputs from multiple samples
kraken_outputs<-list.files(path=input_dir,pattern="^kraken")
df<-list()
for(i in 1:length(kraken_outputs))
{
 df[[i]]<-read.delim(paste0(input_dir,"/",kraken_outputs[i]),comment.char="#",header=FALSE)
 colnames(df[[i]])<-c("pers_reads","reads","taxReads","n_minimizers","n_distinct_minimizers","rank","taxID","taxName")
 if(dim(df[[i]])[1]!=0)
 {
  df[[i]]$SAMPLE<-kraken_outputs[i]
  df[[i]]<-na.omit(df[[i]])
  df[[i]]<-df[[i]][grepl("S",as.character(df[[i]]$rank)) & df[[i]]$reads>n_reads & df[[i]]$n_distinct_minimizers>n_unique_kmers,]
 }
 else{next}
}
merged<-Reduce(rbind,df)
merged$taxName<-trimws(as.character(merged$taxName))
merged<-merged[order(-merged$pers_reads),]
print(head(merged,20))

unique_samples<-unique(merged$SAMPLE)
unique_species<-unique(merged$taxName)
unique_taxids<-unique(merged$taxID)

#Computing abundance matrix
abundance_matrix<-matrix(NA,ncol=length(unique_samples),nrow=length(unique_species))
for(i in 1:length(unique_species))
{
 for(j in 1:length(unique_samples))
 {
  if(length(merged[merged$taxName==unique_species[i] & merged$SAMPLE==unique_samples[j],]$reads)==0)
  {
   abundance_matrix[i,j]<-0
  }
  else
  {
   abundance_matrix[i,j]<-merged[merged$taxName==unique_species[i] & merged$SAMPLE==unique_samples[j],]$reads
  }
 }
}
rownames(abundance_matrix)<-unique_species
colnames(abundance_matrix)<-unique_samples
abundance_matrix<-abundance_matrix[order(rownames(abundance_matrix)),]
print(head(abundance_matrix,20))

write.table(abundance_matrix,file=paste0(output_dir,"/kraken2_abundance_matrix.txt"),col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")

