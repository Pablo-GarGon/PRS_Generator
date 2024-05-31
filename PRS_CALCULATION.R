args=commandArgs(trailingOnly = TRUE)

library("reshape2")

file_betas <- args[1]
file_info <-  args[2]
file_genotypes <-  args[3]
output <- args[4]
output.dir <- args[5]
filter.r2.maf <- as.logical(args[6])
r2.threshold <- as.numeric(args[7])
maf.threshold <- as.numeric(args[8])
eadb.correction <- as.logical(args[9])

# Set work dir to output dir
setwd(output.dir)

# Load summary statistics:
betas<-read.table(file_betas, h=T, sep="\t", stringsAsFactors=F)
betas$Beta <- as.numeric(betas$Beta)

# load genotypes
genotypes<-read.table(file_genotypes, h=T, stringsAsFactors=F)
names(genotypes)<-sub("X","chr",names(genotypes))

# imputation quality
info<-read.table(file_info, h=T, stringsAsFactors=F)

# Create ID variable (gsub : for . and add ref allele after _)
info<-cbind(info,colsplit(info$SNP, ":", names = c("CHR", "POS", "REF", "ALT"))[,c(1,2)])

names(info)<-c("ID_genotype_file","REF_allele","ALT_allele", "Freq_ALT_allele","MAF","callrate","Rsq","Genotyped","LooRsq", "EmpR","EmpRsq","Dose0","Dose1", "CHR","POS")
info$ID_dosage_file<-paste0(gsub(":", ".", info$ID_genotype_file), "_", info$REF_allele, sep="")

# merge the infofile and the betas 

info<-merge(info,betas, by.x="ID_genotype_file", by.y="ID")
print(paste(dim(info)[1],"SNPs present in genotype file"))

row.names(info) <- info$ID_genotype_file


# ID1 <- ID_genotype_file

# Flip SNPs when neccesary:

info$FLIP_SNP<-NA
info$FLIP_SNP[which(info$REF_allele==info$A1_summary & info$ALT_allele==info$A2_summary)]<-0
info$FLIP_SNP[which(info$REF_allele==info$A2_summary & info$ALT_allele==info$A1_summary)]<-1

# start building the PRS


data <- genotypes
data[,"Polygenic_score"]<-0

bad.SNPs <- info[info$Rsq < r2.threshold | info$MAF < maf.threshold,]


print(paste(dim(bad.SNPs)[1]," SNPs should be excluded in PRS calculation due to low imputation quality ( R2 < ",r2.threshold,") or MAF ( MAF<",maf.threshold,"):"), sep="")
if (dim(bad.SNPs)[1] < 20 & dim(bad.SNPs)[1] > 0){
  print(bad.SNPs$ID_genotype_file)
}else{print("SNPs not shown (n>20)")}


#Eliminar los bad SNPs si esta variable es true
if(filter.r2.maf==T){
  data <- data[,!names(data) %in% bad.SNPs$ID_dosage_file]
  info <- info[!info$ID_dosage_file %in% bad.SNPs$ID_dosage_file,]
  
}


# Flip SNP dosages based on previous step and calculate PRS:

for (i in 1:dim(info)[1]){
  snp<-info$ID_dosage_file[i]
  beta<-info$Beta[i]
  info$Beta_freq[i]<-sum(data[,snp])/(dim(data)[1]*2)
  if(info$FLIP_SNP[i] == 0){data[,paste0("Polygenic_score")]<-data[,paste0("Polygenic_score")]+beta*data[,snp]}
  if(info$FLIP_SNP[i] == 1){data[,paste0("Polygenic_score")]<-data[,paste0("Polygenic_score")]+beta*(2-data[,snp])}
  
}

print(paste("PRS calculated for", dim(info)[1], "SNPs"))


if(eadb.correction==T){
  
  correction <- dim(betas)[1]/sum(betas$Beta)
  data$Polygenic_score <- data$Polygenic_score*correction
  print("PRS corrected by dim(betas)[1]/sum(betas$Beta)")
  
}

output.prs <- paste("PRS_", output, ".txt", sep="")
output.info <- paste("REPORT_PRS_", output, ".txt", sep="")
output.histogram <- paste("Histogram_PRS_", output, ".png", sep="")

write.table(data, output.prs, row.names=F, quote=F, sep="\t")
write.table(info, output.info, col.names=T, row.names=F, quote=F,sep="\t")

png(output.histogram)
hist(data$Polygenic_score, main="PRS distribution", xlab="PRS")
dev.off()



