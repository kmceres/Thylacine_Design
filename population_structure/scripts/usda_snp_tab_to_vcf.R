library(dplyr)
library(tidyr)
read_snp_table <- function(path_to_csv, is_csv=T){
  if (is_csv == T){
    tab <- read.csv(path_to_csv)
  } else{
    tab <- readxl::read_excel(path_to_csv)
  }
  tab_t <- t(tab)
  colnames(tab_t) <- tab_t[1,]
  tab_t <- tab_t[-1,]
  return(tab_t)
}

get_snp_df<- function(tab_t){
  snpdf = as.data.frame(tab_t)
  snpdf$id = noquote(rownames(snpdf))
  snpdf = snpdf %>% separate(id, into = c("NC", "ref","4", "pos")) %>% select(-c(NC, ref, "4"))
  snpdf$pos = as.numeric(snpdf$pos)
  snpdf = snpdf %>% arrange(pos) %>% relocate(pos,MQ,annotations)#relocate(pos,`map-quality`,annotations)
  return(snpdf)
}

df_to_vcf <- function(snpdf,haploid=T){
  npos <- nrow(snpdf)
  vcf_t <- tibble(`#CHROM`=rep("NC_002945",npos), POS=rep(NA,npos), ID=rep(".",npos), 
                  REF=rep(NA,npos), ALT=rep(NA,npos), QUAL=rep(".",npos),	
                  FILTER=rep(".",npos), INFO=rep(".",npos), FORMAT= rep("GT",npos),
                  MQ=rep(NA,npos))
  for(r in 1:npos){
    vcf_t$POS[r]=snpdf$pos[r]
    vcf_t$REF[r] = snpdf$root[r]#snpdf$reference_call[r]
    t=as.data.frame(table(t(snpdf[r,5:ncol(snpdf)])))
    vcf_t$ALT[r] = paste(t$Var1[t$Var1!=vcf_t$REF[r]], collapse=",")
    vcf_t$MQ[r] = paste("MQ",snpdf$`map-quality`[r],sep="=")
  }
  snp_data <- snpdf[,5:ncol(snpdf)]
  if (haploid == F){
    for(r in 1:nrow(snp_data)){
      for(c in 1:ncol(snp_data)){
        snp_data[r,c]=ifelse(snp_data[r,c]==vcf_t$REF[r],"0|0","1|1")
                            #   paste(snp_data[r,c],snp_data[r,c],sep="|")
      }
    }
  } else{
    for(r in 1:nrow(snp_data)){
      for(c in 1:ncol(snp_data)){
        snp_data[r,c]=ifelse(snp_data[r,c]==vcf_t$REF[r],"0","1")
      }
    }
  }
  vcf_t <- bind_cols(vcf_t, snp_data)
  
  return(vcf_t)
}

write_out_vcf <- function(vcf_t,path_to_file,haploid=T){
  sink(path_to_file)
  cat("##fileformat=VCFv4.2")
  cat("\n")
  cat('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype",MQ="MappingQuality">')
  cat("\n")
  writeLines(paste(colnames(vcf_t),collapse ="\t"))
  for(i in 1:nrow(vcf_t)){
    writeLines(paste(vcf_t[i,],collapse="\t"))
  }
  sink()
}



snp_tab_to_fasta <- function(snp_tab, path_to_file){
  sink(path_to_file)
  for (i in 1:nrow(snp_tab)){
    name = noquote(paste(">",snp_tab[i,1],sep=""))
    cat(name,sep="\n")
    ln <- paste0(unlist(as.character(snp_tab[i,2:ncol(snp_tab)])),sep="", collapse="")
    cat(ln)
    cat("\n")
  }
  sink()
}


snp_tab_to_seq_table <- function(snp_tab, path_to_file){
  sink(path_to_file)
  for (i in 1:(nrow(snp_tab)-2)){
    name = as.character(snp_tab[i,1])
    cat(name)
    cat("\t")
    ln <- paste0(unlist(as.character(snp_tab[i,2:ncol(snp_tab)])),sep="", collapse="")
    cat(ln)
    cat("\n")
  }
  sink()
}

fasta_to_seq_table <- function(fasta, path_to_file){
  sink(path_to_file)
  for (i in 1:(nrow(fasta)-2)){
    name = as.character(fasta[i,1])
    cat(name)
    cat("\t")
    ln <- paste0(unlist(as.character(fasta[i,2:ncol(fasta)])),sep="", collapse="")
    cat(ln)
    cat("\n")
  }
  sink()
}
