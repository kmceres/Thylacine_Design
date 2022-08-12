##### R scripts for calculations
#https://www.nature.com/articles/s41598-019-55253-0

# reading the aligned sequence data:
#  the data have to be formatted in tab-separated text with two colums,
#  (name of sequence) \t (aligned sequence)
##  TableS2.txt: human, TableS3.txt: lion, TableS4.txt: bacteria

snps_to_boolean_vector <- function(path_to_file){
  sites <- read.table(path_to_file)
  sites<-as.matrix(sites)
  dim(sites)

  ### finding the size of data
  n_sample <- dim(sites)[1]
  n_seq <- nchar(sites[2,2])

  ### translation of the sequence to bollean vectors
  bool <- array(0, dim=c(n_sample, 5*n_seq))

  colnames(bool) <- c(paste("A_", 1:n_seq, sep=""),paste("T_", 1:n_seq, sep=""),paste("G_", 1:n_seq, sep=""),paste("C_", 1:n_seq, sep=""),paste("N_", 1:n_seq, sep=""))
  rownames(bool) <- sites[ ,1]
  
  
  for (s in 1:n_sample){
    se <- sites[s, 2]
    se <- tolower(se)
    
    for (le in  1:n_seq){
      base <- substr(se, le,le)
      
      if(base =="a") {
        bool[s, le] <-1
      } else {
        
        if(base =="t") {
          bool[s, le+n_seq] <-1
        } else {
          
          if(base =="g") {
            bool[s, le+n_seq*2] <-1
          } else {
            
            if(base =="c") {
              bool[s, le+n_seq*3] <-1
            } else {
              
              bool[s, le+n_seq*4] <-1
            } 
            }}}
    }}
  return(bool)  
}

prepare_pca <- function(bool){
  ## centering : the center can be replaced to certain group
  center<- apply(bool, 2, mean)
  diffs<-sweep(bool, 2, center)
  diffs <- diffs/(2^0.5)   
  # compensating the doubled counts in Euclidean distance metrics
  # checking distribution of the distances
  dist<- (apply(diffs^2, 1, sum))^0.5
  qqnorm(dist)
  hist(dist)
  apply(bool, 1, sum) 
  return(diffs)
}

 # here you can verify the translation
# they should show identical values same as n_seq

############ PCA 





# 
# 
# 
# ### PCA core
# res_svd <- svd(diffs)  #
# str(res_svd)
# Left <- res_svd$u		# the left singular vector
# Right <- res_svd$v		# the right singular vector
# sqL <- diag(res_svd$d)		# diagonal matrix of the singular values
# 
# ### calculatinf of pc's 
# sPC_nuc  	<-	 Right %*% sqL / (n_sample^0.5)
# sPC_sample	 <-	 Left %*% sqL/ (n_seq^0.5)
# 
# rownames(sPC_nuc)<- colnames(bool) 
# rownames(sPC_sample)<- rownames(bool) 
# 
