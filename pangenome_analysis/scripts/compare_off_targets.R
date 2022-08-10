library(optparse)

option_list = list(
  make_option(c("-f", "--file"), type="character", default="potential_primers.csv", 
              help="path to pangenome dir", metavar="character"),
  make_option(c("-t","--off_target"), type="character", default=NULL, 
              help="path to off target genome directory", metavar="numeric"),
  make_option(c("-s","--max_subs"), type="numeric", default=5, 
              help="path to off target genome directory", metavar="numeric"),
  make_option(c("-o", "--out"), type="character", default="primers_no_off_target_matches.csv", 
              help="output file name [default= %default]", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.null(opt$off_target)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (off target directory)", call.=FALSE)
}


library(dplyr)
library(abdiv)
library(seqinr)
primers = read.csv(opt$file)
off_target_dir= opt$off_target 
max_subs = opt$max_subs

file_list <- list.files(path=off_target_dir)

new_cols = matrix(NA, nrow=nrow(primers), ncol=length(file_list))
for (f in 1:length(file_list)){
  file = paste0(off_target_dir, file_list[f])
  off <- read.fasta(file,seqonly = T)
  for(row in 1:nrow(primers)){
    dif = agrep(primers$seq[row],off, max.distance = c(substitutions=max_subs, insertions = 0, deletions = 0))
    new_cols[row,f] = ifelse(identical(dif, integer(0)),0,1) 
  }
}

off_df <- as.data.frame(new_cols)
off_df <- off_df %>% mutate(n_fail = rowSums(off_df))
out <- primers[off_df$n_fail == 0,]
write.csv(out, opt$out ,row.names = F)
