library(optparse)

option_list = list(
  make_option(c("-a", "--alignment"), type="character", default=NULL, 
              help="path to alignment file", metavar="character"),
  make_option(c("-c", "--hierbaps_results"), type="character", default=NULL, 
              help="path to hierbaps_results", metavar="character"),
  make_option(c("-p", "--ploidy"), type="numeric", default=1, 
              help="ploidy, default haploid (1)", metavar="numeric"),
  make_option(c("-w", "--window"), type="numeric", default=300, 
              help="window, default 300 bp", metavar="numeric")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
print(opt)
if (is.null(opt$alignment)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (alignment)", call.=FALSE)
}
if (is.null(opt$hierbaps_results)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (hierbaps results)", call.=FALSE)
}

file = opt$alignment 
hb.out  <- read.csv(opt$hierbaps_results)
ploidy <- opt$ploidy 
window <- opt$window
# load libraries
library(rhierbaps)
library(adegenet)
library(hierfstat)
library(ggplot2)
library(dplyr)
library(stringr)

snp.matrix <- load_fasta(file)
colnames(snp.matrix) <- paste("SNP",1:ncol(snp.matrix),sep="_")
factor_l=rep(NA,ncol(snp.matrix))
for (c in 1:ncol(snp.matrix)){
  lev=levels(factor(snp.matrix[,c]))
  factor_l[c] = ifelse("-" %in% lev & length(lev)==2, "del","snp")
}

filt_snps <- toupper(snp.matrix[,factor_l=="snp"])
filt_snps[filt_snps=="-"] <- NA

obj <- df2genind(filt_snps, pop=hb.out$level.1, ploidy = ploidy)
obj.hf = genind2hierfstat(obj)

is_diploid = ifelse(ploidy == 2,T,F)

genet.dist(obj.hf, method='Fst', diploid=is_diploid)
stats <- basic.stats(obj.hf, diploid = is_diploid)
per_loc_stats <- stats$perloc
best_stats <- per_loc_stats %>% filter(Fst > 0.5)

per_loc_stats$POS <- as.numeric(str_remove(rownames(per_loc_stats), "SNP_"))
fst <- ggplot(per_loc_stats)+geom_point(aes(x=POS,y=Fst))+theme_bw()

window_len = window
fst_df <- tibble(window_start = numeric(), window_end = numeric(),
                 mean_Fst = numeric(), n_high_Fst = numeric())
for (w in 1:max(per_loc_stats$POS)){
  next_window_start = w
  next_window_end = w + window_len - 1
  fst_vals <- per_loc_stats$Fst[per_loc_stats$POS >= next_window_start & 
                                  per_loc_stats$POS <= next_window_end]
  n_high <- length(fst_vals[fst_vals > 0.5])
  row = tibble(window_start=next_window_start,
               window_end = next_window_end,
               mean_Fst=mean(fst_vals), n_high_Fst=n_high)
  fst_df <- bind_rows(fst_df, row)
}

fst_df_max <- fst_df %>% mutate(max_n = max(n_high_Fst))%>%
  filter(n_high_Fst == max_n)

fst_w <- ggplot()+geom_point(data=per_loc_stats, aes(x=POS,y=Fst))+
  geom_vline(xintercept = fst_df_max$window_start[10])+
  geom_vline(xintercept = fst_df_max$window_end[10])+theme_bw()


ggsave("fst_window.png",fst_w,width=8,height=2)





