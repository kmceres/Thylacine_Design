#!/usr/bin/env Rscript
library(optparse)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="path to pangenome dir", metavar="character"),
  make_option(c("-t","--target_size"), type="numeric", default=300, 
              help="size of genome target region", metavar="numeric"),
  make_option(c("-p","--primer_size"), type="numeric", default=20, 
              help="size of each primer", metavar="numeric"),
  make_option(c("-o", "--out"), type="character", default="potential_primers.csv", 
              help="output file name [default= %default]", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file)", call.=FALSE)
}
base_dir = opt$file

# load packages
library(dplyr)
library(stringr)
library(adegenet)
library(rhierbaps)
library(seqinr)
library(hierfstat)
library(ggplot2)
library(R453Plus1Toolbox)
library(tidyr)


## Read in gene presence absence matrix
filename = paste0(base_dir, "gene_presence_absence_roary.csv")
genepa <- read.csv(filename)

## Retain strict core
strict_core <- genepa %>% filter(No..isolates == genepa$No..isolates[1])
strict_core <- strict_core %>% mutate(clean_gene = str_remove(Gene, "~.*")) %>%
  mutate(clean_gene=str_remove(clean_gene, "_.*")) %>%
  relocate(Gene,clean_gene)



stat_df <- tibble(gene = character(),
                  Hs = numeric(),
                  POS = numeric(),
                  n_snps = numeric(),
                  gene_len = numeric())
for (gene in 1:nrow(strict_core)){
  filename = paste0(base_dir, "/aligned_gene_sequences/")
  filename = paste0(filename,strict_core$Gene[gene])
  filename = paste0(filename, ".aln.fas")
  aln <-read.alignment(filename, format="fasta")
  gene_length = nchar(aln$seq[1])
  obj <- alignment2genind(aln)
  obj$pop <- rep(as.factor("1"),aln$nb)
  obj.hf = genind2hierfstat(obj)
  
  stats <- basic.stats(obj.hf, diploid = F)
  pl_stats <- stats$perloc
  div <- summary(obj)
  pl_stats$POS <- as.numeric(names(div$loc.n.all))
  out <- pl_stats %>% dplyr::select(Hs, POS) %>% 
    mutate(gene = strict_core$Gene[gene])%>%
    relocate(gene) %>%
    mutate(n_snps = nrow(pl_stats))%>% 
    mutate(gene_len = gene_length) 
  stat_df <- bind_rows(stat_df, out)
}

sum <- stat_df %>% group_by(gene) %>% filter(gene_len > 300) %>%
  summarise(mean_Hs = mean(Hs), mean_n_snps = mean(n_snps), gene_len=mean(gene_len)) %>%
  arrange(mean_Hs) %>% mutate(frac_snp = mean_n_snps / gene_len) %>% 
  mutate(median_frac_snp = quantile(frac_snp)[3]) %>%
  mutate(median_mean_Hs = quantile(mean_Hs)[3])

best_targets <- sum %>% filter(mean_Hs < median_mean_Hs & frac_snp < median_frac_snp)

best <- stat_df %>% filter(gene %in% best_targets$gene)


window_len = 300
window_df <- tibble(gene = character(), 
                    window_start = numeric(), 
                    window_end = numeric(), 
                    mean_Het = numeric(),
                    n_high_Het = numeric())

gene_list <- unique(best$gene)
for (gene in 1:length(gene_list)){
  gene_targ <- gene_list[gene]
  df <- best %>% filter(gene == gene_targ)
  hs_df <- tibble(gene=character(),window_start = numeric(), 
                  window_end = numeric(),
                  mean_Hs = numeric(), n_high_Hs = numeric())
  for (w in 1:(df$gene_len[1]-window_len)){
    next_window_start = w
    next_window_end = w + 300 - 1
    Hs_vals <- df$Hs[df$POS >= next_window_start & 
                           df$POS <= next_window_end]
    n_high <- length(Hs_vals[Hs_vals > 0.3])
    row = tibble(window_start=next_window_start,
                 window_end = next_window_end,
                 mean_Hs=mean(Hs_vals), n_high_Hs=n_high)
    hs_df <- bind_rows(hs_df, row)
  }
  hs_df <- hs_df %>% filter(mean_Hs != "NaN")
  hs_df_min <- hs_df %>% mutate(gene = gene_targ) %>% 
    mutate(min_n = min(n_high_Hs))%>%
    filter(n_high_Hs == min_n) %>% 
    mutate(mean_Het = mean_Hs) %>%
    mutate(n_high_Het = n_high_Hs) %>%
    dplyr::select(gene, window_start, window_end, mean_Het, n_high_Het)
  
  window_df <- bind_rows(window_df, hs_df_min)
}


window_stats = tibble(gene = character(), window_start = numeric(),
                      window_end = numeric(), GC = numeric(), 
                      complexity = numeric(), seq = character())
for (gene in 1:length(gene_list)){
  gene_targ <- gene_list[gene]
  file = paste0(base_dir, "/aligned_gene_sequences/")
  file = paste0(file,gene_targ)
  file = paste0(file, ".aln.fas")
  aln = read.alignment(file, "fasta")
  seq = aln$seq[1]
  window_df_seg <- window_df %>% filter(gene == gene_targ)
  for (w in 1:nrow(window_df_seg)){
    sub_st <- substr(seq, window_df_seg$window_start[w], window_df_seg$window_end[w])
    seq_str <- s2c(sub_st)
    gc <- GC(seq_str)
    seq_comp <- complexity.dust(DNAStringSet(sub_st))
    row = tibble(gene=gene_targ, window_start = window_df_seg$window_start[w],
                 window_end = window_df_seg$window_end[w], 
                 GC = gc, complexity = seq_comp,
                 seq = toupper(sub_st))
    window_stats <- bind_rows(window_stats, row)
  }
  
}

window_stats_f <- window_stats %>% filter(GC >= .40 & GC <= .60) %>% 
  filter(complexity < 10) %>% filter(str_detect(seq,"-")==F)

window_stats_ff <- window_stats_f %>% group_by(gene) %>% 
  mutate(min_complexity = min(complexity)) %>% 
  filter(complexity == min_complexity)

window_size = 20
short_window_df  = tibble(gene = character(),
                          pos_start = numeric(),
                          pos_end = numeric(),
                          seq = character(),
                          target_seq = character(),
                          n_snp = numeric(),
                          pos_of_snps = numeric(),
                          het = numeric(),
                          GC = numeric(),
                          complexity = numeric())
for (r in 1:nrow(window_stats_ff)){
  win_start = window_stats_ff$window_start[r]
  win_end = window_stats_ff$window_end[r]
  targ_gene = window_stats_ff$gene[r]
  targ_seq = window_stats_ff$seq[r]
  new_df <- best %>% filter(gene == targ_gene) %>% 
    filter(POS >= win_start & POS <= win_end) %>% 
    mutate(targ_seq = window_stats_ff$seq[r])
  win_length = nchar(new_df$targ_seq[1])
  w = 1
  while (w < (win_length - window_size)){
    new_start = w
    new_end = new_start +  window_size - 1
    translate_start = win_start + w - 1
    translate_end = translate_start + window_size - 1
    sub_st = substr(new_df$targ_seq[1], new_start, new_end)
    new_win = new_df %>% filter(POS >= translate_start & POS <= translate_end)
    n_snps_in_win = ifelse(any(c(translate_start:translate_end) %in% new_df$POS),nrow(new_win),0)
    mean_het = ifelse(any(c(translate_start:translate_end) %in% new_df$POS), mean(new_win$Hs), 0)
    pos_w_snps = paste(new_win$POS, collapse=",")
    seq_str <- s2c(sub_st)
    gc <- GC(seq_str)
    seq_comp <- complexity.dust(DNAStringSet(sub_st))
    row = tibble(gene = new_df$gene[1],
                 pos_start = translate_start,
                 pos_end = translate_end,
                 seq = toupper(sub_st),
                 target_seq = targ_seq,
                 n_snp = n_snps_in_win,
                 pos_of_snps = pos_w_snps,
                 het = mean_het,
                 GC = gc,
                 complexity = seq_comp
                 )
    short_window_df <- bind_rows(short_window_df, row)
    w = w + 1
  }
}
  
short_window_df_f <- short_window_df %>% 
  filter(GC >= 0.4 & GC <= 0.6) %>%
  filter(complexity < 5) 

write.csv(short_window_df_f, opt$out, row.names = F)

short_window_df_long <- short_window_df_f %>% 
  gather("start_end","POS", -c(gene, seq,target_seq,n_snp,het,GC,complexity))


primer_locs <- ggplot(short_window_df_long)+
  geom_line(aes(x=POS,y=n_snp,group=seq,col=complexity))+facet_wrap(~gene)

out_name = paste0(base_dir, "/primer_locations.png")
ggsave(out_name,primer_locs)
