library(optparse)

option_list = list(
  make_option(c("-a", "--alignment"), type="character", default=NULL, 
              help="path to alignment file", metavar="character"),
  make_option(c("-t", "--tree"), type="character", default=NULL, 
              help="path to tree file", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.null(opt$alignment)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (alignment)", call.=FALSE)
}
if (is.null(opt$tree)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (tree)", call.=FALSE)
}

library(rhierbaps)
library(ape)
library(ggtree)
library(hierfstat)
library(phytools)
library(ggplot2)

file=opt$alignment
tree=opt$tree

# phylo clustering
snp.matrix <- load_fasta(file)
hb.results <- hierBAPS(snp.matrix, max.depth = 3, n.pops = 20, quiet = TRUE)
write.csv(hb.results$partition.df, "hierbaps_out.csv",row.names = F)
head(hb.results$partition.df)
tree <-  phytools::read.newick(tree)
gg <- ggtree(tree, layout = "circular")
gg <- gg %<+% hb.results$partition.df
gg <- gg + geom_tippoint(aes(color = factor(`level 1`))) +geom_tiplab(size=0.2)
gg <- gg + theme(legend.position = "right")
out <- gg + labs(color="Group")

ggsave("hb_tree.png",gg, height=10,width=10)

