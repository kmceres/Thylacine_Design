library(optparse)

option_list = list(
  make_option(c("-a", "--alignment"), type="character", default=NULL, 
              help="path to alignment file", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.null(opt$alignment)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (alignment)", call.=FALSE)
}
library(seqinr)
library(factoextra)
library(ggplot2)
library(abdiv)
library(ggpubr)

# read in PCA prep scripts
source("scripts/snps_to_boolean_pca.R")
source("scripts/usda_snp_tab_to_vcf.R")

file=opt$alignment
# read in alignment
snps<- seqinr::read.fasta(file)

# prepare data for PCA
snps_table <- tibble(id = rep(NA, length(snps)), dna=rep(NA,length(snps)))
for (s in 1:length(snps)){
  seq = snps[[s]]
  snps_table$id[s] = attr(seq, "name")
  snps_table$dna[s] = paste(unlist(as.character(toupper(seq[1:length(seq)]))),collapse="")
}
write.table(snps_table, "data/temp.snps.txt",sep="\t", row.names = F, col.names = F, quote = F)

b <- snps_to_boolean_vector("data/temp.snps.txt")
d <- prepare_pca(b)

# run pca
pca <- prcomp(d, center = F, scale. = F) # prepare_pca already centers and scales
var <- fviz_eig(pca) # see which dimesnions explain the most variance
eig <- get_eig(pca) # eigenvalues and variance %

# prepare data for plotting
pcs <- as.data.frame(pca$x)
pcsdf <- pcs
# optional, use kmeans to get optimal # clusters
opt_n_clusters <- fviz_nbclust(pcs,kmeans,k.max=20) # get optimal number of clusters = 10

# plot PCs
p12 <- ggplot(pcsdf)+geom_point(aes(x=PC1,y=PC2))
p13 <- ggplot(pcsdf)+geom_point(aes(x=PC1,y=PC3))
p23 <- ggplot(pcsdf)+geom_point(aes(x=PC2,y=PC3))
pc_plots <- ggarrange(p12,p13,p23, nrow=1)
ggsave("PC_plot.png",width=8,height=3)
