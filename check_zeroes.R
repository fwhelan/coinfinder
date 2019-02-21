library(getopt)
library(ape)

#Get call input
spec <- matrix(c('path', 'a', 1, "character",
                 'phylogeny', 't', 1, "character"), byrow=TRUE, ncol=4)
opt <- getopt(spec)
setwd(opt$path)

#Read in tree
tree1 <- read.tree(opt$phylogeny)
'0' %in% tree1$edge.length
