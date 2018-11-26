library(caper)
library(phytools)
library(getopt)
library(future)

#Get call input
spec <- matrix(c('path', 'p', 1, "character",
                 'phylogeny', 't', 1, "character",
                 'gene_pa', 'g', 1, "character",
                 'cores', 'c', 1, "integer",
                 'output', 'o', 1, "character"), byrow=TRUE, ncol=4)
opt <- getopt(spec)
setwd(opt$path)

#Read in
outstr <- paste(opt$output, "_nodes_in.csv", sep="")
genes  <- read.csv(outstr) #"coincident_nodes_in.csv")
genes[,length(genes)] <- NULL #remove the last comma in csv file
genepa <- read.csv(opt$gene_pa, header=T, row.names=1)
tree <- read.tree(opt$phylogeny)
#Make annot table
genepa[,1:14] <- NULL
rownames(genepa) <- gsub(" ", "", rownames(genepa))
annot <- genepa[(rownames(genepa) %in% names(genes)),]
genepa <- NULL
annot <- t(annot)
annot[annot!=""] <- "1"
annot[annot==""] <- "0"
annot <- as.data.frame(annot)
annot$Id <- rownames(annot)
#Double check that each column has 2 states
remove <- vector()
for(a in 1:length(colnames(annot))) {
  if (length(table(annot[a])) == 1) {
    print("Encountered the following gene that only exists in one state in the data.")
    print(colnames(annot)[a])
    print("I'm removing it, but you should be cautious as to what this means...")
    remove <- append(remove, a)
  }
}
#Remove any columns with only one state
if (length(remove)>0) {
  for(a in 1:length(remove)) {
    annot[remove[a]] <- NULL
  }
}
#Root phylogeny for input into comparative.data
treeRt <- midpoint.root(tree)
treeRt <- makeLabel(treeRt)
#Make comparative data object
dataset <- comparative.data(phy = treeRt,
                            data = annot,
                            names.col = Id,
                            vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)
#Estimate D (Fritz & Purvis 2010) for all columns in annot
#Detect available cores and decrease the number of parallel runs to that if its < opt$cores
availcores <- availableCores()
cores <- 1
if (availcores < opt$cores) {
	cores <- availcores
} else {
	cores <- opt$cores
}
parallelCluster <- parallel::makeCluster(cores)
mkWorker <- function(dataset) {
  library(caper)
  #Make sure each value is passed
  force(dataset)
  #Define function
  calcD <- function(dataset, binvarry) {
    if(binvarry != "Id") {
      result <- eval(parse(text=paste("caper::phylo.d(data=dataset, binvar=",binvarry,", permut=1000)", sep="")))
    }
  }
  #Define & return worker function
  worker <- function(binvarry) {
    calcD(dataset, binvarry)
  }
  return(worker)
}
results <- parallel::parLapply(parallelCluster, colnames(annot), mkWorker(dataset))
#Update genes to be a data.frame
genes2 <- as.data.frame(cbind(names(genes), rep(NA,length(genes))), stringsAsFactors=FALSE)
colnames(genes2) <- c("ID", "Result")
rownames(genes2) <- genes2$ID
genes2$ID <- NULL
for(a in 1:length(results)) {
  if(!is.null(results[[a]])) {
      genes2[results[[a]]$binvar,1] <- results[[a]]$DEstimate
  }
}
genes2$ID <- rownames(genes2)
genes2 <- genes2[,c(2,1)]
outstr <- paste(opt$output, "_nodes.csv", sep="")
write.table(genes2, outstr, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE) #"coincident_nodes.csv"
