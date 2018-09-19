library(caper)
library(phytools)

setwd("~/pCloud\ Drive/pCloud\ Sync/postdoc-research/18.01-Networks-pangenomes/18.01-PES/roary_out_s_e_n/")
#Read in
genes  <- read.csv("coincident_nodes_in.csv")
genes[,length(genes)] <- NULL #remove the last comma in csv file
genepa <- read.csv("gene_presence_absence.csv", header=T, row.names=1) #TODO
tree <- read.tree("PES-core_gene_aln.newick") #TODO
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
parallelCluster <- parallel::makeCluster(parallel::detectCores()) #**TODO update this to a user defined set of cores
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
setwd("~/Documents/software/coinfinder/coinfinder-plus/coinfinder-plus/")
write.table(genes2, "coincident_nodes.csv", sep=",", col.names=FALSE, quote=FALSE)
