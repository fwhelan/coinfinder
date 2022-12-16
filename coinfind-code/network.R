library(phytools)
library(igraph)
library(dplyr)
library(cowplot)
library(data.table)
library(ggtree)
library(ggraph)
library(getopt)

#Get call input
spec <- matrix(c('path', 'p', 1, "character",
                 'phylogeny', 't', 1, "character",
                 'gene_pa', 'g', 1, "character",
                 'output', 'o', 1, "character"), byrow=TRUE, ncol=4)
opt <- getopt(spec)
setwd(opt$path)

#Read in
nodstr <- paste(opt$output, "_nodes.tsv", sep="")
nodes  <- read.table(nodstr, header=TRUE, sep="\t", quote=NULL)
colnames(nodes) <- c("alphas", "D")
nodes$alphas <- make.names(nodes$alphas)
edgstr <- paste(opt$output, "_edges.tsv", sep="")
edges <- read.table(edgstr, header=TRUE, sep="\t", quote=NULL)
colnames(edges) <- c("alpha1", "alpha2", "p")
edges$alpha1 <- make.names(edges$alpha1)
edges$alpha2 <- make.names(edges$alpha2)
edges$p <- as.numeric(as.character(edges$p))
#genepa <- read.csv(opt$gene_pa, header=T, row.names=1)
tree <- read.newick(opt$phylogeny)
#Remove any nodes that do not have a value for D and edges which don't have a p-value
nodes <- nodes[complete.cases(nodes),]
edges <- edges[complete.cases(edges),]
edges <- subset(edges, (edges$alpha1 %in% nodes$alphas == TRUE) & (edges$alpha2 %in% nodes$alphas == TRUE))

#Define CCs and order by D-value
g <- graph_from_data_frame(d = edges, directed = FALSE)
CCs <- data.frame(components(g)$membership)
colnames(CCs) <- c("CC")
CCs$alphas <- rownames(CCs)
#Remove any spaces in the alpha names
CCs$alphas <- gsub(" ", "", CCs$alphas)
#Add CC information to nodes table
nodes <- merge(nodes,CCs, by="alphas")
#Reorder by D, then by CC
ord    <- nodes %>% arrange(D) %>% mutate(ord=order(D, decreasing=TRUE))
ord2   <- ord %>% group_by(CC) %>% summarize(min_ord=min(ord))
ord.CC <- merge(ord,ord2)
ord.CC <- ord.CC %>% arrange(min_ord,ord)
#Colour by connectedness
ord.CC.cnts <- count(ord.CC, CC)
tmp <- merge(ord.CC, ord.CC.cnts)
tmp <- tmp %>% arrange(min_ord)
tmp2 <- setDT(tmp)[, unique(n), by = CC]
CC.size <- tmp2$V1
#Define colour array based on the size of the components
colour.array <- c("#7bb4ff","#7c0051","#76c3ff","#001106","#01ca10","#e56800","#0005a7","#0a3700","#0245ec","#e3db00","#671a00","#00bf72","#ff9dd6","#f386ff","#d57600","#20000a","#ace967","#840019","#0094f7","#006695","#ff5bf7","#b9007e","#ff89c7","#360093","#001733","#005a69","#d5dca2","#61cb00","#00338b","#26a600","#fdd268","#ffbdca","#b9e590","#ff711b","#8ce9c7","#331300","#a40053","#02d6cf","#640015","#c90039","#fecdb7","#230025","#008c51","#ea21dd","#7a63ff","#b43b00","#0068a7","#01cf41","#f6d633","#93a200","#61f0b9","#ffbf91","#990015","#c4a8ff","#02d9b8","#f2007e","#ac8b00","#d0e072","#ffa368","#02b096","#02c1cc","#ff8de1","#fc0072","#dea3ff","#c37d00","#03c4f4","#ff5dad","#8990ff","#d2d9d9","#b744f6","#380006","#5a3700","#002a23","#392ad1","#238100","#ff673c","#9bddff","#b2e3c0","#89ef68","#00726c","#004376","#5d9fff","#210052","#009e40","#ffac48","#9d1fd7","#d7d5ef","#7cf22c","#b7ddef","#a8a7ff","#008634","#8c3300","#8400a7","#783900","#0099b3","#ff3628","#386700","#003c45","#7f7100","#001f5c","#b4a200","#272200","#ff9ea3","#ff325c","#e19100","#c80024","#ff6a61","#3c5300","#b9d400","#fc2109","#3a3b00","#024fd6","#5d7b00","#a90009","#006027","#6ab200","#016989","#4f005e","#016dbf","#29000b","#002700","#659700","#ffc545","#d7d0ff","#910085","#0064db","#ffa9c8","#dc0034","#620028","#ff5664","#006842","#441500","#ff7588","#9ce5da","#ff55cd","#c27fff","#004f31","#e2da91","#7bb4ff","#7c0051","#76c3ff","#001106","#01ca10","#e56800","#0005a7","#0a3700","#0245ec","#e3db00","#671a00","#00bf72","#ff9dd6","#f386ff","#d57600","#20000a","#ace967","#840019","#0094f7","#006695","#ff5bf7","#b9007e","#ff89c7","#360093","#001733","#005a69","#d5dca2","#61cb00","#00338b","#26a600","#fdd268","#ffbdca","#b9e590","#ff711b","#8ce9c7","#331300","#a40053","#02d6cf","#640015","#c90039","#fecdb7","#230025","#008c51","#ea21dd","#7a63ff","#b43b00","#0068a7","#01cf41","#f6d633","#93a200","#61f0b9","#ffbf91","#990015","#c4a8ff","#02d9b8","#f2007e","#ac8b00","#d0e072","#ffa368","#02b096","#02c1cc","#ff8de1","#fc0072","#dea3ff","#c37d00","#03c4f4","#ff5dad","#8990ff","#d2d9d9","#b744f6","#380006","#5a3700","#002a23","#392ad1","#238100","#ff673c","#9bddff","#b2e3c0","#89ef68","#00726c","#004376","#5d9fff","#210052","#009e40","#ffac48","#9d1fd7","#d7d5ef","#7cf22c","#b7ddef","#a8a7ff","#008634","#8c3300","#8400a7","#783900","#0099b3","#ff3628","#386700","#003c45","#7f7100","#001f5c","#b4a200","#272200","#ff9ea3","#ff325c","#e19100","#c80024","#ff6a61","#3c5300","#b9d400","#fc2109","#3a3b00","#024fd6","#5d7b00","#a90009","#006027","#6ab200","#016989","#4f005e","#016dbf","#29000b","#002700","#659700","#ffc545","#d7d0ff","#910085","#0064db","#ffa9c8","#dc0034","#620028","#ff5664","#006842","#441500","#ff7588","#9ce5da","#ff55cd","#c27fff","#004f31","#e2da91");#,"#7bb4ff");
node.colour <- NULL
for(i in c(1:length(CC.size))){
  #if(i > length(colour.array)) {
  #  print("Error: Add more colours to colour.array")
  #  quit()
  #}
  node.colour <- c(node.colour,rep(colour.array[(i%%length(colour.array))+1],times=CC.size[i]))
}
node.order <- unique(ord.CC$alphas)
names(node.colour) <- node.order
node.order <- factor(node.order, levels=node.order)

#Write CCs to coincident_components.csv
CC_out <- data.frame(id = character(0), stringsAsFactors = FALSE)
#CCs$alpha <- make.names(CCs$alpha)
#for(i in c(1:nrow(CCs))) {
#  if (is.na(CC_out[CCs$CC[i],1])) {
#    CC_out[CCs$CC[i],1] <- CCs$alphas[i]
#  } else {
#    CC_out[CCs$CC[i],1] <- paste(CC_out[CCs$CC[i],1], ",", CCs$alphas[i], sep="")
#  }
#}
ord.CC$alphas <- make.names(ord.CC$alphas)
for(i in c(1:nrow(ord.CC))) {
	  curCC <- which(tmp2$CC == ord.CC$CC[i])
  if (is.na(CC_out[curCC,1])) {
	      CC_out[curCC,1] <- ord.CC$alphas[i]
    } else {
	        CC_out[curCC,1] <- paste(CC_out[curCC,1], ",", ord.CC$alphas[i], sep="")
      }
}
outstr <- paste(opt$output, "_components.tsv", sep="")
write.table(CC_out, file=outstr, sep="\t", quote=FALSE, row.names=TRUE, col.names=FALSE)

#Create annot
#genepa[,1:14] <- NULL
#rownames(genepa) <- gsub(" ", "", rownames(genepa))
#annot <- genepa[(rownames(genepa) %in% as.character(names(node.colour))),]
#genepa <- NULL
con = file(opt$gene_pa, "r")
header=readLines(con, n=1) #read in line
header.sp = strsplit(header,"\",\"") #split line into list
header.sp[[1]][1] = gsub("\"", "", header.sp[[1]][1]) #remove first \" from line
header.sp[[1]][length(header.sp[[1]])] = gsub("\"", "", header.sp[[1]][length(header.sp[[1]])]) #and last
header.sp[[1]] <- make.names(header.sp[[1]])
annot <- matrix(ncol=length(header.sp[[1]]))
colnames(annot) <- header.sp[[1]]
flag <- 1
while(TRUE) {
  line=readLines(con, n=1) #read in line
  if (length(line)==0) {
    break #file done
  }
  line.sp = strsplit(line,"\",\"") #split line into list
  line.sp[[1]][1] = gsub("\"", "", line.sp[[1]][1]) #remove first \" from line
  line.sp[[1]][length(line.sp[[1]])] = gsub("\"", "", line.sp[[1]][length(line.sp[[1]])]) #and last
  line.sp[[1]] <- make.names(line.sp[[1]])
  if (line.sp[[1]][1] %in% names(node.colour)) { #its a gene of interest, keep around
    if (flag == 1) {
      annot[1,] <- as.character(line.sp[[1]])
      flag <- 0
    } else {
      annot <- rbind(annot, as.character(line.sp[[1]]))#, stringsAsFactors=FALSE)
    }
  }
}
close(con)
annot <- as.data.frame(annot)
rownames(annot) <- annot[,1]
annot[,1:14] <- NULL

annot <- as.data.frame(t(annot), stringsAsFactors = FALSE)
annot.cp <- annot
for(a in 1:ncol(annot)) {
  annot[,a][annot[,a]!="X"] <- as.character(colnames(annot)[a])
  annot.cp[,a][annot[,a]!="X"] <- as.character(node.colour[a])
}
setcolorder(annot, as.character(names(node.colour)))
setcolorder(annot.cp, as.character(names(node.colour)))
#setcolorder(annot, as.character(node.order))

#Output annot as a csv file
write.table(annot.cp, file=paste(opt$output, "_heatmap-raw-data.csv", sep=""), sep=",", quote=FALSE, row.names=FALSE)

#Draw heatmap
heatmap.breaks <- colnames(annot)
heatmap.breaks <- factor(heatmap.breaks, levels=node.order)
#node.colour[length(node.colour)+1] <- paste0("name", "black")
node.colour  <- setNames(c(node.colour, "white"), c(names(node.colour), "name"))

tree$tip.label <- make.names(tree$tip.label)
p.tree <- ggtree(tree) + geom_tiplab(size=7, hjust=0) #hjust=0.1
arc.breaks = as.numeric(c(0, 0.005, 1))
#Split the heatmap into sections that equal roughly 500 genes each;
#however, must be sure not to split a cluster into parts so must look
#at the size of the clusters and split the input at a location close
#to multiples of 500 based on the cluster size.
CC.cumsum <- cumsum(CC.size)
countie <- 1
a <- 0
while (countie < length(colnames(annot))) {
	#for(i in seq(0,length(colnames(annot)), by=500)) {
	#j <- min(i+500, length(colnames(annot)))
	i <- countie
  	j.loc <- which.min(abs(CC.cumsum - (countie+500)))
  	j <- CC.cumsum[j.loc]
	#However, if there is an element with >500 elements heatmap creation will spin. Inforce j > i.
	if (j<i) {
		j <- CC.cumsum[j.loc+1]
	}
	#However, however, if the distance between i and j is > 1000, we won't be able to see anything
	#if (j-i > 1000) {
	#	j <- i+500
	#}
  	countie <- j+1
	p.heat <- gheatmap(p.tree, annot[i:j], offset=0, width=8, font.size=10, colnames_angle=-90, hjust=0) +
	  guides(fill=FALSE) +
	  scale_fill_manual(breaks=heatmap.breaks, values=node.colour) +
	  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
	  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
	#Create arced network
	#grp <- graph_from_data_frame(d = subset(edges, (edges$alpha1 %in% node.order[i:j] & edges$alpha2 %in% node.order[i:j])), vertices = node.order[i:j], directed = FALSE)
	#p.arc <- ggraph(grp, layout="linear") +
	#  geom_edge_arc(width=1.5, alpha = 1, curvature=-1) + #width=E(grp)$p, label=p aes(edge_colour=(1-E(grp)$p))
	#  #scale_edge_colour_gradient2(low = "#f0f0f0",
       	#  #                           mid = "gray",
        #  #                           high = "black",
        #  #                           midpoint = 0.005,
        #  #	                      #trans = "log",
        #  #	                      breaks=arc.breaks,
        #  #	    	              labels=arc.breaks
	#  #			      ) +
  	#geom_node_point(color="gray") +
  	#geom_node_text(aes(label = node.order[i:j]), angle = 90, hjust=0.8) +
  	#theme_graph() +
  	#theme(legend.position="bottom") + #bottom #none
  	#theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
	#Output
	#p.blank <- ggplot(data.frame()) + geom_point() + xlim(0, 10) + ylim(0, 100) + theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
	#p.first <- plot_grid(p.blank,p.arc,nrow=1, rel_widths=c(1/7,1))
	#p.out   <- plot_grid(p.heat,p.first,ncol=1, axis="l", rel_heights=c(1,1/4)) #scale=c(1,0.9)
	outstr <- paste(opt$output, "_heatmap", a, ".pdf", sep="")
	a <- a + 1
	pdf(outstr,height=58,width=55)
	#print(p.out)
	print(p.heat)
	dev.off()
}
#p.heat <- gheatmap(p.tree, annot, offset=1, width=8, font.size=2, colnames_angle=-90, hjust=0) + #offset=0.009
#         guides(fill=FALSE) +
#         scale_fill_manual(breaks=heatmap.breaks, values=node.colour) +
#         theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#         theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
#outstr <- paste(opt$output, "_heatmap.pdf", sep="")
#pdf(outstr,height=58,width=55)
#print(p.heat)
#dev.off()

#Draw network
#grp <- graph_from_data_frame(d = edges, vertices = node.order, directed = FALSE)
#net.layout <- create_layout(grp, layout='igraph', algorithm='fr') #dh, graphopt
#node.colour.order <- sort(node.colour)
#p.net <- ggraph(net.layout) +
#  geom_edge_link(aes(colour=(1-E(grp)$p)), width=2, alpha = 0.8) + #width=(1-E(grp)$p), 
#  scale_edge_color_gradient2(low = "#f0f0f0",
#                            mid = "gray",
#                            high = "gray",
#                            midpoint = 0.005,
#                            #trans = "log",
#                            breaks=arc.breaks,
#                            labels=arc.breaks
#			    ) +
#  geom_node_point(aes(color=node.colour), size=10) +
#  scale_color_manual(values=unique(node.colour.order)) +
#  geom_node_text(aes(label = node.order), vjust=-1.2) +
#  theme_graph() +
#  theme(legend.position="none") + #bottom
#  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
#outstr <- paste(opt$output, "_network.pdf", sep="")
#pdf(outstr,height=58,width=55)
#print(p.net)
#dev.off()
