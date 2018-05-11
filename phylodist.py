import ete3
import itertools

#--Addition to coinfinder--#
#Input:
# (a) edgeList_B1 and edgeList_B2 (2 elements for which we are testing for co-occurrence)
#	A1--B1
#	A3--B1
#	A1--B2
#	A2--B2 etc.
# (b) newick formatted phylogeny of the As
#
#Output: an array of pairwise distances between all combinations of B1 and B2 edges (maybe a 2d array with labelled A pairs?)

def calc( *thelist ):
	phylo = ete3.Tree(str(thelist[0]))
	count = thelist[1]
	edgesB1 = thelist[2:2+count]
	edgesB2 = thelist[2+count:]
	#Define distanceList
	distList = list()
	#For each unique combination between edgesB1 and edgesB2, calculate the pairwise phylogenetic distance
	for perm in itertools.combinations(edgesB1, edgesB2):
		try:
			nodeA = phylo&(str(perm[0]))
		except:
			print("Cannot find node {} in phylogeny.".format(perm[0]))
			return(perm[0])
		try:
			nodeB = phylo&(str(perm[1]))
		except:
			print("Cannot find node {} in phylogeny.".format(perm[1]))
			return(perm[1])
		distList.append(nodeA.get_distance(nodeB))
		#add catch for nodes not in the tree
	print(str(distList))
	return(distList)
