import ete3
import itertools

#--Addition to coinfinder--#
#Input:
# (a) edgeList_B1 and edgeList_B2 (2 elements for which we are testing for co-occurrence)
#       A1--B1
#       A3--B1
#       A1--B2
#       A2--B2 etc.
# (b) newick formatted phylogeny of the As
#
#Output: an array of pairwise distances between all combinations of B1 and B2 edges (maybe a 2d array with labelled A pairs?)
def calc( *thelist ):
	phylo = ete3.Tree(str(thelist[0]))
	edges = []
	dists = list()
	maxi = 0
	#for nodeA in phylo.traverse("postorder"):
	#	edges.append(nodeA.name)
	edges = phylo.get_leaf_names()
	for perm in itertools.combinations(edges, 2):
		nodeA = phylo&(str(perm[0]))
		nodeB = phylo&(str(perm[1]))
		curdist = nodeA.get_distance(nodeB)
		dists.append(nodeA.name)
		dists.append(nodeB.name)
		dists.append(curdist)
		#dists[nodeA.name,nodeB.name] = float(416.5)
		#dists['Toronto','Hamilton'] = 416
		#if (curdist > maxi):
		#        maxi = curdist
	return(dists)
