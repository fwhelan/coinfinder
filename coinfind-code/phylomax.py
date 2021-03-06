import ete3
import itertools

#Addition to coinfinder#
#Input:
# (a) edgeList_B1 and edgeList_B2 (2 elements for which we are testing for cooccurrence)
#       A1B1
#       A3B1
#       A1B2
#       A2B2 etc.
# (b) newick formatted phylogeny of the As
#
#Output: an array of pairwise distances between all combinations of B1 and B2 edges (maybe a 2d array with labelled A pairs?)
def calc( *thelist ):
	phylo = ete3.Tree(str(thelist[0]), format=1)
	edges = {}
	dists = list()
	#maxi = 0
	#for nodeA in phylo.traverse("postorder"):
	#	edges.append(nodeA.name)
	edges = {leaf.name : {} for leaf in phylo}
	for perm in itertools.combinations(list(edges.keys()), 2):
		curdist = phylo.get_distance(perm[0],perm[1])
		edges[perm[0]][perm[1]]=curdist
		#nodeA = phylo&(str(perm[0]))
		#nodeB = phylo&(str(perm[1]))
		#curdist = nodeA.get_distance(nodeB)
		#dists.append(nodeA.name)
		#dists.append(nodeB.name)
		#dists.append(curdist)
		#dists[nodeA.name,nodeB.name] = float(416.5)
		#dists['Toronto','Hamilton'] = 416
		#if (curdist > maxi):
		#        maxi = curdist
	else:
		for nodeA in edges:
			for nodeB in edges[nodeA]:
				dists.extend((nodeA,nodeB,edges[nodeA][nodeB]))
	return(dists)
