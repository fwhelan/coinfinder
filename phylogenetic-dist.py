import itertools
import ete3
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

#For now, dummy data:
edgeList_B1 = ['A1', 'A3']
edgeList_B2 = ['A1', 'A2']

#Get local copies of edgeLists from coinfinder proper
#---todo----

#Read in phylogenetic tree from file
treeFile = "test.phylo.newick" #get phylogenetic tree file name from coinfinder proper----todo----
phylo = ete3.Tree(treeFile)

#Define distanceList
distList = list()

#For each unique combination between edgeList_B1 and edgeList_B2, calculate the pairwise phylogenetic distance
for perm in itertools.product(edgeList_B1, edgeList_B2):
	#print(perm[0] + perm[1])
	nodeA = phylo&perm[0]
	nodeB = phylo&perm[1]
	#print("The distance calcs:", nodeA, nodeB, nodeA.get_distance(nodeB))
	distList.append(nodeA.get_distance(nodeB))
print(distList)
