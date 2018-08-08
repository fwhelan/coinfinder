#import ete3
#import itertools

def calc( *thelist ):
	phylo = ete3.Tree(str(thelist[0]), format=8)
	edges_union = thelist[1:]
	#Find the common ancestor internal node in phylo of the edges in edges_union and calculate its distance to the root node
	ancestor = phylo.get_common_ancestor(edges_union)
	dist = ancestor.get_distance(phylo.get_tree_root())
	name = ancestor.name
	if (ancestor.is_root()):
		return("Root")
	return(str(name))
