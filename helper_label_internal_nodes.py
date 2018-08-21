from ete3 import Tree
import argparse

parser = argparse.ArgumentParser(description='File name containing tree.')
parser.add_argument('treein', type=str, help='File name containing tree.')
args = parser.parse_args()

tree = Tree(args.treein)

edge = 0
for node in tree.traverse():
   if not node.is_leaf():
      node.name = "NODE_%d" %edge
      edge += 1

print(tree.write(format=3))
