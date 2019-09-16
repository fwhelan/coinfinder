import csv
import re

#Open test coinfinder dataset with left and middle defined; add right to .out of it
left.mid = open('genome-gene.csv', "r")
left.mid.right = open('genome-gene.csv.out', "w")

#Open roary output file for reading
with open('gene_presence_absence.csv') as infile:
	filereader = csv.reader(infile)
	#Skip header line
	next(filereader, None)
	for line in filereader:
		
