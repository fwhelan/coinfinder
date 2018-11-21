import csv
import re
import subprocess
import argparse

#Parse arguments to get input file for reading
parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i", type=str, required=True)
args = parser.parse_args()

#Open edge output files
edge = open('concident-input-edges.csv', "w")
#Coinfinder won't like to see a header in the edges file
#Open roary output/coinfinder input file for reading
with open(args.input) as genefile:
    genereader = csv.reader(genefile)
    #Skip header
    next(genereader, None)
    for line in genereader:
        #Ensure no tabs in geneID
        gene = re.sub(r'\t',r'',line[0])
        #Grab genome position, but only if it is set in the input
        genePos = 0
        if((isinstance(line[6],int)) and (isinstance(line[7],int))):
            genePos = (int(line[6])-1)*100000 + int(line[7])
        #Grab all gene identifiers per genome
        for i in range(14,len(line)):
            #If not empty..
            if (line[i] != ""):
                #Catch any multi-gene entries per genome
                splity = line[i].split()
                for j in range(0,len(splity)):
                    try:
                        genome = re.search('(.+)_.*', splity[j]).group(1)
                        geneID = re.search('.+_(.*)', splity[j]).group(1)
                    except AttributeError:
                        genome = ''
                        print("AttributeError caught, genome not defined.")
                        print("line[i], {}".format(line[i]))
                        print("{}".format(line))
                        break
                    #Ensure no tabs in genomeID
                    genomeID = re.sub(r'\t',r'',genome)
                    edge.write("{}\t{}\t{}\n".format(gene,genomeID,genePos))
edge.close()
