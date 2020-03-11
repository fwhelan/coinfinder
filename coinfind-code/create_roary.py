import csv
import subprocess
import argparse
import os
import sys

#Motive: take an input file of the format geneFamily \t genome and output a roary style gene_presence_absence.csv.

#Parse arguemnts to get input file for reading
parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i", type=str, required=True)
args = parser.parse_args()

#Ensure roary-style output file doesn't already exist
#exists = os.path.isfile('gene_presence_absence.csv')
#if (exists):
#    print("gene_presence_absence.csv already exists; I don't want to overwrite it.")
#    print("Exiting...")
#    quit()
#Ensure input file is of the format geneFamily \t genome \n
infile = open(args.input, 'r')
line = infile.readline()
linearr = line.split("\t")
if len(linearr) != 2:
    print("Error: input file "+args.input+" does not appear to be in the format geneFamily<tab>genome")
    print("Did you mean to include the -I flag to indicate a Roary formatted input file?")
    print("Exiting...")
    quit()
#Open roary-style output file
roary = open('gene_presence_absence.csv', "w")
#Sort input file by gene ID
command="sort -t '\t' -k1,1 "
ret = subprocess.call("sort -t '\t' -k1,1 "+args.input+" > sorted.tmp", shell=True)
#Collect and count unique genomes
ret = subprocess.Popen('cut -f 2 sorted.tmp | sort | uniq', stdout=subprocess.PIPE, shell=True)
(uniq_gens, err) = ret.communicate()
uniq_gens = uniq_gens.decode('ascii').strip()
gen_list = uniq_gens.split("\n")
gen_list = [ i.strip() for i in gen_list ]
loc_len = len(gen_list)
#Make and populate genome->location hash
hash_count = 0
header = ""
#header = "\""
gen_hash = {}
for gen in gen_list:
    gen_hash[gen] = hash_count
    hash_count+=1
    header = header+",\""+gen+"\""
#Write header line to roary
roary.write("\"Gene\",\"Non-unique Gene name\",\"Annotation\",\"No. isolates\",\"No. sequences\",\"Avg sequences per isolate\",\"Genome Fragment\",\"Order within Fragment\",\"Accessory Fragment\",\"Accessory Order with Fragment\",\"QC\",\"Min group size nuc\",\"Max group size nuc\",\"Avg group size nuc\""+header+"\n")
#While input file isn't empty
#with open(args.input) as f:
with open("sorted.tmp") as f:
    line = f.readline()
    while True:
        if not line:
            break #EOF
        #Make/empty genomeloc array
        genomeloc = ["" for x in range(loc_len)]
        try:
            geneID = line.split("\t")[0]
        except:
            print("Cannot create geneID from line: " + line)
            exit()
        while (line.split("\t")[0] == geneID):
            #Get location for corresponding genome
            gen = (line.split("\t")[1]).strip()
            x = gen_hash[gen]
            #Append current entry to array
            if (genomeloc[x] == ""):
                genomeloc[x] = gen+"_"+geneID
            else:
                genomeloc[x] = genomeloc[x]+" "+gen+"_"+geneID
            line = f.readline()
        bulk = ("".join([',"{}"'.format(genomeloc[n]) for n in range(len(genomeloc))]))
        roary.write("\""+geneID+"\",\"\",\"\",\"\",\"\",\"\",\"\",\"\",\"\",\"\",\"\",\"\",\"\",\"\""+bulk+"\n")
roary.close()
