import csv
import subprocess
import argparse
import os

#Motive: take an input file of the format geneFamily \t genome and output a roary style gene_presence_absence.csv.

#Parse arguemnts to get input file for reading
parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i", type=str, required=True)
args = parser.parse_args()

#Ensure roary-style output file doesn't already exist
exists = os.path.isfile('gene_presence_absence.csv')
if (exists):
    print("gene-presence-absence.csv already exists; I don't want to overwrite it.")
    print("Exiting...")
    quit()
#Open roary-style output file
roary = open('gene-presence-absence.csv', "w")
#Sort input file by gene ID
ret = subprocess.call("sort "+args.input+" > sorted.tmp", shell=True)
#Collect and count unique genomes
ret = subprocess.Popen('cut -f 2 sorted.tmp | sort | uniq', stdout=subprocess.PIPE, shell=True)
(uniq_gens, err) = ret.communicate()
uniq_gens = uniq_gens.decode('ascii').strip()
gen_list = uniq_gens.split("\n")
loc_len = len(gen_list)
#Make and populate genome->location hash
hash_count = 0
header = ""
gen_hash = {}
for gen in gen_list:
    gen_hash[gen] = hash_count
    hash_count+=1
    header = header+",\""+gen+"\""
#Write header line to roary
roary.write("\"Gene\",\"Non-unique Gene name\",\"Annotation\",\"No. isolates\",\"No. sequences\",\"Avg sequences per isolate\",\"Genome Fragment\",\"Order within Fragment\",\"Accessory Fragment\",\"Accessory Order with Fragment\",\"QC\",\"Min group size nuc\",\"Max group size nuc\",\"Avg group size nuc\""+header+"\n")
#While input file isn't empty
with open(args.input) as f:
    line = f.readline()
    while True:
        #Make/empty genomeloc array
        genomeloc = ["" for x in range(loc_len)]
        if not line:
            break #EOF
        geneID = line.split("\t")[0]
        while (line.split("\t")[0] == geneID):
            #Get location for corresponding genome
            gen = (line.split("\t")[1]).strip()
            x = gen_hash[gen]
            #Append current entry to array
            if (genomeloc[x] == ""):
                genomeloc[x] = gen+"_"+geneID
            else:
                geneomeloc[x] = genomeloc[x]+" "+geneID+"_"+gen
            line = f.readline()
        bulk = ("".join([',"{}"'.format(genomeloc[n]) for n in range(len(genomeloc))]))
        roary.write("\""+geneID+"\",\"\",\"\",\"\",\"\",\"\",\"\",\"\",\"\",\"\",\"\",\"\",\"\",\"\""+bulk+"\n")
roary.close()
