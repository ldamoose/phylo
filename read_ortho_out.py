# input: folder with Orthogroups.csv; outfolder path; folder with proteins
# output: folder with multifasta of SINGLE COPY ORTHOGROUPS, presence-absence matrix 
#########################
def orthos2fasta(outfolder, suffix, ort_dict, prot_dict, suff):    
    makedir(path.abspath(outfolder) + "/" + suffix)
    for og in ort_dict.keys():
        outfile = open(path.abspath(outfolder) + "/" + suffix + "/" + og + suff, "w")
        group = ort_dict[og]
        for org in group.taxa:
            c = 1
            if suffix == "single":
                outfile.write(">" + org + "\n")
                outfile.write(str(prot_dict[org][group.taxa_dict[org][0].split(" ")[0]].seq) + "\n")
            else:
                for seq in group.taxa_dict[org]:
                    outfile.write(">" + org + "_" + str(c) + "\n")
                    outfile.write(str(prot_dict[org][seq.strip()].seq) + "\n")
                    c += 1
        outfile.close()

class orthogroup:
    # initializes istance,dictionary of taxa:orthomembers, and status (single copy or not)
    def __init__(self,taxa):
        self.taxa = taxa
        self.taxa_dict = {t:[] for t in self.taxa}
        self.single = "empty"
        self.common = "empty"
        self.taxa_list = []
    # fills dictionary of taxa:orthomembers
    def fill_taxa_dict(self,og_members):
        c = 0                                   #initialize single copy counter
        a = 0                                   #initialize common orthogroup counter
        for i in range(len(og_members)):
            t, m = self.taxa[i], og_members[i]  #organism and respective members of the og
            if m != m:                          #if m is 'nan'
                m = "NA"
            else:
                self.taxa_list += [t]
                a += 1                          #counts how many taxa have the orthogroup
                m = m.split(",")                #split the members -comma separated-
                if len(m) == 1:                 #if single member -> keep count
                    c += 1
            self.taxa_dict[t] = m         
    # this variable stores single copy orthogroup status (T/F)
        if c == len(self.taxa):                 #if all taxa have a single copy -> this is a single copy orthogroup
            self.single = True
        else:
            self.single = False
        if a == len(self.taxa):
            self.common = True
        else:
            self.common = False

def makedir(path):
    import os
    path = os.path.abspath(path)
    if not os.path.exists(path):
        os.makedirs(path)
    
        
def main():
    # read file
    infile = open(args.infile, "r")      
    orthos = pd.DataFrame.from_csv(args.infile, sep='\t', header=0)
    # define taxa
    taxa = list(orthos.columns)
    if args.molecule.startswith("n"):
        suff = ".fna"
        o = "prot"
        i = "nuc"
    elif args.molecule.startswith("p"):
        suff = ".faa"
        o = ""
        i = ""
    taxa = [x.replace(o,i) for x in taxa]
    #og_id, og_members
    c = 0
    sco_dict = {}   # Single Copy Orthogroups dictionary
    com_dict = {}   # Common Orthogroups dictionary
    # nucleotides or proteins:
    # open output file for P/A matrix
    makedir(args.outfolder)
    outfile = open(path.abspath(args.outfolder) + "/OrthoFinder_PAmatrix.tsv", "w")
    outfile.write("\t".join(["OrthoGroupID"] + taxa) + "\n")
    for row in orthos.itertuples():
        og = orthogroup(taxa)
        row = list(row)
        og_id, og_members = row[0], row[1:]
        # fill the internal dictionary and get info about single copy group or not
        og.fill_taxa_dict(og_members)
        # write out line of P/A matrix
        outfile.write("\t".join([og_id] + [str(len(og.taxa_dict[t])) if og.taxa_dict[t] is not "NA" else "0" for t in taxa]) + "\n")
        # select single copy orthogroups
        if og.single:
            sco_dict[og_id] = og
        # select common orthogroups
        if og.common:
            #if len( og.taxa_list) != len(taxa):
            com_dict[og_id] = og
    outfile.close()
    # create directory and multifasta of Single Copy Orthogroups
    # first load all proteins on memory
    proteomes = [path.abspath(args.proteins) + "/" + t + suff for t in taxa]
    try:
        prot_dict = {path.basename(org.replace(suff,"")).replace(o,i) : SeqIO.to_dict(SeqIO.parse(org.replace(o,i), "fasta")) for org in proteomes}
    except OSError:
        print "check the path of your protein files!"
    # then for each single copy orthogroup, create a multifasta file in a subfolder of the output folder
    orthos2fasta(args.outfolder + "/", "single", sco_dict, prot_dict, suff)
    orthos2fasta(args.outfolder + "/", "common", com_dict, prot_dict, suff)
        
    


if __name__ == "__main__":
    from os import path, system, makedirs
    import argparse
    import pandas as pd
    from Bio import SeqIO
    import numpy as np
    import glob
    parser = argparse.ArgumentParser()
    parser.add_argument("--infile","-i", help = "Orthogroups.csv file path")
    parser.add_argument("--proteins","-p", help = "Path of protein files folder")
    parser.add_argument("--outfolder", "-o", help = "output folder path")    
    parser.add_argument("--molecule", "-m", default="p",help = "p(rotein) or n(ucleotide)")
    args = parser.parse_args()
    main()
    
