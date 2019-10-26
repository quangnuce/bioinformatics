import os, shutil, glob
import csv
import pandas as pd
import json
import itertools
import gffutils
import pyfaidx
from Bio.Seq import Seq
from Bio import SeqIO
from sys import argv,stderr
from multiprocessing import Pool

#python combine2MLinput_kmer.py [Fasta folder] [Output folder] [k] [Prokka output folder] [RGI output folder]
FASTA_IN=argv[1]
output=argv[2]
k=int(argv[3])
PROKKA_OUT=argv[4]
RGI_OUT=argv[5]
KMC_OUTPUT=output+'/KMC_OUTPUT'
KMER_OUTPUT=output+'/KMER_OUTPUT'
def countingKmer(fasta_folder,kmc_output,k):
    #generate kmer counting files from fasta files using kmc_dump
    print('counting kmer in ',kmc_output,' folder')
    if not os.path.exists(kmc_output):
        os.makedirs(kmc_output)
    for root, dirs, files in os.walk(fasta_folder):
        for _file in files:
            if _file.endswith(('.fna','.fasta','.fn')):
                print ('Found file in: ' , str(root),_file)
                newname = os.path.splitext(os.path.basename(_file))[0]
                print ('Basename ' , newname);
                cmd='kmc -k'+str(k)+' -ci1 -fa -fm '+fasta_folder+'/'+_file+ ' '+kmc_output+'/'+newname+' '+ kmc_output+'/'+newname+'_tmp'
                os.system(cmd)
                cmd_dump='kmc_dump '+kmc_output+'/'+newname+' '+kmc_output+'/'+newname+'.'+str(k)+'.kmrs'
                os.system(cmd_dump)
                if os.path.exists(kmc_output+'/'+newname+'.kmc_suf'):
                    os.remove(kmc_output+'/'+newname+'.kmc_suf')
                if os.path.exists(kmc_output+'/'+newname+'.kmc_pre'):
                    os.remove(kmc_output+'/'+newname+'.kmc_pre')

def generateKmerInput(fasta_folder,kmc_output,ml_input,k):
    countingKmer(fasta_folder,kmc_output,k)
    #generate all kmer
    print("Generate whole genome ",k,"-mer files")
    if not os.path.exists(ml_input+"/kmer"):
        os.makedirs(ml_input+"/whole_genome_kmers")

    # init all posible k-mer
    kmerHash = {}
    bases=['A','C','G','T']
    #k=10
    kmers=[''.join(p) for p in itertools.product(bases, repeat=k)]
    #hash_kmer={}
    for i in range(len(kmers)):
        kmerHash[kmers[i]]=0;
    #mark the kmer if it exist in genome collection
    for root, dirs, files in os.walk(kmc_output):
        for _file in files:
            if _file.endswith('.kmrs'):
                f = open(str(root)+'/'+_file)
                for i in f:
                    i = i.strip().split('\t')
                    kmerHash[i[0]] = 1
                f.close()
    #make kmer bag hold all kmer exist in genome collection, save the index file
    file1 = open(ml_input+"/whole_genome_kmers"+"/WholeGenomeIndex.txt","w")
    kmerBag={}
    j=0
    for i in kmerHash:
        if kmerHash[i]==1:
            j+=1
            kmerBag[i]=j
            file1.write(str(i)+'\t'+str(j)+'\n')
    file1.close()
    # generate csv file for each strain
    for root, dirs, files in os.walk(kmc_output):
        for _file in files:
            if _file.endswith('.kmrs'):
                generateSample(str(root)+'/'+_file,kmerBag,ml_input+'/whole_genome_kmers/'+str(k))
def generateKmerAMRInput(prokka_out,rgi_out,output,k):
    #generate all kmer
    print("Generate amr k-mer files")
    if not os.path.exists(output+"/amr_genome_kmers"):
        os.makedirs(output+"/amr_genome_kmers")

    print ('k',k)
    # open up the amr myfiles, make list of amr genes id for each strain
    print("getting amr gene for each train...")
    list_amr_by_strains={}
    for root, dirs, files in os.walk(rgi_out):
        for _file in files:
            if _file.endswith('.json'):
                with open (rgi_out+'/'+_file, 'r') as f:
                    jsonreader =json.load(f)
                    basename=os.path.basename(_file)
                    strainname = basename.rsplit('.', 1)[0]
                    amr_genes=[]
                    for item in jsonreader.values():
                        if isinstance(item,dict):
                            for i, v in item.items():
                                if isinstance(v,dict):
                                    amr_genes.append(v['query_from'])
                    list_amr_by_strains[strainname]=amr_genes
    #print(list_amr_by_strains)
    print("getting done, found ",len(list_amr_by_strains))
    #read fnn file to collect dna sequence by amr genes, save in new fasta file
    print("creating fasta file contain amr genes for each train")
    if not os.path.exists(output+"/amr_genome_kmers/fasta"):
        os.makedirs(output+"/amr_genome_kmers/fasta")
    list_fna_files_by_strain={}
    for root, dirs, files in os.walk(prokka_out):
        for _file in files:
            if _file.endswith('.ffn'):
                #strainname = os.path.basename(str(root)).rsplit('.', 1)[0]#remember to delete split when re-running pipeline
                strainname = os.path.basename(str(root))
                ofile = open(output+"/amr_genome_kmers/fasta/"+strainname+".fn", "w")
                for record in SeqIO.parse(str(root)+'/'+_file, "fasta"):
                    seq=record.seq
                    #print(record.description)
                    if record.description in list_amr_by_strains[strainname]:
                        #list_amr_gene_id.append(record.description)
                        #list_amr_gene_seq.append(record.seq)
                        SeqIO.write(record,ofile,"fasta")
                ofile.close()
    print("create done, save in ",output,"/amr_genome_kmers/fasta/")
    #counting kmer with amr fasta
    print("counting kmer")

    countingKmer(output+"/amr_genome_kmers/fasta",output+"/amr_genome_kmers/kmerdump",k)
    print("counting done, create bag of kmer for amr genes")
    kmerHash = {}
    bases=['A','C','G','T']
    print('k=',k)
    kmers=[''.join(p) for p in itertools.product(bases, repeat=int(k))]
    for i in range(len(kmers)):
        kmerHash[kmers[i]]=0;
    #mark the kmer if it exist in genome collection
    for root, dirs, files in os.walk(output+"/amr_genome_kmers/kmerdump"):
        for _file in files:
            if _file.endswith('.kmrs'):
                f = open(str(root)+'/'+_file)
                for i in f:
                    i = i.strip().split('\t')
                    kmerHash[i[0]] = 1
                f.close()
    file1 = open(output+"/amr_genome_kmers/AMRGenomeIndex.txt","w")
    kmerBag={}
    j=0
    for i in kmerHash:
        if kmerHash[i]==1:
            j+=1
            kmerBag[i]=j
            file1.write(str(i)+'\t'+str(j)+'\n')
    file1.close()
    #print('bag of kmer done, save in ',output,"/amr_genome_kmers/AMRGenomeIndex.txt")
    # generate csv file for each strain
    print("generate ml input for each strain")
    for root, dirs, files in os.walk(output+"/amr_genome_kmers/kmerdump"):
        for _file in files:
            if _file.endswith('.kmrs'):
                generateSample(str(root)+'/'+_file,kmerBag,output+'/amr_genome_kmers/'+str(k))
def generateKmerNonAMRInput(prokka_out,rgi_out,output,k):
    #generate all kmer
    print("Generate amr k-mer files")
    if not os.path.exists(output+"/non_amr_genome_kmers"):
        os.makedirs(output+"/non_amr_genome_kmers")

    print ('k',k)
    # open up the amr myfiles, make list of amr genes id for each strain
    print("getting amr gene for each train...")
    list_amr_by_strains={}
    for root, dirs, files in os.walk(rgi_out):
        for _file in files:
            if _file.endswith('.json'):
                with open (rgi_out+'/'+_file, 'r') as f:
                    jsonreader =json.load(f)
                    basename=os.path.basename(_file)
                    strainname = basename.rsplit('.', 1)[0]
                    amr_genes=[]
                    for item in jsonreader.values():
                        if isinstance(item,dict):
                            for i, v in item.items():
                                if isinstance(v,dict):
                                    amr_genes.append(v['query_from'])
                    list_amr_by_strains[strainname]=amr_genes
    #print(list_amr_by_strains)
    print("getting done, found ",len(list_amr_by_strains))
    #read fnn file to collect dna sequence by amr genes, save in new fasta file
    print("creating fasta file contain amr genes for each train")
    if not os.path.exists(output+"/non_amr_genome_kmers/fasta"):
        os.makedirs(output+"/non_amr_genome_kmers/fasta")
    list_fna_files_by_strain={}
    for root, dirs, files in os.walk(prokka_out):
        for _file in files:
            if _file.endswith('.ffn'):
                #strainname = os.path.basename(str(root)).rsplit('.', 1)[0]#remember to delete split when re-running pipeline
                strainname = os.path.basename(str(root))
                ofile = open(output+"/non_amr_genome_kmers/fasta/"+strainname+".fn", "w")
                for record in SeqIO.parse(str(root)+'/'+_file, "fasta"):
                    seq=record.seq
                    #print(record.description)
                    if record.description in list_amr_by_strains[strainname]:
                        print('ignore ',record.description)
                    else:
                        SeqIO.write(record,ofile,"fasta")
                ofile.close()
    print("create done, save in ",output,"/non_amr_genome_kmers/fasta/")
    #counting kmer with amr fasta
    print("counting kmer")

    countingKmer(output+"/non_amr_genome_kmers/fasta",output+"/non_amr_genome_kmers/kmerdump",k)
    print("counting done, create bag of kmer for non amr genes")
    kmerHash = {}
    bases=['A','C','G','T']
    print('k=',k)
    kmers=[''.join(p) for p in itertools.product(bases, repeat=int(k))]
    for i in range(len(kmers)):
        kmerHash[kmers[i]]=0;
    #mark the kmer if it exist in genome collection
    for root, dirs, files in os.walk(output+"/non_amr_genome_kmers/kmerdump"):
        for _file in files:
            if _file.endswith('.kmrs'):
                f = open(str(root)+'/'+_file)
                for i in f:
                    i = i.strip().split('\t')
                    kmerHash[i[0]] = 1
                f.close()
    file1 = open(output+"/non_amr_genome_kmers/NonAMRGenomeIndex.txt","w")
    kmerBag={}
    j=0
    for i in kmerHash:
        if kmerHash[i]==1:
            j+=1
            kmerBag[i]=j
            file1.write(str(i)+'\t'+str(j)+'\n')
    file1.close()
    #print('bag of kmer done, save in ',output,"/amr_genome_kmers/AMRGenomeIndex.txt")
    # generate csv file for each strain
    print("generate ml input for each strain")
    for root, dirs, files in os.walk(output+"/non_amr_genome_kmers/kmerdump"):
        for _file in files:
            if _file.endswith('.kmrs'):
                generateSample(str(root)+'/'+_file,kmerBag,output+'/non_amr_genome_kmers/'+str(k))
def generateSample(kmcfile,featureHash,ml_input):
    #print("Call function generateSample")
    if not os.path.exists(ml_input):
        os.makedirs(ml_input)
    f = open(kmcfile)
    print(kmcfile)
    # init kmc contigs hash
    kmcContigs = {}
    # for each kmr in file, map kmr to count in hash
    for i in f:
        i = i.strip().split('\t')
        kmcContigs[i[0]] = i[1]
    f.close()
    kmer=[]
    count=[]
    for i in featureHash:
        kmer.append(i)
        if i in kmcContigs:
            count.append(kmcContigs[i])
        else:
            count.append(0)
    frameKmer= pd.DataFrame({"kmer":kmer})
    frameCount=pd.DataFrame({"count":count})
    frameKmer=frameKmer.join(frameCount)
    frameKmer.to_csv(ml_input+'/'+os.path.splitext(os.path.basename(kmcfile))[0]+'.kmer.csv')
def main():
    generateKmerInput(FASTA_IN,KMC_OUTPUT,output,k)
    generateKmerAMRInput(PROKKA_OUT,RGI_OUT,output,k)
    generateKmerNonAMRInput(PROKKA_OUT,RGI_OUT,output,k)
if __name__== "__main__" :
    main()
