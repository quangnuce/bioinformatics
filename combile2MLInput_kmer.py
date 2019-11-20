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
#generate kmer counting files from fasta files using kmc_dump
def countingKmer(fasta_folder,kmc_output,k):
    '''
        Generate kmer counting files from fasta files using kmc_dump
        :fasta_folder: the folder contains fasta files
        :kmc_output: the folder stored dump file
        :k: k in k-mer
        :return: dump files
    '''
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
                cmd_dump='kmc_dump '+kmc_output+'/'+newname+' '+kmc_output+'/'+newname+'.kmrs'
                os.system(cmd_dump)
                if os.path.exists(kmc_output+'/'+newname+'.kmc_suf'):
                    os.remove(kmc_output+'/'+newname+'.kmc_suf')
                if os.path.exists(kmc_output+'/'+newname+'.kmc_pre'):
                    os.remove(kmc_output+'/'+newname+'.kmc_pre')

def generateKmerAllInput(fasta_folder,output,k):
    '''
        Generate kmer counting file for all kmer exist in whole genome
        :fasta_folder: the folder contains fasta files
        :kmc_output: the folder stored dump file
        :k: k in k-mer
        :output: root folder to store exported files (kmer counting file and kmer index file)
    '''
    output_for_all_kmer=output+"/whole_genome_kmers"
    if not os.path.exists(output_for_all_kmer):
        os.makedirs(output_for_all_kmer)
    countingKmer(fasta_folder,output_for_all_kmer+"/kmerdump",k)
    #generate all kmer
    print("Generate whole genome ",k,"-mer files")

    #get all kmer exist in
    kmerBag=makeKmerBag(k,output_for_all_kmer+"/kmerdump",output_for_all_kmer+"/WholeGenomeIndex.txt")
    # generate csv file for all strain
    allkmer_df=pd.DataFrame({"kmer":list(kmerBag.keys())});
    print(allkmer_df)
    for root, dirs, files in os.walk(output_for_all_kmer+"/kmerdump"):
        for _file in files:
            if _file.endswith('.kmrs'):
                allkmer_df=allkmer_df.merge(countKmerOnFile(str(root)+'/'+_file,kmerBag), how='left', on='kmer')
    allkmer_df.to_csv(output_for_all_kmer+'/whole_genome'+'.kmer.csv')
def generateKmerAMRInput(prokka_out,rgi_out,output,k):
    '''
        Generate kmer counting file for kmers exist in AMR genes only. List AMR genes (annoted by RGI) for each strain is used to
        build a new fasta file.

        :prokka_out: the output folder of prokka
        :rgi_out: the output folder of RGI
        :k: k in k-mer
        :output: root folder to store exported files (kmer counting file and kmer index file)
    '''
    print("Generate amr k-mer files")
    output_for_amr_kmer=output+"/amr_genome_kmers"
    if not os.path.exists(output_for_amr_kmer):
        os.makedirs(output_for_amr_kmer)
    # open up the amr myfiles, make list of amr genes id for each strain
    print("getting amr gene for each train...")
    list_amr_by_strains=makeDictForAMRGene(rgi_out)
    #print(list_amr_by_strains)
    print("getting done, found ",len(list_amr_by_strains))
    #read fnn file to collect dna sequence by amr genes, save in new fasta file
    print("creating fasta file contain amr genes for each train")
    if not os.path.exists(output+"/amr_genome_kmers/fasta"):
        os.makedirs(output_for_amr_kmer+"/fasta")
    list_fna_files_by_strain={}
    for root, dirs, files in os.walk(prokka_out):
        for _file in files:
            if _file.endswith('.ffn'):
                strainname = os.path.basename(str(root)).rsplit('.', 1)[0]#remember to delete split when re-running pipeline
                #strainname = os.path.basename(str(root))
                ofile = open(output_for_amr_kmer+"/fasta/"+strainname+".fn", "w")
                for record in SeqIO.parse(str(root)+'/'+_file, "fasta"):
                    seq=record.seq
                    #print(record.description)
                    if record.description in list_amr_by_strains[strainname]:
                        #list_amr_gene_id.append(record.description)
                        #list_amr_gene_seq.append(record.seq)
                        SeqIO.write(record,ofile,"fasta")
                ofile.close()
    print("create done, save in ",output_for_amr_kmer+"/fasta/")
    #counting kmer with amr fasta
    print("counting kmer")

    countingKmer(output_for_amr_kmer+"/fasta",output_for_amr_kmer+"/kmerdump",k)
    print("counting done, create bag of kmer for amr genes")
    kmerBag=makeKmerBag(k,output_for_amr_kmer+"/kmerdump",output_for_amr_kmer+'/amr_kmer_index.txt')

    #print('bag of kmer done, save in ',output,"/amr_genome_kmers/AMRGenomeIndex.txt")
    # generate csv file for each strain
    # generate csv file for all strain
    amrkmer_df=pd.DataFrame({"kmer":list(kmerBag.keys())});
    print("generate ml input for each strain")
    for root, dirs, files in os.walk(output_for_amr_kmer+"/kmerdump"):
        for _file in files:
            if _file.endswith('.kmrs'):
                amrkmer_df=amrkmer_df.merge(countKmerOnFile(str(root)+'/'+_file,kmerBag), how='left', on='kmer')
    amrkmer_df.to_csv(output_for_amr_kmer+'/amr_genome'+'.kmer.csv')
def generateKmerNonAMRInput(prokka_out,rgi_out,output,k):
    '''
        Generate kmer counting file for kmers exist in non AMR genes only.
        :prokka_out: the output folder of prokka
        :rgi_out: the output folder of RGI
        :k: k in k-mer
        :output: root folder to store exported files (kmer counting file and kmer index file)
    '''
    #generate non amr kmer
    print("Generate non amr k-mer files")
    output_for_non_amr_kmer= output+"/non_amr_genome_kmers"
    if not os.path.exists(output_for_non_amr_kmer):
        os.makedirs(output_for_non_amr_kmer)

    print ('k',k)
    # open up the amr myfiles, make list of amr genes id for each strain
    print("getting amr gene for each train...")
    list_amr_by_strains=makeDictForAMRGene(rgi_out)
    #print(list_amr_by_strains)
    print("getting done, found ",len(list_amr_by_strains))
    #read fnn file to collect dna sequence by non amr genes (except amr gene), save in new fasta file
    print("creating fasta file contain amr genes for each train")
    if not os.path.exists(output_for_non_amr_kmer+"/fasta"):
        os.makedirs(output_for_non_amr_kmer+"/fasta")
    list_fna_files_by_strain={}
    for root, dirs, files in os.walk(prokka_out):
        for _file in files:
            if _file.endswith('.ffn'):
                strainname = os.path.basename(str(root)).rsplit('.', 1)[0]#remember to delete split when re-running pipeline
                #strainname = os.path.basename(str(root))
                ofile = open(output_for_non_amr_kmer+"/fasta/"+strainname+".fn", "w")
                for record in SeqIO.parse(str(root)+'/'+_file, "fasta"):
                    seq=record.seq
                    #print(record.description)
                    if record.description in list_amr_by_strains[strainname]:
                        print('ignore ',record.description)
                    else:
                        SeqIO.write(record,ofile,"fasta")
                ofile.close()
    print("create done, save in ",output_for_non_amr_kmer)
    #counting kmer with amr fasta
    print("counting kmer")

    countingKmer(output_for_non_amr_kmer+"/fasta",output_for_non_amr_kmer+"/kmerdump",k)
    print("counting done, create bag of kmer for non amr genes")
    kmerBag=makeKmerBag(k,output_for_non_amr_kmer+"/kmerdump",output_for_non_amr_kmer+'/non_amr_kmer_index.txt')
    nonamrkmer_df=pd.DataFrame({"kmer":list(kmerBag.keys())});
    print("generate ml input for all strain")
    for root, dirs, files in os.walk(output_for_non_amr_kmer+"/kmerdump"):
        for _file in files:
            if _file.endswith('.kmrs'):
                nonamrkmer_df=nonamrkmer_df.merge(countKmerOnFile(str(root)+'/'+_file,kmerBag), how='left', on='kmer')
    nonamrkmer_df.to_csv(output_for_non_amr_kmer+'/non_amr_genome'+'.kmer.csv')
def countKmerOnFile(kmcfile,featureHash):
    #print("Call function generateSample")
    #if not os.path.exists(ml_input):
    #    os.makedirs(ml_input)
    f = open(kmcfile)
    print(kmcfile)
    #get genomeid from kmcfile __name__
    genomeid=os.path.splitext(os.path.basename(kmcfile))[0]
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
    frameCount=pd.DataFrame({genomeid:count})
    frameKmer=frameKmer.join(frameCount)
    return frameKmer
    #frameKmer.to_csv(ml_input+'/'+os.path.splitext(os.path.basename(kmcfile))[0]+'.kmer.csv')
def makeKmerBag(k,kmc_folder,index_file_output):
    kmerHash = {}
    bases=['A','C','G','T']
    #print('k=',k)
    kmers=[''.join(p) for p in itertools.product(bases, repeat=int(k))]
    for i in range(len(kmers)):
        kmerHash[kmers[i]]=0
    #mark the kmer if it exist in genome collection
    for root, dirs, files in os.walk(kmc_folder):
        for _file in files:
            if _file.endswith('.kmrs'):
                f = open(str(root)+'/'+_file)
                for i in f:
                    i = i.strip().split('\t')
                    kmerHash[i[0]] = 1
                f.close()
    file1 = open(index_file_output,"w")
    kmerBag={}
    j=0
    for i in kmerHash:
        if kmerHash[i]==1:
            j+=1
            kmerBag[i]=j
            file1.write(str(i)+'\t'+str(j)+'\n')
    file1.close()
    return kmerBag;
def makeDictForAMRGene(rgi_out):
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
    return list_amr_by_strains
def main():
    generateKmerAllInput(FASTA_IN,output,k)
    generateKmerAMRInput(PROKKA_OUT,RGI_OUT,output,k)
    generateKmerNonAMRInput(PROKKA_OUT,RGI_OUT,output,k)
if __name__== "__main__" :
    main()
