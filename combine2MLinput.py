import os, shutil, glob
import csv
import pandas as pd
import json
from sys import argv,stderr
#python combine2MLinput.py [roary_out] [RGI_PATH] [K-mer dump folder] [ArrIndr] [ML_INPUT]

ROARY_OUT=argv[1]
RGI_OUT=argv[2]
KMER_DUMP=argv[3]
arrIndr=argv[4]
ML_INPUT=argv[5]
def generatePangenomeInput(roary_out,rgi_out,mlinput):
    #open file gene_presence_absence.csv, make input files for all genes and accessory genes
    strain=[]
    data = pd.read_csv(roary_out+'/gene_presence_absence.csv', sep=',',na_filter=False, error_bad_lines=False, index_col=False, dtype='unicode')
    strain = data.columns[14:]
    print(strain)
    export_fullgenes_tb = pd.DataFrame({'Strains': strain})
    export_accessory_tb = pd.DataFrame({'Strains': strain})
    export_amr_tb = pd.DataFrame({'Strains': strain})
    export_amr_accessory_tb = pd.DataFrame({'Strains': strain})
    #print(strain)
    #open cgi files, collect AMR gene into a set and make a hash map
    set_of_amr_aro=set()
    hashmap={}
    for root, dirs, files in os.walk(rgi_out):
        for _file in files:
            if _file.endswith('.json'):
                with open (rgi_out+'/'+_file, 'r') as f:
                    jsonreader =json.load(f)
                    basename=os.path.basename(_file)
                    strainname = basename.rsplit('.', 1)[0]
                    for item in jsonreader.values():
                        if isinstance(item,dict):
                            for k, v in item.items():
                                if isinstance(v,dict):
                                    set_of_amr_aro.add(v['ARO_name'].strip().lower())
                                    hashmap[strainname+v['ARO_name'].strip().lower()]=1
    print('number of arm genes:',len(set_of_amr_aro))
    print(set_of_amr_aro)
    #init amr-genes table
    #print (strain)
    for gene in set_of_amr_aro:
        genes=[]
        for s in strain:
            key=s+gene;
            #print(s)
            if key in hashmap.keys() :
                #print(hashmap[key])
                genes.append(1)
                #print('add 0')
            else:
                #print('add 1')
                genes.append(0)
        newc= pd.DataFrame({gene:genes})
        export_amr_tb=export_amr_tb.join(newc)
    print('export_amr_tb:',export_amr_tb)
    #generate export files for whole genes and accessory gene
    print('process pangenome')
    for index, row in data.iterrows():
        if(index>0):
            genes=[]
            #print('process at ',index)
            for s in strain:
                #print(s,":",row[s])
                if(row[s]==''):
                    genes.append(0)
                    #print('add 0')
                else:
                    #print('add 1')
                    genes.append(1)
            #print('new gen:',gene)
            newc= pd.DataFrame({row['Gene']: genes})
            #print('new gen:',newc)
            export_fullgenes_tb=export_fullgenes_tb.join(newc)
            if(row['Accessory Fragment']!=''):
                export_accessory_tb=export_accessory_tb.join(newc)
                isFound=False
                for gene in set_of_amr_aro:
                    if row['Gene'].strip().lower() ==gene or gene.find(row['Gene'].strip().lower()) != -1 :
                        #print("found in amr genes:",row['Gene'])
                        isFound=True
                if isFound:
                    export_amr_accessory_tb=export_amr_accessory_tb.join(newc)
    print('export_accessory_tb:',export_accessory_tb)
    print('export_amr_accessory_tb:',export_amr_accessory_tb)

    export_fullgenes_tb.to_csv(mlinput+'/AllGenes.csv')
    export_accessory_tb.to_csv(mlinput+'/AccessoryGenes.csv')
    export_amr_tb.to_csv(mlinput+'/AMRGenes.csv')

def generateKmerInput(kmerdump,ArrIndr,ml_input):
    #generate all kmer
    print("Generate k-mer files")
    if not os.path.exists(ml_input+"/kmer"):
        os.makedirs(ml_input+"/kmer")
    f = open(ArrIndr)
    # init feature hash
    featureHash = {}
    # for each line in file, get the feature and index
    # and set it in the hash
    for i in f:
        i = i.strip().split('\t')
        if len(i) < 2:
            continue
        featureHash[i[1]] = i[0]
    f.close()
    # open up the contigs file
    for root, dirs, files in os.walk(kmerdump):
        for _file in files:
            if _file.endswith('.kmrs'):
                generateSample(str(root)+'/'+_file,featureHash,ml_input)

    # create array to hold contigs

def generateSample(kmcfile,featureHash,ml_input):
    #print("Call function generateSample")
    f = open(kmcfile)

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
        kmer.append(featureHash[i])
        if featureHash[i] in kmcContigs:
            count.append(kmcContigs[featureHash[i]])
        else:
            count.append(0)
    frameKmer= pd.DataFrame({"kmer":kmer})
    frameCount=pd.DataFrame({"count":count})
    frameKmer=frameKmer.join(frameCount)
    frameKmer.to_csv(ml_input+'/kmer/'+os.path.splitext(os.path.basename(kmcfile))[0]+'.kmer.csv')
def main():
    generatePangenomeInput(ROARY_OUT,RGI_OUT,ML_INPUT)
    generateKmerInput(KMER_DUMP,arrIndr,ML_INPUT)
if __name__== "__main__" :
    main()
