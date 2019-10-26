import os, shutil, glob
import csv
import pandas as pd
import json
import itertools
import gffutils
import pyfaidx
from sys import argv,stderr
#python combine2MLinput_GeneCluster.py [PROKKA_OUT] [roary_out] [RGI_PATH] [ML_INPUT]
PROKKA_OUT=argv[1]
ROARY_OUT=argv[2]
RGI_OUT=argv[3]
ML_INPUT=argv[4]

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
def main():
    generatePangenomeInput(ROARY_OUT,RGI_OUT,ML_INPUT)
if __name__== "__main__" :
    main()
