import os, shutil, glob
from sys import argv,stderr

#python pipeline.py [fna_folder] [output_folder]

FASTA_IN=argv[1]
output=argv[2]
PROKKA_OUTPUT=output+'/PROKKA_OUT'
GFF_OUTPUT=output+'/PROKKA_GFF'
ROARY_OUTPUT=output+'/ROARY_OUT'
RGI_OUTPUT=output+'/RGI_OUTPUT'
KMC_OUTPUT=output+'/KMC_OUTPUT'
#RGI_PATH=argv[3]
def annotateProtein(fasta_folder,prokka_folder):
    #annotate each fasta file with prokka, the output is stored in the folder with the folder name is fasta filename
    for root, dirs, files in os.walk(fasta_folder):
        for _file in files:
            if _file.endswith(('.fna','.fasta','.fn')):
                # If we find it, notify us about it and copy it it to C:\NewPath\
                print ('Found file in: ' , str(root))
                #strain name is filename
                newname = os.path.splitext(os.path.basename(_file))[0]
                cmd='prokka --kingdom Bacteria --outdir '+prokka_folder+'/'+newname+' --genus Listeria --locustag '+newname+' '+str(root)+'/'+ _file
                os.system(cmd)
def geneClustering(prokka_folder,gff_folder,roary_folder):
    #copy all gff files in output of Prokka into a folder for clustering
    if not os.path.exists(gff_folder):
        os.makedirs(gff_folder)
    for root, dirs, files in os.walk(prokka_folder):
        for _file in files:
            if _file.endswith('.gff'):
              #print ('Found file in: ' , str(root),',',_file)
              newname = os.path.basename(str(root))+'.gff'
              print ('new  file: ' , newname)
              shutil.copy(os.path.abspath(str(root) + '/' + _file), gff_folder+'/'+newname)
    myCmd = 'roary -f '+ROARY_OUTPUT+' -e -n -v '+gff_folder+'/*.gff'
    os.system(myCmd)
def annotateAMRgenes(PROKKA_OUTPUT,RGI_OUTPUT,RGI_PATH):
    #annotate AMR genes base on protein-annotation files of prokka
    for root, dirs, files in os.walk(PROKKA_OUTPUT):
        for _file in files:
            if _file.endswith('.faa'):
              newname =os.path.basename(str(root))
              myCmd = 'rgi main --input_sequence '+PROKKA_OUTPUT+'/'+newname+'.fna/'+_file+' --output_file '+RGI_OUTPUT+'/'+newname+' --input_type protein --local --clean'
              print ('cmd: ' , myCmd)
              os.system(myCmd)
def countingKmer(fasta_folder,kmc_output):
    #generate kmer counting files from fasta files using kmc_dump
    if not os.path.exists(kmc_output):
        os.makedirs(kmc_output)
    for root, dirs, files in os.walk(fasta_folder):
        for _file in files:
            if _file.endswith(('.fna','.fasta','.fn')):
                #print ('Found file in: ' , str(root),_file)
                newname = os.path.splitext(os.path.basename(_file))[0]
                #print ('Basename ' , newname);
                cmd='kmc -k10 -ci1 -fa -fm '+fasta_folder+'/'+_file+ ' '+kmc_output+'/'+newname+' '+ kmc_output+'/'+newname+'_tmp'
                os.system(cmd)
                cmd_dump='kmc_dump '+kmc_output+'/'+newname+' '+kmc_output+'/'+newname+'.10.kmrs'
                os.system(cmd_dump)
                os.remove(kmc_output+'/'+newname+'.kmc_suf')
                os.remove(kmc_output+'/'+newname+'.kmc_pre')

def main():
    annotateProtein(FASTA_IN,PROKKA_OUTPUT)
    geneClustering(PROKKA_OUTPUT,GFF_OUTPUT,ROARY_OUTPUT)
    annotateAMRgenes(PROKKA_OUTPUT,RGI_OUTPUT,RGI_PATH)
if __name__== "__main__" :
    main()
