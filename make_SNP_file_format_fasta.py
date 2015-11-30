#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON get_baytes
#     				Ce script permet de mettre les SNPs au format fasta avec des [] au niveau des SNPs
#				
#  					Yan Holtz, yan1166@hotmail.com
#-------------------------------------


import os
import sys
import re
try:
  import argparse
except ImportError:
  print"oops, the import /argparse/ didn't work"


parser = argparse.ArgumentParser(description= '\n\nCe script permet de mettre les SNPs au format fasta avec des [] au niveau des SNPs\n\n')
parser.add_argument('-fasta1', required=True, help='	fasta de EPO')
parser.add_argument('-fasta2', required=True, help='	fasta de EPO')
parser.add_argument('-SNP', required=True, help='	fichier de SNP')
parser.add_argument('-out', required=True, help='	fichier de sortie')


args = parser.parse_args()
fasta1=args.fasta1
fasta2=args.fasta2
SNP=args.SNP
output=args.out







#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 1 : Collecte des infos des SNPs
print "\n\nSTEP 1 : dictionnaire des positions des SNPs"
nombre=0 
tmp=open(output,"w")
dico_SNP=dict()


for line in open(SNP):
	nombre=nombre+1
	line=line.strip()
	line=line.split()
	contig=line[6].replace(">","")
	if contig.startswith("Traes"):
		contig=contig.split("|")[0]
	position=int(line[7]) 
	SNP_name=contig+"@"+str(position)
	
	#Maintenant il va falloir récupérer les 2 allèles du SNP... pas facile.
	allele1=""
	allele2=""
	first = "yes"
	for i in range(9,len(line),2):
		all=line[i][0]
		if all == "-" :
			continue
		if first == "yes" :
			allele1=all
			first="no"
		if all != allele1 :
			allele2=all
	
	to_add=	str(position)+"\t"+allele1+"\t"+allele2
	if contig not in dico_SNP:
		dico_SNP[contig]=list()
	dico_SNP[contig].append(to_add)
	

print "nbr de SNPs détectés = "+str(nombre)
	
	
	
	
	

#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 2 : Collecte de l'ensemble des contigs blé tendre dans un dictionnaire
print "\n\nSTEP 2 : dictionnaire des contigs de blé tendre"
dico_EPO=dict() 
nombre=0

for line in open(fasta1):
	line=line.strip()
	if line.startswith(">"):
		nombre=nombre+1
		contig=line.replace(">","")
		if contig.startswith("Traes"):
			contig=contig.split("|")[0]
		dico_EPO[contig]=""
		
	else:
		dico_EPO[contig]=dico_EPO[contig]+line


print "\n\n      ---> done !     "+str(nombre)+"     Contigs ont été récuré dans le fichier fasta"





#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 2 PRIM: Collecte de l'ensemble des contigs blé tendre dans un dictionnaire
print "\n\nSTEP 2 : dictionnaire des contigs de blé tendre"
nombre=0

for line in open(fasta2):
	line=line.strip()
	if line.startswith(">"):
		nombre=nombre+1
		contig=line.replace(">","")
		if contig.startswith("Traes"):
			contig=contig.split("|")[0]
		dico_EPO[contig]=""
		
	else:
		dico_EPO[contig]=dico_EPO[contig]+line


print "\n\n      ---> done !     "+str(nombre)+"     Contigs ont été récuré dans le fichier fasta numéro 2"


#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------










#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 3 : modification des fasta

name=output+"_short"
tmp2=open(name,"w")
nbr_snp=0

print "\n\nSTEP 3 : modification des fasta"

for contig in dico_EPO:
	sequence=dico_EPO[contig]

	#Si je n'ai pas de SNPs sur le contig
	if contig not in dico_SNP :
		tmp.write(">"+contig+"\n"+sequence+"\n")

	#je modifie le fasta SNP par SNP quand j'ai des SNPs dessus
	sequence=list(sequence)
	if contig in dico_SNP :
		SNP_name=contig
		for snp in range(0,len(dico_SNP[contig]))  :
		
			nbr_snp+=1
			
			#Je récuprere la position du SNP et je l'ajoute au nom du SNP
			position = int(dico_SNP[contig][snp-1].split()[0] )
			SNP_name=SNP_name+"@"+str(position)
			
			#Je récupere les 2 alleles et je modifie le fasta d origine.
			allele1 = dico_SNP[contig][snp-1].split()[1] 
			allele2 = dico_SNP[contig][snp-1].split()[2] 
			
			sequence[position-1]="["+str(allele1)+"/"+str(allele2)+"]"
			
			#Je sors un autre fasta, avec une séquence de 100 base par SNP entourant le SNP.
			seq_short = list(dico_EPO[contig])
			seq_short[position-1]="["+str(allele1)+"/"+str(allele2)+"]"
			seq_short = seq_short[position-50 : position+49]
			seq_short="".join(seq_short)
			toprint=">"+contig+"@"+str(position)+"\n"+seq_short+"\n"
			tmp2.write(toprint)

		#Maintenant que j ai traité tous les SNP, je peux l'imprimer dans mon fichier de sortie.
		sequence="".join(sequence)
		tmp.write(">"+SNP_name+"\n"+sequence+"\n")	
	
	
				


print "\n\n      ---> done !    "+str(nbr_snp)+"   SNPs ont été testés et ajouté!!"









