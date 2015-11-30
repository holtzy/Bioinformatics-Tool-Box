#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON Make_my_PFAS
#     				This script permits to build a pfas file from a fasta + a SNP file
#  					Yan Holtz, yan1166@hotmail.com
#-------------------------------------

import os
import sys
import re
try:
  import argparse
except ImportError:
  print"oops, the import /argparse/ didn't work"


parser = argparse.ArgumentParser(description= 'This script permits to build a pfas file from a fasta + a SNP file')
parser.add_argument('-fasta', required=True, help='fasta file')
parser.add_argument('-SNP', required=True, help='SNP file file')
parser.add_argument('-all', required=False, help='yes or no, do you want all contigs in the pfas, or only the ones containing a SNP ?')
parser.add_argument('-out', required=True, help='SNP file file')


args = parser.parse_args()

fasta=args.fasta
SNP=args.SNP
all=args.all
out=args.out

#Je vais écrire dans mon fichier de sortie
tmp=open(out , "w")


#------------------------------------------------------------------------------------------------------
print "\n\n--------\n\nBienvenu ! C'est cool un peu d'activité !\n\n------------\n\n"		
#------------------------------------------------------------------------------------------------------




#------------------------------------------------------------------------------------------------------
print "\n\n--------\n\nETAPE 1 : dictionnaire du fasta mon ami ! \n\n------------\n\n"		

dico_fasta=dict() 
nombre=0

for line in open(fasta):
	line=line.strip()
	if line.startswith(">"):
		nombre=nombre+1
		contig=line.replace(">","")
		if contig.startswith("Traes"):
			contig=contig.split("|")[0]
		dico_fasta[contig]=""
		
	else:
		dico_fasta[contig]=dico_fasta[contig]+line
#------------------------------------------------------------------------------------------------------







#------------------------------------------------------------------------------------------------------
print "\n\n--------\n\nETAPE 2 : lecture des SNPs mon ami ! \n\n------------\n\n"		

nombre=0 
dico_pfas=dict()
dico_contigs_with_SNP=dict()

for line in open(SNP):
	nombre=nombre+1
	line=line.strip()
	line=line.split()
	contig=line[6].replace(">","")
	if contig.startswith("Traes"):
		contig=contig.split("|")[0]
	position=int(line[7]) 
	
	sequence=dico_fasta[contig]
	dico_contigs_with_SNP[contig]="yes"
	
	#Maintenant il va falloir créer une séquence pour chaque allele du SNP
	ind=0
	for i in range(9,len(line),2):
		
		#Quelle est le num de l'indiv que j'étudie ?
		ind+=1 
		
		#Quelles sont ses 2 allèles ? Attention, en cas de données manquantes, je met un "N"
		if line[i][0] == "-" :
			all1="-" ; all2="-"
		else:
			all1=line[i][0]
			all2=line[i][1]
		
		#Quelle va être le nom de la séquence de sortie ?
		name1=contig+"|ind_"+str(ind)+"|Allele_1"
		name2=contig+"|ind_"+str(ind)+"|Allele_2"
		
		#A t'on déja comencé à taffer sur ce contig ?
		#Si c'est le cas, on récupère la séquence que l'on a déja commencé à modifier 
		if name1 in dico_pfas :
			new_seq_allele1=dico_pfas[name1]
			new_seq_allele2=dico_pfas[name2]
		#Sinon on prend la séquence du fasta d'origine
		else :
			new_seq_allele1=sequence 
			new_seq_allele2=sequence 
		
		#Je transfome ma séquence en liste pour pouvoir transformer un élément.
		new_seq_allele1=list(new_seq_allele1) 
		new_seq_allele2=list(new_seq_allele2)
		
		#Et puis on rajoute les données du SNP actuel à cette séquence
		new_seq_allele1[position - 1] = all1
		new_seq_allele2[position - 1] = all2
		
		#Et on retransforme en str la séquence
		new_seq_allele1="".join(new_seq_allele1)
		new_seq_allele2="".join(new_seq_allele2)

		#Et on réactualise le dico pfas
		dico_pfas[name1]=new_seq_allele1
		dico_pfas[name2]=new_seq_allele2

print "nbr de SNPs détectés = "+str(nombre)
print "nbr de génotypes dans le fichier de snp = "+str(ind)


#------------------------------------------------------------------------------------------------------




#------------------------------------------------------------------------------------------------------
# OPTION ALL : Si je veux afficher aussi les contigs sans SNP dans le pfas de sortie

if all=="yes":
	for contig in dico_fasta:
		if contig not in dico_contigs_with_SNP:
			for ind in range(1,ind+1):
				name1=contig+"|ind_"+str(ind)+"|Allele_1"
				name2=contig+"|ind_"+str(ind)+"|Allele_2"
				dico_pfas[name1]=dico_fasta[contig]
				dico_pfas[name2]=dico_fasta[contig]
	



#------------------------------------------------------------------------------------------------------
print "\n\n--------\n\nETAPE 3 : impression du fichier de sortie \n\n------------\n\n"		

#Pour la sortie, je vais essayer d'obtenir un pfas ordonné : toutes les séquences d'un contigs groupées.

for contig in dico_fasta:
	for ind in range(1,ind+1):
		name1=contig+"|ind_"+str(ind)+"|Allele_1"
		name2=contig+"|ind_"+str(ind)+"|Allele_2"
		
		if name1 in dico_pfas : 
			to_print= ">"+name1+"\n"+dico_pfas[name1]+"\n"+">"+name2+"\n"+dico_pfas[name2]+"\n"
			tmp.write(to_print)

