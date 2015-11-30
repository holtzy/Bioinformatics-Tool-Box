#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON 
#     				A partir des fichiers Info_SNP, je vais trouver et sortir le meilleur SNP de chaque contig et le baytes qui va avec
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


parser = argparse.ArgumentParser(description= 'Va résumer toutes les infos de divers fichiers pour chaque contig EPO.')
parser.add_argument('-ecart', required=True, help=' ecart minimum entre 2 SNPs')
parser.add_argument('-baytes', required=True, help=' fichier contenant la totalité des baytes ')
parser.add_argument('-out', required=True, help=' fichier Info_SNP, dis notamment si le SNP a un baytes ou non. ')

args = parser.parse_args()
ecart=args.ecart
baytes=args.baytes
out=args.out








#-------------------------------
#PREMIERE LECTURE

dico_contig=dict()
dico_snp_to_print=dict()

#Je lis mes baytes un par un.
for line in open(baytes):

	if line.startswith(">"):
		line=line.replace(">","")
		snp=line.split("|All")[0]
		position=snp.split("@")[1]
		contig=snp.split("@")[0]

		#de base je garde le SNP, sauf si il respecte certaines conditions auquel cas je le laisse tomber
		to_keep="yes"
		
		#Si ce SNP a déja été rencontré, on ne le traitera pas une deuxième fois (allèle Dic2 uis Silur)
		if snp in dico_snp_to_print:
			continue
		
		#Si le contig a déja 2 SNP, je n'en prendrai pas un de plus
		if contig in dico_contig:
			nb_snp=len(dico_contig[contig].split(","))
			if nb_snp >= 2 :
				to_keep="no"	

		#Si le contig a déja un seul SNP, mais que le SNP actuel est trop proche de celui, alors on laisse tomber
		if contig in dico_contig:
			nb_snp=len(dico_contig[contig].split(","))
			if nb_snp == 1:
				pos_snp1=int(dico_contig[contig])
				if int(position) <= (int(pos_snp1) + int(ecart)) and int(position) >= (int(pos_snp1) - int(ecart))  :
					to_keep="no"
					
					
					

					
					
		#J'enregistre les SNPs sélectionnés dans les 2 dicos
		if to_keep=="yes" :
			
			#dico des snps sélectionnés
			dico_snp_to_print[snp]="ok"
			
			#dico repertoriant les snps de chaque contig
			if contig not in dico_contig:
				dico_contig[contig]=position
			else:
				dico_contig[contig]=dico_contig[contig]+","+str(position)






#-------------------------------
#récupération des baytes séectionnés

tmp=open(out,"w")

for line in open(baytes):

	if line.startswith(">"):
		__line=line.replace(">","")
		snp=__line.split("|All")[0]


	if snp in dico_snp_to_print:
		to_print="yes"
	else:
		to_print="no"
		
		
	if to_print=="yes":
		tmp.write(line)






	

























