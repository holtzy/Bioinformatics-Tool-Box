#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON Merge_IDXSTAT_Files
#     				You have several samtools IDXstats ? (e.g. one per individual) --> Merge them together !
#
# 					Yan Holtz, yan1166@hotmail.com
#-------------------------------------

import os
import sys
import re
try:
  import argparse
except ImportError:
  print"oops, the import /argparse/ didn't work"


parser = argparse.ArgumentParser(description= '\n\nYou have several samtools IDXstats ? (e.g. one per individual) --> Merge them together !\n\n')
parser.add_argument('-list', required=True, help='list of IDXstat files to merge together')


args = parser.parse_args()

list=args.list




#------------------------------------------------------------------------------------------------------

print "\n\n--------\n\nBienvenu ! C'est cool un peu d'activité !\n\n------------\n\n"		

dico=dict()

#pour l'initialiation :
deb="yes"
entete="contig"

for fichier in open(list):
	fichier=fichier.strip()
	print fichier
	line_number=0
	name=fichier.replace("idxstats_resultat_mapping_sur_EPO_","")
	name=name.replace(".bam","")
	entete=entete+"\t"+name
	
	
	#Initialisation = on rempli les noms du dico
	if deb=="yes" :	
		for line in open(fichier):
			line=line.split()
			contig=line[0]
			dico[contig]=""
		deb="no"
			
	#Ensuite on rempli pour chaque fichier
	for line in open(fichier):
		line=line.split()
		contig=line[0]
		nbr_reads=line[2]
		dico[contig]=str(dico[contig])+"\t"+str(nbr_reads)
			
#Puis j'imprime
tmp=open("Bilan_idxstats.txt","w")
tmp.write(entete+"\n")
for contig in dico:
	tmp.write(contig+dico[contig]+"\n")
tmp.close
		
		



print "\n\n--------\n\nOK, c'est gagné !!!\n\n------------\n\n"		

#------------------------------------------------------------------------------------------------------

