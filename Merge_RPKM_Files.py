#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON Merge_RPKM_Files
#     				permet de merger des tableaux contenant les données RPKM
#  					Yan Holtz, yan1166@hotmail.com
#-------------------------------------

import os
import sys
import re
try:
  import argparse
except ImportError:
  print"oops, the import /argparse/ didn't work"


parser = argparse.ArgumentParser(description= 'permet de couper un pfas contig par contig')
parser.add_argument('-list', required=True, help='list of RPKM files to merge')


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
	name=fichier.replace("RPKM_","")
	entete=entete+"\t"+name
	
	
	#Initialisation = on rempli les noms du dico
	if deb=="yes" :	
		for line in open(fichier):
			line_number+=1
			if line_number>2:
				line=line.split()
				contig=line[0]
				dico[contig]=""
		deb="no"
		line_number=0
			
	#Ensuite on rempli pour chaque fichier
	for line in open(fichier):
		line_number+=1
		if line_number>2:
			line=line.split()
			contig=line[0]
			RPKM=line[7]
			dico[contig]=str(dico[contig])+"\t"+str(RPKM)
			
#Puis j'imprime
tmp=open("Bilan_RPKM.txt","w")
tmp.write(entete+"\n")
for contig in dico:
	tmp.write(contig+dico[contig]+"\n")
tmp.close
		
		



print "\n\n--------\n\nOK, c'est gagné !!!\n\n------------\n\n"		
	
		
		
		
		
		
		
		
		
		
		
		
