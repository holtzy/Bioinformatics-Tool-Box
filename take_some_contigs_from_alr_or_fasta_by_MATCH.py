#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON Take_some_contigs_from_alr_or_fasta_by_MATCH.py
#     				Prend un fichier fasta ou un fichier alr, et récupère les contigs qui matchent sur un élément de la liste
#  					Yan Holtz, yan1166@hotmail.com
#-------------------------------------

import os
import sys
import re
try:
  import argparse
except ImportError:
  print"oops, the import /argparse/ didn't work"


parser = argparse.ArgumentParser(description= 'Prend un fichier fasta ou un fichier alr, et récupère les contigs qui matchent sur un élément de la liste')
parser.add_argument('-fic', required=True, help=' fichier fasta ou fichier alr dans lequel on veut récupérer un contig ')
parser.add_argument('-liste', required=True, help=' liste des éléments sur lesquels les contigs peuvent matcher')
parser.add_argument('-out', required=True, help=' fichier de sortie')

args = parser.parse_args()

fic=args.fic
liste=args.liste
output=args.out

import re


#------------------------------------------------------------------------------------------------------

#Je créé une liste contenant les séquences que je veux garder.
seq_to_keep=dict()
for line in open(liste):	
	line=line.strip()
	line=line.replace(">","")
	seq_to_keep[line]=""
#------------------------------------------------------------------------------------------------------



	
	
	
#------------------------------------------------------------------------------------------------------

#Et je conserve uniquement les contigs présents dans la liste
tmp=open(output,"w")

for line in open(fic):
	line=line.strip()

	#Si la ligne commence par un ">", alors je vais checker si un élément de ma liste va gréper.
	if line.startswith(">"):
		toprint="no"
		for element in seq_to_keep:
			if re.search(element , line) is not None :
				toprint="yes"

		
	if toprint=="yes":
		tmp.write(line+"\n")
		
tmp.close



print "\n\nOk, c'est gagné, les contigs ont été récupérés !\n\n"

#------------------------------------------------------------------------------------------------------











