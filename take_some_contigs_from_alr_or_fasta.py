#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON Take_some_contigs_from_ALR_or_Fasta
#     				Prend un fichier fasta ou un fichier alr, et récupère les contigs que l'on place dans une liste.
#  					Yan Holtz, yan1166@hotmail.com
#-------------------------------------

import os
import sys
import re
try:
  import argparse
except ImportError:
  print"oops, the import /argparse/ didn't work"


parser = argparse.ArgumentParser(description= 'permet d\'utiliser les infos d\'un fichier .alr et d\'un .gen pour détecter des snps.')
parser.add_argument('-fic', required=True, help=' fichier fasta ou fichier alr dans lequel on veut récupérer un contig ')
parser.add_argument('-liste', required=True, help=' liste des contigs à récupérer')
parser.add_argument('-out', required=True, help=' liste des contigs à récupérer')

args = parser.parse_args()

fic=args.fic
liste=args.liste
output=args.out




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
	if line.startswith(">"):
		contig=line.replace(">","")
	if contig in  seq_to_keep:
		toprint="yes"
	else:
		toprint="no"
		
	if toprint=="yes":
		tmp.write(line+"\n")
		
tmp.close



print "\n\nOk, c'est gagné, les contigs ont été récupérés !\n\n"

#------------------------------------------------------------------------------------------------------












