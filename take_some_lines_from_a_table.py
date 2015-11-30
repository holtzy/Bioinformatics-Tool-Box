#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON Take_some_lines_from_a_table
#     				Prend un fichier contenant un tableau, et récupère toutes les lignes pour lesquels les champs d'une colonne donnée sont présent dans une liste
#  					Yan Holtz, yan1166@hotmail.com
#-------------------------------------


import os
import sys
import re
try:
  import argparse
except ImportError:
  print"oops, the import /argparse/ didn't work"


parser = argparse.ArgumentParser(description= 'permet de récupérer dans un tableau, toutes les lignes dont une colonne précise est présente dans un fichier liste.')
parser.add_argument('-fic', required=True, help=' fichier fasta ou fichier alr dans lequel on veut récupérer un contig ')
parser.add_argument('-column', required=True, help=' colonne dans laquelle on cherche un élément de la liste')
parser.add_argument('-liste', required=True, help=' liste des contigs à récupérer')
parser.add_argument('-out', required=True, help=' liste des contigs à récupérer')

args = parser.parse_args()

fic=args.fic
column=args.column
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
column=int(column)

for line in open(fic):
	line=line.strip()
	save=line
	line=line.split()
	contig=line[column-1]
	
	if contig.startswith(">"):
		contig=contig.replace(">","")

	if contig in  seq_to_keep:
		tmp.write(save+"\n")
		
tmp.close



print "\n\nOk, c'est gagné, les entrées recherchées ont été récupérés !\n\n"
#------------------------------------------------------------------------------------------------------












