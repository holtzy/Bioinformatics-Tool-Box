#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON SPLIT_FASTA_BY_CONTIGS
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
parser.add_argument('-keepOthers', action='store_true', help="if used sequenced not in the list will be kept", default=False)

args = parser.parse_args()

fic=args.fic
liste=args.liste
output=args.out
keepOthers=args.keepOthers

#------------------------------------------------------------------------------------------------------

#Je créé une liste contenant les séquences que je veux garder.
seq_to_keep=set()
for line in open(liste):	
	line=line.strip()
	line=line.replace(">","")
	seq_to_keep.add(line)
	
#------------------------------------------------------------------------------------------------------

#Et je conserve uniquement les contigs présents dans la liste
tmp=open(output,"w")

for line in open(fic):
	line=line.strip()
	if line.startswith(">"):
		contig=line.replace(">","")
	toprint="no"
	if contig in  seq_to_keep:
		if keepOthers==False :
			toprint="yes"
	else:
		if keepOthers==True :
			toprint="yes"
		
	if toprint=="yes":
		tmp.write(line+"\n")
		
tmp.close



print "les contigs ont été récupérés !\n"












