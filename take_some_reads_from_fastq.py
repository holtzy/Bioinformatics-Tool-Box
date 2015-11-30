#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON Take_some_reads_from_fastq
#     				Read a Fastq File and keep only reads specified in a list in the output
#  					Yan Holtz, yan1166@hotmail.com
#-------------------------------------

import os
import sys
import re
try:
  import argparse
except ImportError:
  print"oops, the import /argparse/ didn't work"


parser = argparse.ArgumentParser(description= '\n\nRead a Fastq File and keep only reads specified in a list in the output\n\n')
parser.add_argument('-fic', required=True, help=' FastQ file you want to filter ')
parser.add_argument('-liste', required=True, help=' List of reads you want to keep')
parser.add_argument('-out', required=True, help=' output')

args = parser.parse_args()

fic=args.fic
liste=args.liste
output=args.out




#------------------------------------------------------------------------------------------------------

#Je créé une liste contenant les séquences que je veux garder.
seq_to_keep=dict()
for line in open(liste):	
	line=line.strip()
	seq_to_keep[line]=""	
#------------------------------------------------------------------------------------------------------



#------------------------------------------------------------------------------------------------------

#Et je conserve uniquement les contigs présents dans la liste
tmp=open(output,"w")
num=3

for line in open(fic):
	num+=1
	if num == 4:
		num=0
		line=line.strip()
		if line.startswith("@"):
			print line
			contig=line.replace("@","")
			contig=contig.split()[0]
		if contig in  seq_to_keep:
			toprint="yes"
		else:
			toprint="no"
		
	if toprint=="yes":
		line=line.replace("\n","")
		tmp.write(line+"\n")
		
tmp.close



print "\n\nOk, c'est gagné, les reads desires  ont été récupérés !\n\n"
#------------------------------------------------------------------------------------------------------












