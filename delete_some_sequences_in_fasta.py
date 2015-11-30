#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON Delete_Some_sequences_in_fasta
#     				Take a Fasta OR ALR file and delete chosen sequences that are present in a given list
#  					Yan Holtz, yan1166@hotmail.com
#-------------------------------------

import os
import sys
import re
try:
  import argparse
except ImportError:
  print"oops, the import /argparse/ didn't work"


parser = argparse.ArgumentParser(description= '\n\nTake a Fasta OR ALR file and delete sequences that are present several times\n\n')
parser.add_argument('-fic', required=True, help=' Input : ALR or FASTA file ')
parser.add_argument('-list', required=True, help=' list of sequence to delete')
parser.add_argument('-out', required=True, help=' Output')

args = parser.parse_args()

fic=args.fic
list=args.list
output=args.out

	
	
	
#------------------------------------------------------------------------------------------------------
# DICTIONNAIRY OF SEQUENCES TO DELETE

dico_to_del=dict()

for line in open(list):
	line=line.strip()
	line=line.replace(">","")
	dico_to_del[line]="deleteme"

#------------------------------------------------------------------------------------------------------






#------------------------------------------------------------------------------------------------------
#LECTURE DU FASTA

tmp=open(output,"w")
toprint="no"

for line in open(fic):
	line=line.strip()
	if line.startswith(">"):
	
		line=line.replace(">","")
		if line not in dico_to_del :
			toprint="yes"
			tmp.write(">"+line+"\n")
		else:
			toprint="no"
	else:
		if toprint=="yes" :
			tmp.write(line+"\n")
	
	

tmp.close



print "\n\nOk, c'est gagné, les contigs ont été récupérés et les contigs indiqués ont été supprimés !! \n\n"
#------------------------------------------------------------------------------------------------------



