#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON Delete_redondant_sequences_in_fasta.py
#     				Take a Fasta OR ALR file and delete sequences that are present several times
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
parser.add_argument('-out', required=True, help=' Output')

args = parser.parse_args()

fic=args.fic
output=args.out

	
	
	
#------------------------------------------------------------------------------------------------------

Mydico=dict()
tmp=open(output,"w")
toprint="no"

for line in open(fic):
	line=line.strip()
	if line.startswith(">"):
		if line not in Mydico :
			toprint="yes"
			Mydico[line]="yes"
			tmp.write(line+"\n")
		else:
			toprint="no"
	else:
		if toprint=="yes" :
			tmp.write(line+"\n")
	
	

tmp.close



print "\n\nOk, c'est gagné, les contigs ont été récupérés et les contigs redondant ont été supprimés !! \n\n"
#------------------------------------------------------------------------------------------------------












