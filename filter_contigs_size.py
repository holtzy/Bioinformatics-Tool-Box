#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON Filter_contigs_size.py
#     				Take a Fasta file and delete sequences that have a size under a chosen threshold
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


parser = argparse.ArgumentParser(description= '\n\nTake a Fasta file and delete sequences that have a size under a chosen threshold\n\n')
parser.add_argument('-fasta', required=True, help='	Input : fasta file that you want to filter')
parser.add_argument('-size', required=True, help='	limit size you want to keep')
parser.add_argument('-out', required=True, help='	output')


args = parser.parse_args()
fic_fasta=args.fasta
size=args.size
out=args.out







#------------------------------------------------------------------------------------------------------

### STEP 0 : Collecte de l'ensemble des contigs dans un dictionnaire, j'en profite pour récupérer la longueur des contigs
tmp=open(out,"w")
tot=0
tot_pass=0
sequence=""

for line in open(fic_fasta):
	line=line.strip()

	if line.startswith(">"):
		
		#Je finalise le contig précédent.
		if len(sequence) > size :
			tmp.write(contig+"\n"+sequence+"\n")
			tot_pass+=1
		sequence=""
		contig=""
		
		#Et je prépare le contig en cours
		tot=tot+1
		contig=line
		
	else:
		sequence=sequence+line
		


print "\n\n-----------\nInitial number of contigs = "+str(tot)+"\nContigs larger than "+str(size)+ " : "+str(tot_pass)+"\n-----------\n\n"

#------------------------------------------------------------------------------------------------------





		
