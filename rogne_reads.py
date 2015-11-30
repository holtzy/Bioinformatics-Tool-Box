#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON Rogne_read.py
#     				Permet de couper les x premieres bases de tous les reads d'un fastq
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


parser = argparse.ArgumentParser(description= '\n\npermet de rogner les x premieres bases de tous les reads d\'un fichier fastq\n\n')
parser.add_argument('-fastq', required=True, help='	fichier fastq')
parser.add_argument('-nbr', required=False, help='	nbr de bases que l\'on veut rogner au DEBUT de la séquence')
parser.add_argument('-out', required=False, help='	fichier de sortie')

args = parser.parse_args()
fic_fastq=args.fastq
nbr=args.nbr
out=args.out





#----------------------------------------------------------------------------------------------------------------------------------------------
tmp=open(out,"w")
num=0
nbr=int(nbr)

#En fait le jeu c'est de rogner une ligne sur 2.
for line in open(fic_fastq) :
	line=line.strip()
	num+=1
	
	if num==2 or num==4 :
		line=line[nbr:]
	
	tmp.write(line+"\n")
	
	if num==4 :
		num=0	

tmp.close()
#----------------------------------------------------------------------------------------------------------------------------------------------













