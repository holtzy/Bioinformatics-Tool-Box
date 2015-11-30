#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON Calculte_GC_content.py
#     				Ce script permet de calculer le taux de GC de toutes les séquences d'un fichier Fasta
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


parser = argparse.ArgumentParser(description= 'get the baytes from durum wheat SNPs')
parser.add_argument('-fasta', required=True, help='	fasta de EPO')
parser.add_argument('-out', required=True, help=' fichier de sortie')


args = parser.parse_args()
fasta=args.fasta
out=args.out




tmp=open(out,"w")




#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### Calcul pour chaque séquence



for line in open(fasta):
	line=line.strip()
	
	if line.startswith(">"):
		name=line.replace(">","")
	
	else :
		tot=0 ; nb=0
		for i in line :
			tot+=1
			if i=="C" or i=="G" :
				nb+=1
					
		taux=float(nb)/float(tot)*100
		
		tmp.write(name+"\t"+str(taux)+"\n")
	
				
				
				
				
				
				
				
				












