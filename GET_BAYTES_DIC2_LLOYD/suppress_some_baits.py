#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT Tranform_bayt_format
#     				Pour transformer les baytes comme les veux sylvain santoni
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


parser = argparse.ArgumentParser(description= 'Va résumer toutes les infos de divers fichiers pour chaque contig EPO.')
parser.add_argument('-baits', required=True, help=' fichier contenant les baytes provenant de Dic2 * Silur')
parser.add_argument('-dell', required=False, help=' fichier contenant les baytes provenant de EPO')
parser.add_argument('-out', required=False, help=' fichier contenant les baytes provenant de EPO')


args = parser.parse_args()
baits=args.baits
dell=args.dell
out=args.out




#DICO DES sequences a supprimer
dico_del=dict()
for i in open(dell) :
	i=i.strip()
	dico_del[i]="-"

	


tmp=open(out , "w")
nb_line=0

for line in open(baits):
	
	line=line.strip()
	if line.startswith(">"):
		snp=line.split("|All")[0]
		snp=snp.replace(">","")
		
		if snp in dico_del:
			toprint="no"
		else:
			toprint="yes"
			
	if toprint == "yes" and nb_line<40000:
		tmp.write(line+"\n")
		nb_line+=1
		
		
		
		
		
		









