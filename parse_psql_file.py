#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT Parse_psql_file.py
#
#     				La sortie de prot4est (fichier .psql) est assez difficile à déchiffrer. Ce script permet de lire ce fichier et de sortir une sortie simple : contig, début CDS, fin CDS, orientation.
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


parser = argparse.ArgumentParser(description= 'La sortie de prot4est (fichier .psql) est assez difficile à déchiffrer. Ce script permet de lire ce fichier et de sortir une sortie simple : contig, début CDS, fin CDS, orientation.')
parser.add_argument('-psql', required=True, help=' fichier psql a traiter')
parser.add_argument('-out', required=False, help=' fichier de sortie : fichier texte : contig,début CDS, fin CDS, orientation')


args = parser.parse_args()
psql=args.psql
out=args.out

#------------------------------------------------------------------------------------------------------














#------------------------------------------------------------------------------------------------------

# C'est parti
tmp=open(out,"w")
nb_contig_tot=0
nb_contig_annote=0 


#Pour chaque ligne du fichier psql d'origine.
for line in open(psql):
	
	nb_contig_tot+=1
	line=line.strip()
	line=line.split(",")
	
	#Quel est le contig de la ligne?
	Contig=line[0]

	#Je reformate le nom pour retrouver le nom de mon fichier fichier fasta
	Contig=re.sub("_p$","",Contig)
	Contig=re.sub("_n$","",Contig)
	Contig=Contig.replace("_Con","|Con")
	Contig=Contig.replace("_like","|like")
	Contig=Contig.replace("_compl","|compl")
	Contig=Contig.replace("_orig","|orig")
	Contig=Contig.replace("_less","|less")

	#Petite partie : au cas ou j'ai des "," dans les noms, je dois recoller les morceaux a posteriori (c'est le cas pour certains noms d'EPO"
	if len(line)==12:
		line[0]=line[0]+","+line[1]
		line[1]=line[2]+","+line[3]
		line.pop(2) ; line.pop(2)
	
	if len(line)==14:
		line[0]=line[0]+","+line[1]+","+line[2]
		line[1]=line[3]+","+line[4]+","+line[5]
		line.pop(2) ; line.pop(2) ; line.pop(2) ; line.pop(2)

	#Quel est le début?
	if line[4]=="":
		start=line[5]
	else:
		start=line[4]
	
	#Quelestlafin?
	if line[8]=="" :
		end=line[7]
	else:
		end=line[8]
	
	#Quelestl'orientation?
	if int(line[6])<0:
		orientation="antisens"
	else:
		orientation="sens"
	
	
	#Impression du fichier de sortie
	tmp.write(Contig+"\t"+str(start)+"\t"+str(end)+"\t"+str(orientation)+"\n")
 

 
print "-----\n\nLe fichier psql a été parsé.\n\tNbr de contigs annotés : "+str(nb_contig_tot)+"\n\n----\n"

#------------------------------------------------------------------------------------------------------

