#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON SPLIT_FASTA_BY_CONTIGS
#     				Prend un fichier fasta, le sépare en plusieurs fichiers fasta contenant chacun un contig
#  					Yan Holtz, yan1166@hotmail.com
#-------------------------------------

import os
import sys
import re
try:
  import argparse
except ImportError:
  print"oops, the import /argparse/ didn't work"


parser = argparse.ArgumentParser(description= '\n\nPrend un fichier fasta, le sépare en plusieurs fichiers fasta contenant chacun un contig. Les Fastas Individuels seront placés dans un dossier CONTIGS_OBTENUS\n\n')
parser.add_argument('-fasta', required=True, help=' fichier fasta ')

args = parser.parse_args()

fic=args.fasta







#------------------------------------------------------------------------------------------------------
print "--------------\n\nLe decoupage de votre fichier Fasta a debuté !!!\n\n--------------\n"


fichier=""
useless=""

for line in open(fic):
	line=line.strip()
	

	if line.startswith('>'):						# Si la ligne est un nom de contig
		contig=line.split("|")[0]					# Je récupère le nom de contig en coupant avec le premier pipe est en enlevant le chevron
		contig=contig.replace(">","")
		
	if contig == fichier:							#Si le nom de contig est le même que le fichier en cours d'écriture
		tmp.write(line+"\n")						#Alors j ajour la ligne au fichier
	else:
		useless=useless+"-"
		print useless
		fichier = contig							#Sinon je commence un nouveau fichier ayant le nom du contig
		tmp=open(contig,"w")			
		tmp.write(line+"\n")						#J'oublie pas d'y ajouter la ligne quand même !
		
		
commande="mkdir CONTIGS_OBTENUS"
os.system(commande)
commande="mv Con* CONTIGS_OBTENUS"
os.system(commande)
	
		
print "\n--------------\n\nC'est gagné, votre fichier a été découpé avec succès !!!\n\n--------------"
	