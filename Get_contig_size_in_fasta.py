#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON Get_Contigs_size.py
#     				Récupere la taille de tous les contigs d'un fichier fasta
#					Les fichiers d'entrées sont :
#						- Fichier fasta des contigs
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
parser.add_argument('-fasta', required=True, help='	séquences fasta des contigs d\'étude')
parser.add_argument('-out', required=True, help='	output')


args = parser.parse_args()
fic_fasta=args.fasta
out=args.out


#------------------------------------------------------------------------------------------------------

### STEP 0 : Collecte de l'ensemble des contigs dans un dictionnaire, j'en profite pour récupérer la longueur des contigs
dico_des_contigs=dict()
tot=0

for line in open(fic_fasta):
	line=line.strip()
	
	if line.startswith(">"):
		tot=tot+1
		contig_name=line.replace(">","")
		
		#Je raccourci le nom si c'est un blé tendre, car sinon c'est inbuvable
		if contig_name.startswith("Traes"):
			contig_name=contig_name.split("|")[0]

		# Et j'initialise la variable longueur	
		longueur_du_contig=0
		
	else:
		longueur_du_contig=longueur_du_contig + len(line)
		dico_des_contigs[contig_name]=str(longueur_du_contig)
		
#Fabrication de l'en tete
entete="contig_name"+"\t"+"longueur_en_pb"+"\t"

#------------------------------------------------------------------------------------------------------
	



#------------------------------------------------------------------------------------------------------

# STEP 1: Enregistrement
tmp=open(out,"w")
tmp.write(entete+"\n")
for contig_name in dico_des_contigs:
	tmp.write(contig_name+"\t"+dico_des_contigs[contig_name]+"\n")
tmp.close()


#------------------------------------------------------------------------------------------------------


print "\n-------\n\nTous les contigs de blé dur ont été répertorié, il y en a "+str(tot)+" \n\n------\n"




