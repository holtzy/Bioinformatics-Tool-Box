#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON CUT_PFAS_IN_CONTIGS
#     				Take a PFAS file : split it contig by contig
#					2 output format are possible
#  					Yan Holtz, yan1166@hotmail.com
#-------------------------------------

import os
import sys
import re
try:
  import argparse
except ImportError:
  print"oops, the import /argparse/ didn't work"


parser = argparse.ArgumentParser(description= '\n\nTake a PFAS file : split it contig by contig\n\n')
parser.add_argument('-pfas', required=True, help=' input PFAS file')
parser.add_argument('-format', required=True, help=' output format : "nabholz" or "initial". With Nabholz, sequence names are renamed, permitting to use the diversity indexes calculation programmes of Benoit Nabholz')


args = parser.parse_args()

fic_pfas=args.pfas
format=args.format



#------------------------------------------------------------------------------------------------------

print "\n\n--------\n\nBienvenu ! Merci de me faire taffer, je commencais à m'ennuyer !\n\n------------\n\n"		



#Initiation
nb_contigs_traites=0
my_dico=dict()




print "\n\n--------\n\nRemplissage\n\n------------\n\n"		


for line in open(fic_pfas):
		line=line.strip()

		#Si la ligne est un nom de séquence
		if line.startswith('>'):
		
			#De quel contig s'agit il ?
			contig=line.split("|ind")[0]
			contig=contig.replace(">","")
			
			#Pour le blé tendre qui a des noms de contig a ralonge, n les raccourci
			if contig.startswith("Traes"):
				contig=contig.split("|")[0]
			
			#numéro de l'allèle ?
			allele=line.split("Allele_")[1]
			
		
			#On change la ligne pour le format Nabholz:			
			if format =="nabholz":
				line=re.sub(r">.*ind",r">EPOR",line)
				line=re.sub(r"\|Allele_",r".",line)
				
				
		#Et maintenant je rempli le dico
		if contig in my_dico:
			my_dico[contig]+=line+"\n"		
		else:
			my_dico[contig]=line+"\n"

#------------------------------------------------------------------------------------------------------







#------------------------------------------------------------------------------------------------------
#Impression du dico.
for contig in my_dico:
	nb_contigs_traites += 1
	fichier_de_sortie=contig+".fasta"	
	tmp=open(fichier_de_sortie , "w")	
	tmp.write(my_dico[contig])



				
print "\n\n--------\n\nTravail terminé, "+str(nb_contigs_traites)+" Contigs ont été extraits  !\n\n--------"	
#------------------------------------------------------------------------------------------------------
		
				
				
				
				
				
				
				
				
				
				
				
				
