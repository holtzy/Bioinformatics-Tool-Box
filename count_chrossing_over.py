#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON Count_chrossing_over.py
#     				Calculate the number of chrossing over for each individual of a genotyping matrix
#
#  					Yan Holtz, yan1166@hotmail.com
#-------------------------------------

# Format du fichier d'entrée : marker / chromo ou LG / position / une colonne par indiv avec A ou B.
	#marker;chromo_map;position_map;V2;V3;V4;V5;V6;V7;V8
	#Cluster_3354|Contig2|original@283;1A;0;B;A;A;B;A;-;B
	#cluster_singlet|EPO_092_2826919364122958018+,...,87210-|original@654;1A;0;-;A;A;B;A;-;B
	#Cluster_10322|Contig1|complementarySeq@393;1A;0.7;B;A;A;B;A;A;B
	#Cluster_3335|Contig1|complementarySeq@1206;1A;0.7;B;-;A;B;A;A;B
	#Cluster_3335|Contig1|likelySeq@1206;1A;0.7;B;-;A;B;-;-;B
	#Cluster_5949|Contig4|original@2220;1A;0.7;B;A;A;B;A;A;B
	#Cluster_7725|Contig5|original@540;1A;0.7;B;A;A;B;A;-;B
	#Cluster_1386|Contig2|original@880;1A;1.3;B;A;A;B;-;A;B
	#Cluster_6015|Contig5|original@225;1A;7.6;-;A;A;B;A;A;B

import os
import sys
import re
try:
  import argparse
except ImportError:
  print"oops, the import /argparse/ didn't work"


parser = argparse.ArgumentParser(description= 'Calculate the number of chrossing over for each individual of a genotyping matrix')
parser.add_argument('-genot', required=True, help='Fichier au format .raw mais ordonné dans le sens de la carte. Un example de fichier d entrée est donné dans le script')
parser.add_argument('-out', required=True, help=' fichier de sortie')


args = parser.parse_args()
fic_in=args.genot
out=args.out




#------------------------------------------------------------------------------------------------------
#Combien y a t'il d'individus dans le fichier + Récupération des noms d'individu
#------------------------------------------------------------------------------------------------------

fic = open(fic_in,'r')
header = fic.readline()
header=header.strip()
header=header.split(";")
size_header=len(header)




#------------------------------------------------------------------------------------------------------
# Récupération des chromosome
#------------------------------------------------------------------------------------------------------
num=0
chromo_list=[]
for line in open(fic_in):
	num+=1
	if num>1:
		line=line.split(";")
		chromo=line[1]
		if chromo not in chromo_list:
				chromo_list.append(chromo)
		
my_header="genotype"
for i in chromo_list:
	my_header=my_header+"\t"+i
	
	






#------------------------------------------------------------------------------------------------------
#Fonction qui lit le fichier pour un individu donné, et inscrit dans l'output le nbr de chrossing over :
#------------------------------------------------------------------------------------------------------

def count_chrossing(indiv):

	#Initialisation des variables
	nb_cross=0
	allele_mem="-"
	chromo_mem="-"
	to_print=header[indiv]
	num=0
	double_cross=0
	
	for line in open(fic_in):
		
		#Il faut passer le header..
		num=num+1
		if num>1:
			
			#Je récupere le chromosome et l'allele du marqueur
			line=line.strip()
			line=line.split(";")
			chromo=line[1]
			allele=line[indiv]
			
			#Si c'est un allele manquant, je passe a la ligne d'après !
			if allele=="-":
				continue
			
			#Si j'ai un changement d'allele, alors je marque un point a nb cross. et double cross vaut 1. Sinon, je réinitiailise double cross
			if allele != allele_mem and allele_mem!="-":
				nb_cross+=1
				double_cross+=1
			else:
				double_cross=0

			#Si j'ai un double crossing over, la variable double cross vaut 2. Dans ce cas la je retranche 2 points à nb_cross, car je ne veux pas compter les doubles cross !
			if double_cross==2:
				nb_cross=nb_cross-2
				double_cross=0
			
			# j'actualise mon allele memoire
			allele_mem=allele
		
			#Si je passe a un autre chromo, j'inscrit le nbr de crossing over dans toprint, et je remet le nbr de chross a 0
			if chromo_mem != chromo and chromo_mem !="-":
				to_print=to_print+"\t"+str(nb_cross)
				nb_cross=0
				
			#J'actualise mon chromosome mémoire
			chromo_mem=chromo
		
	# J'ajoute le dernier chromosome !
	to_print=to_print+"\t"+str(nb_cross)

	# J'imprime dans ma sortie ces données
	tmp.write(to_print+"\n")





#------------------------------------------------------------------------------------------------------
# Passage de la fonction
#------------------------------------------------------------------------------------------------------

tmp=open(out,"w")
tmp.write(my_header+"\n")

for indiv in range(3,size_header):
	count_chrossing(indiv)




print("work is done" )













