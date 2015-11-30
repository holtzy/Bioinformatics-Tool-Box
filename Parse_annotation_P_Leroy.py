#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT Parse_annotation_P_Leroy.py
#
#     				Philippe leroy m'a envoyé le résultat du blast de EPO sur blé tendre avec les position d'exon correspondantes. le pb, c'est qu'il faut trouver le best blast + faire des changements de formats, pour obtenir mon format position exons classique.
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


parser = argparse.ArgumentParser(description= 'Philippe leroy m\'a envoyé le résultat du blast de EPO sur blé tendre avec les position d\'exon correspondantes. le pb, c\'est qu\'il faut trouver le best blast + faire des changements de formats, pour obtenir mon format position exons classique.')
parser.add_argument('-input', required=True, help=' fichier de philippe leroy Position_exons')
parser.add_argument('-out', required=False, help=' fichier de sortie : fichier texte : contig,position des exons')


args = parser.parse_args()
input=args.input
out=args.out

#------------------------------------------------------------------------------------------------------




#------------------------------------------------------------------------------------------------------
print "-----\n\nLecture du fichier\n\n----\n"

dico_intron=dict()
dico_blast=dict()
nombre=0
nbr_uniq=0


#A chaque nouvelle ligne je récupère toutes les infos : contig, scores de blast, et position des introns
for line in open(input):
	nombre=nombre+1
	line=line.strip()
	line=line.split()
	
	#Je récupère le nom du contig de manière conventionnelle.
	contig=line[1]
	contig=contig.replace("_Con","|Con")
	contig=contig.replace("_like","|like")
	contig=contig.replace("_compl","|compl")
	contig=contig.replace("_orig","|orig")
	contig=contig.replace("_less","|less")
	
	#Je récupère les infos de la ligne pour ce contig
	identity=line[2]
	coverage=line[3]
	positions=line[4]
	

	#Si on recontre ce contig pour la premiere fois, on récupère la position de ses introns.
	if contig not in dico_intron:
		nbr_uniq+=1
		dico_intron[contig]=positions
		
		#Je stocke les valeurs de blast dans un dico à part (ca me servira si ce contig se représente)
		dico_blast[contig]=identity+","+coverage
	
	
	#Si le contig a déja été repertorié, il faut alors choisir son MEILLEUR blast et les positions d'introns correspondantes.
	else:
	
		#Si le nouveau contig a un meilleur % d'identité, alors on le garde d'office, il remplace alors l'ancien dans les dicos
		if identity > dico_blast[contig].split(",")[0]:
			dico_intron[contig]=positions
			dico_blast[contig]=identity+","+coverage
						
		#Si c est pile égal, alors on regarde la distance de couverture.
		if identity == dico_blast[contig].split(",")[0]:
			if coverage > dico_blast[contig].split(",")[1]:
				dico_intron[contig]=positions
				dico_blast[contig]=identity+","+coverage
								



print "-----\n\nDone ! Nombre de ligne dans le fichier d'entrée : "+str(nombre)+"\nNombre de contig concernés au total : "+str(nbr_uniq)+"\n\n----\n"

#------------------------------------------------------------------------------------------------------






#------------------------------------------------------------------------------------------------------
#Impression du fichier de sortie

tmp=open(out,"w")

for contig in dico_intron : 
	tmp.write(contig+"\t"+dico_intron[contig]+"\n")
	
#------------------------------------------------------------------------------------------------------

