#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON Find_Blas_family
#     				Consider a blast result.
#					The contig 1 blast with the contig 4. And the contig 4 with the contig 18
#					Then you have a family which groups together the contigs 1,4 and 18
#					This script permits to find the family
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


parser = argparse.ArgumentParser(description=  '\n\nConsider a blast result.The contig 1 blast with the contig 4. And the contig 4 with the contig 18. Then you have a family which groups together the contigs 1,4 and 18. This script permits to find the family\n\n')
parser.add_argument('-blast', required=True, help='resultat de blastn. Basivally, column 1 and 2 are sufficient. This input must be filtered before with your threshold')
parser.add_argument('-out', required=True, help=' output file')


args = parser.parse_args()
fic_blast=args.blast
out=args.out







#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
num=0
nligne=0
contig_trouves=[]
liste_des_groupes=[]
numero_apartenance=dict()

#Je parse mon fichier
for line in open(fic_blast):
	nligne+=1
	line=line.strip()
	line=line.split()
	contig1=line[0]
	contig2=line[1]
	
	if contig1==contig2:
		continue

	# Si les 2 contigs ont déja été traités, alors rien d'autre à faire
	if contig1 in contig_trouves and contig2 in contig_trouves:
		continue
	
	# Si aucun des 2 n'a été trouvés, je vais devoir faire un nouveau groupe
	if contig1 not in contig_trouves and contig2 not in contig_trouves:
		
		#Ma nouvelle valeur de num est mon numéro de groupe
		num+=1
		#Je note que mes 2 contigs appartiennent a mon groupe numéro num
		numero_apartenance[contig1]=num
		numero_apartenance[contig2]=num
		
		#Je crée une nouvelle liste avec seulement ces 2 contigs et je l'ajoute à ma liste des groupes trouvées
		groupe=[]
		groupe.append(contig1)
		groupe.append(contig2)
		liste_des_groupes.append(groupe)
	
		#Je n'oublie pas de noté que ces 2 contigs ont déja été trouvés
		contig_trouves.append(contig1)
		contig_trouves.append(contig2)


	#Si seul le contig 1 a été trouvé, alors je dois ajouter le contig 2 au même groupe que le 1
	if contig1 in contig_trouves and contig2 not in contig_trouves:
		
		#Je retrouve le numéro de groupe du contig 1
		val=numero_apartenance[contig1]
		
		#Et j'ajoute le contig 2 à ce groupe
		liste_des_groupes[val-1].append(contig2)
		
		#Je n'oublie pas d'indiquer que le contig2 a été trouvé, et je lui donne son numéro de groupe
		contig_trouves.append(contig2)
		numero_apartenance[contig2]=val
	

	#Réciproque Si seul le contig 2 a été trouvé
	if contig2 in contig_trouves and contig1 not in contig_trouves:
		
		#Je retrouve le numéro de groupe du contig 1
		val=numero_apartenance[contig2]
		
		#Et j'ajoute le contig 2 à ce groupe
		liste_des_groupes[val-1].append(contig1)

		#Je n'oublie pas d'indiquer que le contig2 a été trouvé, et je lui donne son numéro de groupe
		contig_trouves.append(contig1)
		numero_apartenance[contig1]=val



print("\n\nNombre de ligne dans le resultat de blast = "+str(nligne))
print("Nombre de famille trouvées = "+str(num)+"\n\n")





#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Impression de la liste des groupes sous une sorte de format fasta

tmp=open(out , "w" )
num=0

for i in range(0,len(liste_des_groupes)):
	num+=1
	tmp.write(">groupe_"+str(num)+"\n")
	for v in liste_des_groupes[i] : 
		tmp.write(v+"\n")
		

		
	
	
	
	
	
	
	
	
	
	
	
	
	
	
