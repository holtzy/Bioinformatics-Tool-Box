#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON Nucleotidic_distance_between individuals
#     				A Partir d'une matrice de SNP, calcule la distance entre accessions.
#					Les fichiers d'entrées sont :
#						- matrice de SNP 
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


My_description="\n\nA partir d\'une matrice de SNP, ce programme va calculer la distance entre individus."+"\n---\n"+"La distance est calculée tel que :  distance = (nb de génotype identiques / nb de cas identiques) + nb de cas différents"+"\n---\n"+"Format d'entrée : une ligne d'entete qui comporte le nom des individus de chaque colonne. Et une seule premiere colonne qui donne le nom du SNP ou marqueur."+"\n---\n"+"La matrice de sortie est une demi matrice supérieure\n\n"

parser = argparse.ArgumentParser(description= My_description)
parser.add_argument('-SNP', required=True, help='	fichier de SNP à analyser. Format : une ligne d\'entete avec nom_SNP + tab + nom_individu1 + tab + nom2 etc.... Puis une ligne par SNP, avec le nom du SNP puis le genotype de chaqueindividu en colonne, au format A/T par exemple et des tiret pour les NA')
parser.add_argument('-out', required=True, help='	output')


args = parser.parse_args()
fic_SNP=args.SNP
out=args.out








#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 0 : Détermination du noombre d'individus

fic = open(fic_SNP,'r')
entete = fic.readline();
nb_indiv=int(len(entete.split("\t")))
fic.close()

print"\n-------\n\nLe fichier de SNP comporte "+str(nb_indiv-1)+" accessions\n-------\n\n"









#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 1 : Calcul des distances
tmp=open(out , "w")


#Fabrication de m'en tete du fichier de sortie. Je transforme l'entete initiale en une nouvelle "to_add"
entete=entete.strip()
entete=entete.split("\t")
to_add="versus"
for i in entete[1:]:
	to_add=to_add+"\t"+i
tmp.write(to_add)




#Je commence les calculs ligne par ligne = accession par accession
for i in range(1,nb_indiv):

	print "traitement de l'individu "+str(i)+" sur "+str(nb_indiv-1)+" individus"

	indiv=entete[i]
	tmp.write("\n" + indiv)
	
	#Pour chaque autre accession
	for j in range(1 , nb_indiv) :
		
		#Si on a déja traité l'accession, alors on met juste un tiret
		if j<i:
			tmp.write("\t-")
		
		#Sinon on fait le calcul de distance
		if j>=i :
			nb_same=0
			nb_different=0
			num_ligne=0
			for line in open(fic_SNP):
				line=line.strip()
				num_ligne+=1
				if num_ligne>1:
					line=line.split("\t")
					a=line[i]
					b=line[j]
					if a!="-" and b!="-":
						if a==b:
							nb_same+=1
						else:
							nb_different+=1
			distance=float(nb_same)/float(nb_same+nb_different)
			tmp.write("\t"+str(distance))

tmp.write("\n")





















