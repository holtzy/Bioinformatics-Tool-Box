# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON 
#     				A partir d'un fichier de SNP format A/A C/C avec les parents en colonne 3 et 4, recode tous les génotype en A (allele Parent1) ou B(allele parent2)
#  					Yan Holtz, yan1166@hotmail.com
#-------------------------------------


import os
import sys
import re
try:
  import argparse
except ImportError:
  print"oops, the import /argparse/ didn't work"


parser = argparse.ArgumentParser(description= 'A partir d un fichier de SNP format A/A C/C avec les parents en colonne 3 et 4, fournit un fichier SNP/allele parent1/ allele parent2')
parser.add_argument('-SNP', required=True, help=' fichier de SNP ')
parser.add_argument('-out', required=True, help=' fichier de sortie ')

args = parser.parse_args()
SNP=args.SNP
out=args.out


#-------------------------------












#-------------------------------


tmp=open(out,"w")
nb_SNP=0
nb_perdu=0
nb_final=0

#Pour chaque ligne de mon fichier de SNP
for line in open(SNP):
	nb_SNP+=1
	line=line.strip()
	line=line.split()

	#Nom du contig :
	contig=line[0]
	if contig.startswith("Traes"):
		contig=contig.split("|")[0]


	#récupération allèle parent 1 et 2
	parent1=line[2]
	parent2=line[3]
	
	#Récupération des allèles des 2 parents : 
	#parent1
	if parent1=="-" :
		a1="-"
		a2="-"
	else :
		a1=parent1[0]
		a2=parent1[2]
	#parent2
	if parent2=="-" :
		b1="-"
		b2="-"
	else :
		b1=parent2[0]
		b2=parent2[2]
	
	
	#Dans certains cas, je ne garde pas le SNP : Si les parents sont égaux / si les 2 parents sont manquants / Si les 2 parents sont hétéro / Si un parent est hétéro et l'autre manquant
	if parent1 == parent2 or (a1!=a2 and b1!=b2) or (parent1=="-" and b1!=b2) or (parent2=="-" and a1!=a2 )  :
		nb_perdu+=1
		continue
	
	# Si les 2 parents sont OK, alors je vais imprimer ma sortie:
	if parent1 != parent2 and (a1==a2 and b1==b2) and parent1!="-" and parent2!="-" :
		nb_final+=1
		tmp.write( contig+"@"+line[1]+" "+a1+a2+ " " + b1+b2 + "\n")
		continue
	
	# Sinon il faut que je récupère les 2 allèles présents dans le SNP pour pouvoir deviner par différence.
	# Donc pour chaque individu hors Dic2 et Silur de ce SNP, je vais récup le génotype :
	for i in range(5,len(line)+1) :
	
		# Je récupère le génotype de l'indiv en question 
		genot=line[i-1]
			
		# si le génotype est manquant ou hétéro, il me donne pas d'info, je continue
		if genot=="-" :
			continue
			
		# Le génotype n'est pas manquant, je peux récupérer ses 2 alleles
		g1=genot[0]
		g2=genot[2]
		
		# Si ce génotype est hétérozygote, il me donne pas d'info, je continue
		if g1!=g2 :
			continue
			
		# CAS 1: génotype du parent 1 est manquant. Soit j'ai trouvé l'allele manquant et j'actualise, soit je ne l'ai pas trouvé et je continue
		if parent1=="-" or a1!=a2 :
			if genot==parent2:
				continue
			if genot!=parent2:
				tmp.write( contig+"@"+line[1]+" "+g1+g2+ " " + b1+b2 + "\n")
				nb_final+=1
				break

		# CAS 2: génotype du parent 2 est manquant. Soit j'ai trouvé l'allele manquant et j'actualise, soit je ne l'ai pas trouvé et je continue
		if parent2=="-" or b1!=b2 :
			if genot==parent1:
				continue
			if genot!=parent1:
				tmp.write( contig+"@"+line[1]+" "+a1+a2+ " " + g1+g2 + "\n")
				nb_final+=1
				break
		
	

	
print "\n\n--------"
print "nbr de SNP dans le fichier d'entrée : "+str(nb_SNP)
print "nbr de SNP perdu car parents manquants / parent identiques / 2 parents hétéro : "+str(nb_perdu)
print "Remaining SNP : "+str(nb_final)
print "--------\n\n"
