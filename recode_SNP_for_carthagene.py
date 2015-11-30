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


parser = argparse.ArgumentParser(description= 'A partir d un fichier de SNP format A/A C/C avec les parents en colonne 3 et 4, recode tous les génotype en A (allele Parent1) ou B(allele parent2')
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
	
	
	#Dans certains cas, je ne garde pas le SNP : Si les parents sont égaux / si les 2 parents sont manquants / Si les 2 parents sont hétéro
	if parent1 == parent2 or (a1!=a2 and b1!=b2) :
		nb_perdu+=1
		continue
		
	#Sinon je vais garder le SNP donc je créé l'objet ligne finale que je vais faire grossir petit a petit
	toprint=contig+"@"+line[1]
	
	#Ensuite pour chaque individu hors Dic2 et Silur de ce SNP, je vais trouver le génotype :
	for i in range(5,len(line)+1) :
	
		#Je récupère le génotype de l'indiv en question 
		genot=line[i-1]
		
		#jinitialise Le génotype final de cet individu par un "-"
		a="-"
		
		#si le génotype est manquant, alors je laisse un tiret	que j'ajoute a toprint et je continue a l'indiv suivant
		if genot=="-" :
			toprint=toprint+"\t"+a 
			continue
			
		#Le génotype n'est pas manquant, je peux récupérer ses 2 alleles
		g1=genot[0]
		g2=genot[2]
		
		#Si ce génotype est hétérozygote, alors je lui laisse un "-"
		if g1!=g2 :
			toprint=toprint+"\t"+a
			continue
		
		#Si les 2 parents sont homozygotes différents, alors c'est simple, je regarde a qui ca correspond
		if a1==a2 and b1==b2 :
			if genot==parent1:
				a="A"
			if genot==parent2 :
				a="B"
						
		#Si parent1 est manquant ou hétérozygote :
		if parent1=="-" or a1!=a2 :
			if genot==parent2 :
				a="B"
			else:
				a="A"
						
						
		#Si parent2 est manquant ou hétérozygote :
		if parent2=="-" or b1!=b2 :
			if genot==parent1 :
				a="A"
			else:
				a="B"
		
		#J'ajoute l'allele de l'individu a la ligne
		toprint=toprint+"\t"+a
	
	#J'joute la ligne au fihier de sortie
	tmp.write(toprint+"\n")
	nb_final+=1

	
print "\n\n--------"
print "nbr de SNP dans le fichier d'entrée : "+str(nb_SNP)
print "nbr de SNP perdu car parents manquants / parent identiques / 2 parents hétéro : "+str(nb_perdu)
print "Remaining SNP : "+str(nb_final)
print "--------\n\n"
