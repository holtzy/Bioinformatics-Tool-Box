#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON 
#     				A partir des fichiers Info_SNP, je vais trouver et sortir le meilleur SNP de chaque contig et le baytes qui va avec
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
parser.add_argument('-SNP_EPO', required=True, help=' fichier Info_SNP ')
parser.add_argument('-baytes', required=True, help=' fichier contenant la totalité des baytes ')
parser.add_argument('-out', required=True, help=' fichier de sortie ')

args = parser.parse_args()
SNP_EPO=args.SNP_EPO
baytes=args.baytes
out=args.out







#-------------------------------
#-------------------------------


#----------1/ Je fais un dico : pour chaque snp présent dans le fichier de tous les baits :

dico_snp_in_baits=dict()

for line in open(baytes):
	line=line.strip()
	
	if line.startswith(">"):
		line=line.replace(">","")
		snp=line.split("|All")[0]
		dico_snp_in_baits[snp]=""

nbr_de_snp=len(dico_snp_in_baits)
print "nombre de SNPs différents présents dans le fichier de Baits : "+str(nbr_de_snp)

	
#-------------------------------
#-------------------------------






#-------------------------------
#-------------------------------
		
#----------2/ Je fais un dico des SNPs présents chez EPO
dico_snp_EPO=dict()
for line in open(SNP_EPO):
	line=line.strip()
	line=line.split()
	contig=line[6].replace(">","")
	He=line[3]
	nb_indiv=line[5]

	#Je raccourci le nom si c'est un blé tendre, car sinon c'est inbuvable
	if contig.startswith("Traes"):
		contig=contig.split("|")[0]
	
	snp_name=contig+"@"+line[7]

	dico_snp_EPO[snp_name]= str(He)+","+str(nb_indiv)

nbr_de_snp=len(dico_snp_EPO)
print "nombre de SNPs présents dans le fichier des SNPs EPOs : "+str(nbr_de_snp)


#-------------------------------
#-------------------------------
		
		
		
		



#-------------------------------
#-------------------------------
		
#----------3/ Je fais un dico : pour chaque contig je met tous les SNPs ayant un bayte et leurs caractéristiques:
	
dico_contig=dict()
for snp in dico_snp_in_baits :
	contig=snp.split("@")[0]	
	position=snp.split("@")[1]

	if snp in dico_snp_EPO:
		presence_EPO="yes"
		He_EPO=dico_snp_EPO[snp].split(",")[0]
		nb_indiv_EPO=dico_snp_EPO[snp].split(",")[1]
	else:
		presence_EPO="no"
		He_EPO="-"
		nb_indiv_EPO="-"

	
	if contig not in dico_contig:
		dico_contig[contig]=str(contig)+"@"+str(position)+"@"+str(He_EPO)+"@"+str(nb_indiv_EPO)
	else :		
		dico_contig[contig] = dico_contig[contig]+",,"+str(contig)+"@"+str(position)+"@"+str(He_EPO)+"@"+str(nb_indiv_EPO)


nbr_de_contig=len(dico_contig)
print "nombre de contigs différents présents dans le fichier des baits : "+str(nbr_de_contig)


#-------------------------------
#-------------------------------











#-------------------------------
#-------------------------------

#----------2/ Maintenant contig par contig je vais sélectionner le meilleur SNP.
good_SNP=dict()

for contig in dico_contig:
	all_SNP=dico_contig[contig].split(",,")

	#Si j'ai qu'un SNP, je le prends
	if len(all_SNP) == 1 :
		all_SNP=all_SNP[0].split("@")
		to_add=(all_SNP[0]+"@"+str(all_SNP[1]))
		good_SNP[to_add] = ""
	
	#Sinon je garde uniquement les SNPs qui ont équivalent EPO. Si aucun, je choisi au pif.
	else:
		
		#Je compte le nbr de SNP ayant un équivalent EPO
		num=0
		for i in all_SNP:
			presence_EPO=i.split("@")[2]
			if presence_EPO=="yes" :
				ayant_EPO = i
				num+=1

		#Si personne n'a d'équivalent EPO, alors je garde le premier qui vient au pif	
		if num==0:
			i=i.split("@")
			to_add=(i[0]+"@"+str(i[1]))
			good_SNP[to_add] = ""

		#Si il y en a un seul qui a un équivalent EPO
		if num==1:
			ayant_EPO=ayant_EPO.split("@")
			to_add=(ayant_EPO[0]+"@"+str(ayant_EPO[1]))
			good_SNP[to_add] = ""

			
		#Et si il y en a plusieurs? Je prend celui qui a le meilleure FIS, et au Pif si il y a des égalités.
		if num>1 :
			
			memoire=0
			for i in all_SNP:
				nb_indiv_EPO=i.split("@")[3]
				if nb_indiv_EPO != "-":
					if nb_indiv_EPO > memoire:
						to_keep=i
						memoire=i
							
			to_keep=to_keep.split("@")
			to_add=(to_keep[0]+"@"+str(to_keep[1]))
			good_SNP[to_add] = ""


print "nbr de paires  bayte a récupérer a raison de 1 par contig : "+str(len(good_SNP))	

#-------------------------------
#-------------------------------









#-------------------------------
#-------------------------------

#----------   RECUPERATION DES BAYTES

#Attention, je peux avoir un SNP dédoublé dans mon fichier de baits, il ne faudra pas le mettre 2 fois dans le fichier de sortie..
already_seen=dict()
tmp=open(out,"w")
nb=0

for line in open(baytes):
	if line.startswith(">"):
		contig=line.replace(">","")
		snp=contig.split("|All")[0]
		snp_plus_allele=line.split("|pos")[0]
		
		if snp in good_SNP and snp_plus_allele not in already_seen :
			toprint="yes"
			already_seen[snp_plus_allele]="yes"
			nb+=1
		else:
			toprint="no"
			
			
	if toprint=="yes":
		tmp.write(line)


print "nombre de séquence dans le fichier de sortie : "+str(nb)


#-------------------------------
#-------------------------------



	
	

























