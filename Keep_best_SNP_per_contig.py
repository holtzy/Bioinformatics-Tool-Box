#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON 
#     				A partir des fichiers Info_SNP, je vais trouver et sortir le meilleur SNP de chaque contig.
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


parser = argparse.ArgumentParser(description= 'Permet de récupérer le meilleur SNP de chaque contig dans un fichier de SNP')
parser.add_argument('-SNP', required=True, help=' fichier de SNP ')
parser.add_argument('-out', required=True, help=' fichier de sortie ')

args = parser.parse_args()
SNP=args.SNP
out=args.out







#-------------------------------
#-------------------------------
		
# DICO des SNPs

#Je créé 2 dictionnaire : dans le dico_snp je met la ligne du SNP. dans le dico_contig, je met tous les SNPs existant pour un contig donné.
dico_snp=dict()
dico_contig=dict()
num=0

for line in open(SNP):
	num+=1
	line=line.strip()
	mem=line
	line=line.split()
	contig=line[6].replace(">","")

	#Je raccourci le nom si c'est un blé tendre, car sinon c'est inbuvable
	if contig.startswith("Traes"):
		contig=contig.split("|")[0]
	
	snp_name=contig+"@"+line[7]

	dico_snp[snp_name]=mem
	
	if contig not in dico_contig:
		dico_contig[contig]=snp_name
	else:
		dico_contig[contig]=dico_contig[contig]+",,"+snp_name		
		
print "\n\n------"
print "nbr de SNP total dans le fichier d'entrée : "+str(num)
print "------\n\n"

		
#-------------------------------
#-------------------------------










#-------------------------------
#-------------------------------

# SELECTION du meilleur SNP.

good_SNP=dict()
tmp=open(out,"w")

for contig in dico_contig :
	all_SNP=dico_contig[contig].split(",,")

	#Si j'ai qu'un SNP, je le prends
	if len(all_SNP) == 1 :
		good_SNP[contig] = dico_snp[all_SNP[0]]

			
	# Et si il y en a plusieurs? Je dois faire un choix en utilisant les paramètres : Couverture, FIS, He
	# Règle de décision ?
	# Si j'ai plusieurs SNPs avec FIS > 0.7, ou aucun, alors je prendrai celui qui a la meilleure couverture.
	# Si j'ai un seul SNPs avec FIS > 0.7 alors je le prend
	else :
			
			#Initialisation des valeurs pour le contig en question
			memoire=0 ; old_seuil_FIS="no" ; old_He=0 ; old_nb_indiv=0
			
			#Pour chaque SNP du contig
			for i in all_SNP:
				
				#Je récupère ses infos
				mem=dico_snp[i]
				line=dico_snp[i].split()
				He=line[3]
				nb_indiv=line[5]
				if line[4]=="-" :
					seuil_FIS="no"
				elif float(line[4]) > 0.7 :
					seuil_FIS="yes"
				else :
					seuil_FIS="no"					
					
				#Je dois comparer ces infos a celles du SNP précédent pour déterminer si ce SNP est meilleur ou pas
				
				# --- Si le SNP passe le seuil de FIS contrairement a tous les SNPs précédent, alors je le garde.
				if seuil_FIS=="yes" and old_seuil_FIS=="no" :
					good_SNP[contig] = mem
					old_seuil_FIS="yes"
					old_nb_indiv=nb_indiv
					old_He=He
					continue

				# --- Si le SNP a un FIS en dessous du seuil mais qu'il en va de même pour tous les autres SNP, alors on va choisir sur tous les autres critères.
				# --- Si il y a eu plusieurs SNPs au dessus du seuil, alors il faut choisir sur les autres critères.
				if ( seuil_FIS=="no" and old_seuil_FIS=="no" )  or  ( seuil_FIS=="yes" and old_seuil_FIS=="yes" ) :
							
					if nb_indiv > old_nb_indiv :
						good_SNP[contig] = mem
						old_nb_indiv=nb_indiv
						
					if nb_indiv==old_nb_indiv :
						if He > old_He:
							good_SNP[contig] = mem
							old_He=He
				
				
				
				

#Impression du fichier de sortie
for i in good_SNP:
	tmp.write(good_SNP[i]+"\n")
	
	

print "\n\n------"
print "Nbr de SNPs récupérés à raison de 1 par contig : "+str(len(good_SNP))	
print "------\n\n"



	
	

























