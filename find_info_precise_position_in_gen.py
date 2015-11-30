#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON Find_Info_precise_position_in_gen.py
#     				Ce script permet de sortir des stats pour des positions données d'un fichier GEN :
#						- Nombre d'homozygotes
#						- Nombre d'hétérozygotes et pourcentage
#						- manquant
#  					Yan Holtz, yan1166@hotmail.com
#-------------------------------------


import os
import sys
import re
try:
  import argparse
except ImportError:
  print"oops, the import /argparse/ didn't work"


parser = argparse.ArgumentParser(description= 'Ce script permet de sortir nb homo, nb hétéro, nb manquant d\'un fichier GEN')
parser.add_argument('-gen', required=True, help=' fichier .gen de résultat de reads2snp')
parser.add_argument('-positions', required=True, help=' Position à étudier : format contig / tabulation / position')
parser.add_argument('-out', required=True, help=' fichier de sortie')


args = parser.parse_args()
gen=args.gen
output=args.out
positions=args.positions








#------------------------------------------------------------------------------------------------------

### STEP -1 : Dictionnaire des positions à étudier
nbr_de_position_au_total=0
nbr_de_contigs_au_total=0
pos_to_check=dict()

for line in open(positions):
	nbr_de_position_au_total+=1
	line=line.split()
 	contig=line[0].replace(">","")
 	
 	if contig not in pos_to_check:
 		nbr_de_contigs_au_total+=1
 		pos_to_check[contig]=line[1]
 		
 	else:
 		pos_to_check[contig]=pos_to_check[contig]+","+line[1]
 
pos_redondante=nbr_de_position_au_total - len(pos_to_check)
print "\n\n\n---"
print str(nbr_de_position_au_total)+" positions sont présentes dans la liste donnée."+"\n"+"Ces positions concernent "+str(nbr_de_contigs_au_total)+" contigs différents au total."+"\n"+"Cependant, "+str(pos_redondante)+" positions données sont redondantes"
print "---\n\n\n"

#------------------------------------------------------------------------------------------------------


 
 








 
 
#------------------------------------------------------------------------------------------------------

### STEP 2  -- Je parse mon fichier ALR a la recherche de mes positions !!

tmp=open(output,"w")
tmp.write("Contig"+"\t"+"Position"+"\t"+"Nb_homo"+"\t"+"Nb_hetero"+"\t"+"Pourc_hetero"+"\t"+"Nb_missing"+"\n")


#Je parcours mon fichier GEN
for line in open(gen) :
	line=line.strip()
	
	#Lorsque je change de contig, je réinitialise tout
	if line.startswith(">") :
		liste=[]
		num=-1
		contig=line.replace(">","")
		if contig in pos_to_check:
			for i in pos_to_check[contig].split(","):
				liste.append(i)

				
	#Sinon je me balade a la recherche de ma position
	else:
		line=line.split("\t")
		num+=1
		
		#Si je suis dans une position ciblée:
		if str(num) in liste:
			vecteur=[]
			nb_homo=0
			nb_hetero=0
			nb_missing=0
			
			#Je récupère les infos des individus un par un
			for col in range(1,len(line)):
				a=line[col].split("|")[0]
				vecteur.append(a)
			
			#Je calcule des stats sur le vecteur obtenu
			for i in vecteur:
				if i=="NN":
					nb_missing+=1

				else:
					if i[0]==i[1]:
						nb_homo+=1
					else:
						nb_hetero+=1
			
			tot=nb_homo+nb_hetero
			if tot !=0:
				pourc_hetero=float(nb_hetero)/float(nb_homo+nb_hetero)*100
			else:
				pourc_hetero="-"
				
			#J'inscris quelques stats dans le fichier de sortie
			tmp.write(contig+"\t"+str(num)+"\t"+str(nb_homo)+"\t"+str(nb_hetero)+"\t"+str(pourc_hetero)+"\t"+str(nb_missing)+"\n")
		
		


print "\n\n\n---"
print("Fin du taff, Merci pour votre visite, votre fichier "+output+" est prêt, bonne lecture")
print "---\n\n\n"

	