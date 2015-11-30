#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON Take_some_positions_from_alr.py
#     				Prend un fichier ALR ou GEN, et récupère toutes les lignes spécifiées dans une liste
#  					Yan Holtz, yan1166@hotmail.com
#-------------------------------------

import os
import sys
import re
try:
  import argparse
except ImportError:
  print"oops, the import /argparse/ didn't work"


parser = argparse.ArgumentParser(description= 'Prend un fichier ALR ou GEN, et récupère toutes les lignes spécifiées dans une liste')
parser.add_argument('-input', required=True, help=' fichier .alr de résultat du mapping')
parser.add_argument('-positions', required=True, help=' Position à récupérer : format contig / tabulation / position')
parser.add_argument('-out', required=True, help=' fichier de sortie')


args = parser.parse_args()
input=args.input
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
 
print "\n\n\n---"
print str(nbr_de_position_au_total)+" positions sont présentes dans la liste donnée."+"\n"+"Ces positions concernent "+str(nbr_de_contigs_au_total)+" contigs différents au total."+"\n"
print "---\n\n\n"

#------------------------------------------------------------------------------------------------------


 
 








 
 
#------------------------------------------------------------------------------------------------------

### STEP 2  -- Je parse mon fichier ALR a la recherche de mes positions !!

tmp=open(output,"w")

#Je parcours mon fichier ALR
for line in open(input) :
	line=line.strip()
	
	#Lorsque je change de contig, je réinitialise tout
	if line.startswith(">") :
		liste=[]
		num=-1
		contig=line.replace(">","")
		if contig in pos_to_check:
			tmp.write(line+"\n")
			for i in pos_to_check[contig].split(","):
				liste.append(i)


	
	#Sinon je me balade a la recherche de ma position
	else:
		#Si je suis dans un contig d'intéret (sinon ca sert a rien)
		if contig in pos_to_check:
	
			num+=1
		
			#Si je suis dans une position ciblée:
			if str(num) in liste or num==0:
				tmp.write(line+"\n")		
		


print "\n\n\n---"
print("Fin du taff, Merci pour votre visite, votre fichier "+output+" est prêt, bonne lecture")
print "---\n\n\n"
