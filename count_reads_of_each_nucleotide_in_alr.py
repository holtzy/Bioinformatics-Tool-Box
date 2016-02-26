#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON Count_reads_of_each_nucleotide_in_alr.py
#     				Script réalisé pour Muriel. Il part d'un fichier alr. Il va récupérer toutes les positions polymorphe et écrire, pour chaque individu, son nombre de A,C,G et T
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


parser = argparse.ArgumentParser(description= 'Ce script permet de sortir couverture moyenne, max et min pour des positions données, + nombre d\'individus avec plus de 0 puis plus de 10 reads d\'un fichier ALR')
parser.add_argument('-alr', required=True, help=' fichier .alr de résultat du mapping')
parser.add_argument('-out', required=True, help=' fichier de sortie')


args = parser.parse_args()
alr=args.alr
output=args.out








#------------------------------------------------------------------------------------------------------
### STEP 0 : réalisation de l'en tête du fichier de sortie

entete="Contig"+"\t"+"Position"
fic = open(alr,'r')
my_text = fic.readline()
my_text = fic.readline()
my_text=my_text.split()
for col in range(2,len(my_text)):
	entete=entete+"\t"+my_text[col]+"_A"+"\t"
	entete=entete+"\t"+my_text[col]+"_C"+"\t"
	entete=entete+"\t"+my_text[col]+"_G"+"\t"
	entete=entete+"\t"+my_text[col]+"_T"+"\t"
fic.close()

#------------------------------------------------------------------------------------------------------






#------------------------------------------------------------------------------------------------------

### STEP 1  -- Je parse mon fichier ALR a la recherche de mes positions !!

tmp=open(output,"w")
tmp.write(entete+"\n")

#Je parcours mon fichier ALR
for line in open(alr) :
	line=line.strip()
	
	#Lorsque je change de contig, je réinitialise tout
	if line.startswith(">") :
		#Alors j'ai un nouveau nom de contig
		contig=line.replace(">","")
		#Je réinitialise ma position
		position=-1
		
	#Ensuite je me balade dans mon contig
	else:
		line=line.split("\t")
		position=position+1
		
		# Je vais récupérer que les positions ou j'ai du polymorphisme !
		if line[1]=="P": 
		
			#j'initialise ma ligne de sortie avec le nom du contig et la position :
			my_line=contig+"\t"+str(position)
	
			#Je récupère les infos des individus un par un : nbr de A, de C, de G et de T.
			for col in range(2,len(line)):
				text=line[col]
				text=text.split("[")[1]
				text=re.sub(r"]",r"",text)
				text=re.sub(r"/",r"\t",text)
				my_line=my_line+"\t"+text
			
			#J'imprime la ligne dans l'output
			tmp.write(my_line+"\n")

			


print "\n\n\n---"
print("Fin du taff, Merci pour votre visite, votre fichier "+output+" est prêt, bonne lecture")
print "---\n\n\n"


#------------------------------------------------------------------------------------------------------

 
 

