#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON Recode_SNP_file
#     				A Partir d'une matrice de SNP , transforme les génotypes de chaque indiv en format 0,1,2
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


My_description=" \n\n Transforme un fichier de génotypage en format A/A en format 0,1,2\n\n"

parser = argparse.ArgumentParser(description= My_description)
parser.add_argument('-SNP', required=True, help='	fichier de SNP à analyser. Format : une ligne d\'en tete, une premiere colonne donnant le Nom du SNP, puis chaque colonne donnant le génotype de l\'individu au format "A/A"\n\n')
parser.add_argument('-out', required=True, help='	output')


args = parser.parse_args()
fic_SNP=args.SNP
out=args.out








#------------------------------------------------------------------------------------------------------

### STEP 0 : Transformation
tmp=open(out , "w")

deb="yes"

#Je lis mon fichier de SNP
for line in open(fic_SNP) :
	line=line.strip()

	#Cas de l'entete que je remet en début de fichier
	if deb=="yes":
		tmp.write(line)
		deb="no"
		
	#C'est parti pour la lecture du fichier
	else:
		
		#J'initialise mes variables au début de chaque ligne. J'affiche le nom du SNP / marqueur dans le fichier de sortie
		nucl1=""
		line=line.split("\t")
		tmp.write("\n"+line[0])
		line=line[1:]
		
		#Puis je transforma pour chaque individu:
		for i in line:

			#Si je suis sur une donnée manquante, j'affiche un tiret :
			if i == "-":
				tmp.write("\t-")
			
			#Sinon je peux taffer :
			else:
			
				#récupération des allèles
				all1=i.split("/")[0]
				all2=i.split("/")[1]
			
				#Si je suis sur un hétérozygote
				if all1 != all2:
					tmp.write("\t1")
				
				#Si je suis sur un homozygote
				else:					
					#Si nucl1 n'existe pas, que je suis sur un homozygote
					if nucl1=="" :
						nucl1 = all1
						tmp.write("\t0")
						continue
	
					#Si je suis sur un homozygote équivalent a nucl1
					if all1==nucl1 :
						nucl1 = all1
						tmp.write("\t0")
	
					#Si je suis sur un homozygote différent de nucl1 et de nucl2
					if all1 != nucl1 and nucl1!= "":
						tmp.write("\t2")
			

#------------------------------------------------------------------------------------------------------






