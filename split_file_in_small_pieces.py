#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON RONDELE
#     				permet de couper un fichier alr en rondele rapidement
#  					Yan Holtz, yan1166@hotmail.com
#-------------------------------------

import os
import sys
import re
try:
  import argparse
except ImportError:
  print"oops, the import /argparse/ didn't work"


parser = argparse.ArgumentParser(description= ' \n\nPermet de couper un ALR ou un FASTA en plusieurs parties rapidement.    ----     Les fichiers obtenus seront dans un répertoire EXTRAIT\n\n')
parser.add_argument('-alr', required=True, help=' fichier .alr ou .fasta d\'entrée')
parser.add_argument('-nb', required=True, help=' nbr de contigs que l\'on veut par extrait')


args = parser.parse_args()

fic_alr=args.alr
nb=args.nb





#------------------------------------------------------------------------------------------------------

#Je crée un dossier EXTRAIT pour recevoir mes extraits
commande="mkdir EXTRAIT"
os.system(commande)


#Et je fais mes extraits
compte=0
num_extrait=1
tmp=open("EXTRAIT/extrait_1","w")

for line in open(fic_alr):
		line=line.strip()
		#Je décide dans quel fichier de sortie je place ma ligne
		if line.startswith('>'):
			compte=compte+1
			if int(compte) == int(nb):
				num_extrait=num_extrait+1
				compte=0
				nom_de_fichier="extrait_"+str(num_extrait)
				tmp=open(nom_de_fichier,"w")
		
		#Et je crée mes fichiers de sortie
		tmp.write(line+"\n")
		
tmp.close
		



#Je bouge tous mes extraits dans le dossier extrait
commande="mv extrait* EXTRAIT"
os.system(commande)

#------------------------------------------------------------------------------------------------------



		
							
