#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON Get_Cds_From_Pfas.py
#     				Permet d'obtenir les CDS et UTR d'un fichier .pfas (sortie du programme geno2pfas) avec l'aide d'un fichier d'annotation (.psql, sortie de prot4est)
#  					Yan Holtz, yan1166@hotmail.com
#-------------------------------------

import os
import sys
import re
try:
  import argparse
except ImportError:
  print"oops, the import /argparse/ didn't work"


parser = argparse.ArgumentParser(description= 'permet d\'utiliser les infos d\'un fichier pfas et d\'un psql pour soustirer les CDS et UTR')
parser.add_argument('-pfas', required=True, help=' fichier .pfas d\'entrée')
parser.add_argument('-annot', required=True, help=' fichier fasta du blé tendre contenant les positions des CDS dans l\'en tête.')


args = parser.parse_args()

fic_pfas=args.pfas
fic_annot=args.annot






#------------------------------------------------------------------------------------------------------
#STEP ONE  : DICTIONNAIRE DES ANNOTATIONS
#En fait je créé un dictionnaire avec comm clef : nom du contig puis 2 champs : longueur de l'UTR5 et longueur de l'UTR3
print "\n\n-----\n\nLe dictionnaire des annotations est en cours de réalisation ...\n\n---\n\n "

dico_size=dict()

for line in open(fic_annot) :
	line.strip()
	
	if line.startswith(">"):
		line=line.split("|")
		
		#### - Je vire les cas chelou :
		#ou je n ai pas d'utr ni de CDS --> l'annotation du contig n'a pas été possible
		#Ou j'ai des UTR décomposé.. Comment c'est possible ?
		if len(line) < 7 or len(line[6].split(";"))>1 or len(line)>10 and len(line[11].split(";"))>1:
			continue
			
		### - Je récupère les principales infos: 
		#--Nom du contig
		contig=line[0].replace(">","")
		
		#--Taille de l'UTR5 : si j'ai juste les pipe sans valeurs, c'est que pas d'UTR5, donc taille=0
		if line[7] == "" and line[6] == "" :
			size_UTR5=0
		else:
			size_UTR5=abs(int(line[7])-int(line[6]))+1
		
		#--Taille de l'UTR3 : attention, si je n'ai pas d'UTR3, alors je n'ai meme pas les pipes à la fin...
		if len(line)<11 :
			size_UTR3=0
		else:
			size_UTR3=abs(int(line[10])-int(line[11]))+1


		#j'enregistre dans le dico
		dico_size[contig]=str(size_UTR5)+","+str(size_UTR3)
#------------------------------------------------------------------------------------------------------











#------------------------------------------------------------------------------------------------------
#STEP 2 : lecture du PFAS et découpe des séquences
print "-----\n\nDécoupage du PFAS ...\n\n---\n\n "


#Je créé mon fichier de sortie
tmp_UTR5=open("UTR5.pfas","w")
tmp_CDS=open("CDS.pfas","w")
tmp_UTR3=open("UTR3.pfas","w")

debut="yes"


for line in open(fic_pfas) :
	line=line.strip()
	
	if line.startswith(">"):
	
	
		######### JE FINALISE LE CONTIG PRECEDENT
		if debut != "yes" :
			longueur=len(sequence)	
			
			#Si le contig a des données d'annotation, je vais l'ajouter a mon fichier de sortie !
			if contig in dico_size:
			
				#je récupère mes tailles d'UTR
				size_UTR5=int(dico_size[contig].split(",")[0])			
				size_UTR3=int(dico_size[contig].split(",")[1])			
			
				#Je coupe mes séquences et je les imprime
				sequence_UTR5=sequence[0:size_UTR5]
				tmp_UTR5.write(entete+"\n"+sequence_UTR5+"\n")
				sequence_CDS=sequence[ size_UTR5 : (longueur-size_UTR3) ]
				tmp_CDS.write(entete+"\n"+sequence_CDS+"\n")
				sequence_UTR3=sequence [ longueur-size_UTR3 : longueur ]
				tmp_UTR3.write(entete+"\n"+sequence_UTR3+"\n")
		
	



		######### JE PREPARE LE CONTIG SUIVANT

		entete=line
		line=line.split("|")
		sequence=""
		debut="no"
		
		### - Je récupère les principales infos: 
		#--Nom du contig
		contig=line[0].replace(">","")
		
	else :
		sequence=sequence+line



#Finalisation des fichiers de sortie		
tmp_UTR5.close()
tmp_CDS.close()
tmp_UTR3.close()


print "\n\n-----\n\nC'est dans la poche !! \n\n-------\n\n "




		
		