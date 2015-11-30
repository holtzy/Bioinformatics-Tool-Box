#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON Extract_CDS_from_PFAS.py
#     				Ce script permet de récupérer les CDS d'un fichier PFAS. Entrée : 
#						-un fichier PFAS : Pour chaque conig, une séquence par individu et par allèle.
#						-un fichier d'annotation, de la forme contig,debut CDS, fin CDS.
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


parser = argparse.ArgumentParser(description= 'get the baits from SNPs, fasta and annotation. ')
parser.add_argument('-pfas', required=True, help='	fichier pfas cible')
parser.add_argument('-position_CDS', required=True, help=' fichier indiquant la position des CDS dans les contigs. Format : contig,debut CDS, fin CDS')
parser.add_argument('-out', required=True, help=' fichier de sortie')


args = parser.parse_args()
pfas=args.pfas
position_CDS=args.position_CDS
output=args.out


#------------------------------------------------------------------------------------------------------











#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 1 : Collecte des positions des CDS
print "\n\nSTEP 1 : dictionnaire des positions des  CDS"


dico_CDS=dict()
nombre=0 


#A chaque nouvelle ligne je récupère toutes les infos : contig, positions des exons
for line in open(position_CDS):
	line.strip()
	line=line.split()
	nombre+=1
	
	#Je récupère les infos : il s'agit des début et fin des CDS
	contig=line[0]
	deb=int(line[1])
	fin=int(line[2])
	orientation=line[3]
	
	#J'enregistre dans un Dico
	dico_CDS[contig]=str(deb)+","+str(fin)+","+str(orientation)

	

print "\n\n nbr de lignes dans le fichier position des CDS : "+str(nombre)
	

#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------










#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------

### STEP 2 : Lecture du fichier PFAS et rognage des séquences.
print "\n\nSTEP 2 : Lecture du fichier PFAS "


# Je créé mes fichier de sortie
name1=output+"_UTR5.pfas" ; tmp_UTR5=open(name1,"w")
name2=output+"_CDS.pfas"  ; tmp_CDS=open(name2,"w")
name3=output+"_UTR3.pfas"  ; tmp_UTR3=open(name3,"w")


# Fonction permettant de faire le reverse complement
def multiple_replace(string, rep_dict):
    pattern = re.compile("|".join([re.escape(k) for k in rep_dict.keys()]), re.M)
    return pattern.sub(lambda x: rep_dict[x.group(0)], string)



# Lecture du fichier PFAS
toprint="no"
for line in open(pfas) :

	line=line.strip()
	
	#Si je suis sur une ligne de nom de contig, alors je récupère les infos et imprime la ligne telle quelle.
	if line.startswith(">") : 
		contig=line.split("|ind")[0]
		contig=contig.replace(">" , "")
		
		if contig in dico_CDS:
			
			#Alors je passe en mode toprint pour affichier la séquence du contig dans les fichiers de sortie.
			toprint="yes"
			
			#Je récupère les Infos.
			start=int(dico_CDS[contig].split(",")[0])
			end=int(dico_CDS[contig].split(",")[1])
			orientation=dico_CDS[contig].split(",")[2]

			#J'inscrit le nom dans les fichiers de sortie
			tmp_UTR5.write(line+"\n") ; tmp_CDS.write(line+"\n") ; tmp_UTR3.write(line+"\n")
		
		else : 
			
			# Alors je sors du mode toprint, pour ne pas enregistrer la séquence dans la sortie
			toprint="no"
			
			
	#Si je suis au niveau d'une séquence, alors je la oupe en 3 morceaux
	elif toprint=="yes" :
		sequence=line
		
		
		if orientation=="sens":

			sequence_UTR5=sequence[0:int(start)]
			sequence_CDS=sequence[int(start)-1:int(end)]
			sequence_UTR3=sequence[int(end):]


		if orientation=="antisens":
			
			#Pas facile a comprendre dans ce sens la. En gros, si les bornes donnée sont 338 puis 1918. Alors la partie de 1918 jusqu'à la fin = taille de l'UTR5. Mais je vais chercher cet UTR5 au début du fasta, sans faire de reverse complement.
			longueur=len(sequence)

			# Je coupe
			sequence_UTR5=sequence[0:longueur-end]
			sequence_CDS=sequence[longueur-end:longueur-start+1]
			sequence_UTR3=sequence[longueur-start+1:]

		
			
		#Puis j'ajoute au fichiers de sortie !
		tmp_UTR5.write(sequence_UTR5+"\n") ; tmp_CDS.write(sequence_CDS+"\n") ; tmp_UTR3.write(sequence_UTR3+"\n")

	

#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------









