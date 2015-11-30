#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON FIND_BEST_BLAST
#     				Take a blastn output, determine for each LEFT sequence the best blast made on the RIGHT sequence. 
#					Warning, the format of the blastn output MUST be  : -outfmt '6 qseqid sseqid qstart qend qlen sstart send slen length pident evalue'
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


parser = argparse.ArgumentParser(description=  '\n\nA partir d\'un résultat de blast, détermine pour chaque séquence de la colonne de GAUCHE, quel est sa meilleure séquence de DROITE. \nATTENTION, le blastn doit etre effectué avec es parametres de sortie -outfmt \'6 qseqid sseqid qstart qend qlen sstart send slen length pident evalue\'\n\n')
parser.add_argument('-blast', required=True, help='resultat de blastn')
parser.add_argument('-similarity', required=True, help=' similarity threshold')
parser.add_argument('-overlap_abs', required=True, help=' overlap threshold (absolute value)')
parser.add_argument('-overlap_rel', required=True, help=' overlap threshold (relative value)')
parser.add_argument('-out', required=True, help=' output file')


args = parser.parse_args()
fic_blast=args.blast
p=args.similarity
overlap_abs=args.overlap_abs
overlap_rel=args.overlap_rel
out=args.out







#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
dico_best_blast=dict()

#Je parse mon fichier de résultat de blast :
for line in open(fic_blast):
	
	#Récup des données :
	infos=line.split()
	id1=infos[0]; id2=infos[1]; deb1=int(infos[2]) ; end1=int(infos[3]) ; len1=int(infos[4]) ;pc_ident=infos[9]; len_align=int(infos[8]) ; evalue=infos[10]	

	#quelle est le % de la taille de la séquence sur laquelle on a le blast ?
	pourc = (end1 - deb1)*100 / len1
	
	
	#Je taff que pour les blasts qui respectent mes conditions
	if pc_ident >= float(p) and len_align > int(overlap_abs) and id1!=id2 and float(pourc) > float(overlap_rel):		
		
		#Si la séquence 1 n'a as été trouvée jusque la, alors c'est forcément son meilleur blast !
		if id1 not in dico_best_blast:	
			dico_best_blast[id1] = id2+","+evalue+","+pc_ident
			continue
		
	
		#Sinon je dois checker si ce nouveau blast est mieux ou moins bien que le précédent (comparaison a l'aide de la evalue)
		if id1 in dico_best_blast :
			old_evalue=dico_best_blast[id1].split(",")[1]
		
			old_pc_ident=dico_best_blast[id1].split(",")[2]

			#Cas 1 : la evalue est mieux que la précédente, alors je remplace mon best blast par celui ci ! :
			if float(old_evalue) > float(evalue) :
				dico_best_blast[id1] = id2+","+evalue+","+pc_ident
				
			#Cas 2 : la evalue est égale à la précédente, alors je prend le meilleure % de similarité :
			if old_evalue == evalue and float(old_pc_ident) < float(pc_ident) :
				dico_best_blast[id1] = id2+","+evalue+","+pc_ident
	
			#Cas 3 : la evalue est moins bien que la précédente, alors je ne fais rien


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Impression de la liste des best blasts

tmp=open(out , "w" )

for seq in dico_best_blast:
	tmp.write(seq+"\t"+dico_best_blast[seq].split(",")[0] + "\n")
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
