#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON
#     				permet, pour une liste de nom donnée, de les modifier dans un fichier .fasta
#  					Yan Holtz, yan1166@hotmail.com
#-------------------------------------


import sys
import re
try:
  import argparse
except ImportError:
  print"oops, the import /argparse/ didn't work"


parser = argparse.ArgumentParser(description= '\n\nPermet de modifier le nom des contigs présent dans la liste d\'entrée dans un fasta\n\n')
parser.add_argument('-liste', required=True, help='Liste des contigs pour lesquels on veut effectuer une modification')
parser.add_argument('-input', required=True, help='Fichier .fasta dans lequel on veut modifier certains noms de contig')
parser.add_argument('-output', required=True, help='Fichier de sortie')

args = parser.parse_args()

liste_des_noms=args.liste
sequence_fasta=args.input
output=args.output


#------------------------------------------------------------------------------------------------------

### FUNCTION 0 : Création d'un dictionnaire à partir de la liste de nom.
def building_hash_table(liste_des_noms):
	"""	A partir de la liste des noms, on fabrique un dictionnaire"""
	id=1
	dictionnaire=dict()
	for line in open(liste_des_noms):
			if line not in dictionnaire:
				line=re.sub(r'(\n)',r"",line)
			 	dictionnaire[line]=id
			 	id=id+1
	#print dictionnaire
	return dictionnaire,id
	
#------------------------------------------------------------------------------------------------------
	
	
	
### FUNCTION 1 : changement du nom des contigs sélectionnés.
def changement_nom_sequence(sequence_fasta, dictionnaire):
	"""pour chaque sequence fasta présente dans le dico, j'effectue une modification"""
	open(output, "w")
	for line in open(sequence_fasta):
		line = line.rstrip()
		cluster=line.replace('>','')
		cluster=re.sub(r'(\|Contig[0-9]*)', r"", cluster)
		if cluster in dictionnaire:
			line=line+"|less_than_5_individuals"
		line=line+'\n'
		tmp=open(output,"a")
		tmp.write(line)
		tmp.close
		
#------------------------------------------------------------------------------------------------------
	
### UTILISATION FUNCTION 0 : creation du dictionnaire
[a,b]=building_hash_table(liste_des_noms)
dictionnaire=a
id=b
print dictionnaire

#------------------------------------------------------------------------------------------------------

### UTILISATION FUNCTION 1 : changement des noms
changement_nom_sequence(sequence_fasta, dictionnaire)


