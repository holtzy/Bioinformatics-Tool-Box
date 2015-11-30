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
parser.add_argument('-in', required=True, help='Liste des contigs pour lesquels on veut effectuer une modification')

args = parser.parse_args()

fic1=args.in




#------------------------------------------------------------------------------------------------------
#Fonctions préliminaires

def multiple_replace(string, rep_dict):
    pattern = re.compile("|".join([re.escape(k) for k in rep_dict.keys()]), re.M)
    return pattern.sub(lambda x: rep_dict[x.group(0)], string)
    		

def reverse_complement(sequence) :
	sequence_reversed=multiple_replace(sequence, {'A':'T', 'T':'A', 'C':'G', 'G':'C'})
	sequence_reversed=sequence_reversed[::-1]
	return(sequence_reversed)

#------------------------------------------------------------------------------------------------------

	



#------------------------------------------------------------------------------------------------------
#Activation des fonctions

sequence=""

for line in open(fic1) :
	line=line.strip()
		
	if line.startswith(">"):
		print line
		
	else:
		sequence=sequence+line
		

sequence=reverse_complement(sequence)
print sequence
		


