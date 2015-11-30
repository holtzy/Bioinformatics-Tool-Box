#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON Trasnform_fasta_in_monoline.py
#     				Take a Fasta file. If sequences of a contig are split in several lines, it will paste them in one line
#  					Yan Holtz, yan1166@hotmail.com
#-------------------------------------

import os
import sys
import re
try:
  import argparse
except ImportError:
  print"oops, the import /argparse/ didn't work"


parser = argparse.ArgumentParser(description= '\n\nTake a Fasta file. If sequences of a contig are split in several lines, it will paste them in one line\n\n')
parser.add_argument('-in', required=True, help=' Fasta input ')

args = parser.parse_args()

fic1=args.in




#------------------------------------------------------------------------------------------------------
name=fic1+"_monoline"
tmp=open(name,"w")
debut="yes"

for line in open(fic1) :
	line=line.strip()
	
	if line.startswith(">") :
		
		#1/ Je finis le contig précédent (sauf dans le cas du premier)
		if debut != "yes":
			tmp.write(sequence+"\n")
			
		
		#2/ J'initialise le suivant.
		tmp.write(line+"\n")
		sequence=""
		debut="no"
		
		
	else :
		sequence=sequence+line


#Attention a ne pas oublier la derniere séquence !
tmp.write(sequence)


#Et j'écrase le fichier d'origine
commande="mv "+name+" "+fic1
os.system(commande)


print "\n----------------------\n"
print "done ! Chaque séquence de votre Fasta ne prend maintenant qu'une seule ligne !"
print "\n----------------------\n"

#------------------------------------------------------------------------------------------------------





