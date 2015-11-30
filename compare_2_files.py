#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------------------------------------------------------------------------------#
#
#    SCRIPT PYTHON COMPARE_2_FILES.PY
#     				Take 2 files. Return the number of lines of fic1 that ARE in fic2, then the number of lines of fic1 that are NOT in FIC2
#					
#  					Author : Yan Holtz, yan1166@hotmail.com
#-------------------------------------------------------------------------------------------------------------#




import os
import sys
import re
try:
  import argparse
except ImportError:
  print"oops, the import /argparse/ didn't work"


parser = argparse.ArgumentParser(description= '\n Take 2 files. Return the number of lines of fic1 that ARE in fic2, then the number of lines of fic1 that are NOT in FIC2\n')
parser.add_argument('-fic1', required=True, help=' file number 1')
parser.add_argument('-fic2', required=True, help=' file number 2 ')


args = parser.parse_args()

fic1=args.fic1
fic2=args.fic2






#------------------------------------------------------------------------------------------------------

print "\n\n--------\n\nBienvenu ! Merci de me faire taffer, je commencais à m'ennuyer !\n\n------------\n\n"		





#Script python permettant de comparer nos SNPs avec ceux de Wong
#Je rentre deux fichiers fichiers 1 et fichiers 2
#Il me donne le nombre de lignes de fichier 1 qui SONT dans fichier 2, puis le nombre de ligne qui n'y SONT PAS.


import sys
import re

dico_fic1=dict()
dico_fic2=dict()

for line in open(fic1):
	dico_fic1[line]=""
for line in open(fic2):
	dico_fic2[line]=""
	
nombre=0
nombre2=0
for i in dico_fic1:
	if i in dico_fic2:
		nombre+=1
	else:
		nombre2+=1
	
print nombre
print nombre2




print "\n\n--------\n\nEnd of Execution\n\n------------\n\n"		

#------------------------------------------------------------------------------------------------------


