#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON Find_Info_precise_position_in_alr.py
#     				Ce script permet de sortir des stats pour des positions données d'un fichier ALR :
#						- couverture moyenne
#						- couverture max et mini
#						- Nombres d'indivs avec couverture >0 puis >10
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
parser.add_argument('-positions', required=False, help=' Position à étudier : format contig / tabulation / position')
parser.add_argument('-out', required=True, help=' fichier de sortie')


args = parser.parse_args()
alr=args.alr
output=args.out
positions=args.positions








#------------------------------------------------------------------------------------------------------

### STEP -1 : Dictionnaire des positions à étudier

#Si des positions sont indiquées, j'étudie ces positions
if positions is not None :

	nbr_de_position_au_total=0
	nbr_de_contigs_au_total=0
	pos_to_check=dict()
	
	for line in open(positions):
		nbr_de_position_au_total+=1
		line=line.split()
		contig=line[0].replace(">","")
		
		if contig not in pos_to_check:
			nbr_de_contigs_au_total+=1
			pos_to_check[contig]=line[1]
			
		else:
			pos_to_check[contig]=pos_to_check[contig]+","+line[1]
	 
	pos_redondante=nbr_de_position_au_total - len(pos_to_check)
	print "\n\n\n---"
	print str(nbr_de_position_au_total)+" positions sont présentes dans la liste donnée."+"\n"+"Ces positions concernent "+str(nbr_de_contigs_au_total)+" contigs différents au total."+"\n"+"Cependant, "+str(pos_redondante)+" positions données sont redondantes"
	print "---\n\n\n"


#Sinon j'étudie toutes les positions.
if positions is None :
	print_all="yes"
	print "Vous avez décidé d'analyser toutes les positions du fichier .alr donné !!!!!"
	pos_to_check=""

#------------------------------------------------------------------------------------------------------


 
 








 
 
#------------------------------------------------------------------------------------------------------

### STEP 2  -- Je parse mon fichier ALR a la recherche de mes positions !!

tmp=open(output,"w")
tmp.write("Contig"+"\t"+"Position"+"\t"+"Polymorphe"+"\t"+"Moyenne"+"\t"+"Max"+"\t"+"Min"+"\t"+"NumOver0"+"\t"+"NumOver10"+"\n")
nbr_reads=""

#Je parcours mon fichier ALR
for line in open(alr) :
	line=line.strip()
	
	#Lorsque je change de contig, je réinitialise tout
	if line.startswith(">") :
		liste=[]
		num=-1
		contig=line.replace(">","")
		if contig in pos_to_check:
			for i in pos_to_check[contig].split(","):
				liste.append(i)

				
	#Sinon je me balade a la recherche de ma position
	else:
		line=line.split("\t")
		num+=1
		
		#Si je suis dans une position ciblée:
		if str(num) in liste or ( print_all=="yes" and num>0 ) :
			vecteur=[]
			
			#Je récupère les infos des individus un par un
			for col in range(2,len(line)):
					if line[1]=="P":
						nbr_reads=line[col].split("[")[0]
						poly="P"
					if line[1]=="M":
						nbr_reads=line[col]
						poly="M"
					vecteur.append(int(nbr_reads))
			
			
			#Je calcule quelques stats sur le vecteur
			maxi=max(vecteur)
			mini=min(vecteur)
			tot=0
			for v in vecteur:
				tot=tot+int(v)
			moyenne=tot/len(vecteur)
			sup_0=len([elem for elem in vecteur if int(elem) > 0])
			sup_10=len([elem for elem in vecteur if int(elem) > 10])
			
			#J'inscris quelques stats dans le fichier de sortie
			tmp.write(contig+"\t"+str(num)+"\t"+poly+"\t"+str(moyenne)+"\t"+str(maxi)+"\t"+str(mini)+"\t"+str(sup_0)+"\t"+str(sup_10)+"\n")
		
		


print "\n\n\n---"
print("Fin du taff, Merci pour votre visite, votre fichier "+output+" est prêt, bonne lecture")
print "---\n\n\n"

	