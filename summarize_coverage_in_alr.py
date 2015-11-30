#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON Summarize coverage in alr
#     				Ce script permet de récupérer les couvertures de chaque position d'un fichier alr pur une liste de contig donnée
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
parser.add_argument('-liste', required=False, help=' liste des contigs a étudier. Si aucun liste n\'est donnée, tous les contigs seront étudiés.')
parser.add_argument('-out', required=True, help=' fichier de sortie')


args = parser.parse_args()
alr=args.alr
output=args.out
liste=args.liste








#------------------------------------------------------------------------------------------------------

### STEP -1 : Dictionnaire des contigs à étudier. Si aucune liste n'est donné, alors je mets tous les contigs dans la liste

contigs_to_check=dict()

if liste is not None :
	for line in open(liste):
 		line=line.strip()
 		contig=line.replace(">","")
 		contigs_to_check[contig]=""
 
if liste is None :
	for line in open(alr) :
		line=line.strip()
		if line.startswith(">") :
			contig=line.replace(">","")
 			contigs_to_check[contig]=""
 			
 			
 			
#------------------------------------------------------------------------------------------------------


 
 




 
 
#------------------------------------------------------------------------------------------------------

### STEP 2  -- Je parse mon fichier ALR a la recherche de mes contigs !!

tmp=open(output,"w")
tmp.write("Contig"+"\t"+"Position"+"\t"+"Moyenne"+"\t"+"Max"+"\t"+"Min"+"\t"+"NumOver0"+"\t"+"NumOver10"+"\n")

nbr_reads=""
to_print=""
num=0

#Je parcours mon fichier ALR
for line in open(alr) :
	line=line.strip()
	num+=1
	
	#Lorsque je change de contig, je réinitialise tout
	if line.startswith(">") :
		num=-1
		contig=line.replace(">","")
		if contig in contigs_to_check:
			to_print="yes"
		else:
			to_print="no"
				
	#Sinon je me balade a la recherche de ma position
	elif to_print=="yes" and num>0:
		line=line.split("\t")
		vecteur=[]
			
		#Je récupère les infos des individus un par un
		for col in range(2,len(line)):
			if line[1]=="P":
				nbr_reads=line[col].split("[")[0]
			if line[1]=="M":
				nbr_reads=line[col]
			vecteur.append(int(nbr_reads))
		
		
		#Je calcule quelques stats sur le vecteur
		maxi=max(vecteur)
		mini=min(vecteur)
		tot=0
		for v in vecteur:
			tot=tot+int(v)
		moyenne=float(tot)/float(len(vecteur))
		moyenne=round(moyenne,1)
		sup_0=len([elem for elem in vecteur if int(elem) > 0])
		sup_10=len([elem for elem in vecteur if int(elem) > 10])
		
		#J'inscris quelques stats dans le fichier de sortie
		tmp.write(contig+"\t"+str(num)+"\t"+str(moyenne)+"\t"+str(maxi)+"\t"+str(mini)+"\t"+str(sup_0)+"\t"+str(sup_10)+"\n")
		
		


print "\n\n\n---"
print("Fin du taff, Merci pour votre visite, votre fichier "+output+" est prêt, bonne lecture")
print "---\n\n\n"

	