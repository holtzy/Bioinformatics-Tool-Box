#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON get_majoritary_nucleotide_of_precise_positions_in_gen.py
#     				Ce script permet de sortir le génotype pour des positions données et pour un individu donné d'un fichier GEN :
#  					Yan Holtz, yan1166@hotmail.com
#-------------------------------------


import os
import sys
import re
try:
  import argparse
except ImportError:
  print"oops, the import /argparse/ didn't work"


parser = argparse.ArgumentParser(description= 'Ce script permet de sortir le génotype pour des positions données et pour un individu donné d un fichier GEN :')
parser.add_argument('-gen', required=True, help=' fichier .gen de résultat de reads2snp')
parser.add_argument('-positions', required=True, help=' Position à étudier : format contig / tabulation / position')
parser.add_argument('-ind', required=True, help=' numéro de colonne de l individu à étudier')
parser.add_argument('-out', required=True, help=' fichier de sortie')


args = parser.parse_args()
gen=args.gen
output=args.out
ind=args.ind
positions=args.positions








#------------------------------------------------------------------------------------------------------

### STEP -1 : Dictionnaire des positions à étudier
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
 
pos_redondante=nbr_de_position_au_total - len(pos_to_check) -1
print "\n\n\n---"
print str(nbr_de_position_au_total)+" positions sont présentes dans la liste donnée."+"\n"+"Ces positions concernent "+str(nbr_de_contigs_au_total)+" contigs différents au total."+"\n"+"Cependant, "+str(pos_redondante)+" positions données sont redondantes"
print "---\n\n\n"

#------------------------------------------------------------------------------------------------------


 
 








 
 
#------------------------------------------------------------------------------------------------------

### STEP 2  -- Je parse mon fichier GEN a la recherche de mes positions !!

tmp=open(output,"w")
tmp.write("Contig"+"\t"+"Position"+"\t"+"Nb_homo"+"\t"+"Nb_hetero"+"\t"+"Pourc_hetero"+"\t"+"Nb_missing"+"\t"+"homo_maj"+"\t"+"nb_A"+"\t"+"nb_C"+"\t"+"nb_G"+"\t"+"nb_T"+"\t"+"all_selected_ind"+"\n")


#Je parcours mon fichier GEN
for line in open(gen) :
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
		if str(num) in liste:
			vecteur=[]
			nb_homo=0
			nb_hetero=0
			nb_missing=0
			nb_A=0
			nb_C=0
			nb_G=0
			nb_T=0
			genot_studied="-"
			
			#Je récupère les infos des individus un par un
			for col in range(1,len(line)):
				a=line[col].split("|")[0]
				vecteur.append(a)
				
				# Récup l'individu étudié
				if col==int(ind):
					genot_studied=a
			
			#Je calcule des stats sur le vecteur obtenu
			for i in vecteur:
				if i=="NN":
					nb_missing+=1

				else:
					if i[0]==i[1]:
						nb_homo+=1
					else:
						nb_hetero+=1
				
				if i=="AA":
					nb_A+=1
				if i=="CC":
					nb_C+=1
				if i=="GG":
					nb_G+=1
				if i=="TT":
					nb_T+=1
			
			tot=nb_homo+nb_hetero
			if tot !=0:
				pourc_hetero=float(nb_hetero)/float(nb_homo+nb_hetero)*100
			else:
				pourc_hetero="-"
				
			# Récup du nucl majoritaire:
			all_nucl=["AA","CC","GG","TT"]
			a=[nb_A, nb_C, nb_G, nb_T]
			m=max(a)
			my_num=[i for i, j in enumerate(a) if j == m]
			nucl=all_nucl[int(my_num[0])]
			# Si 2 nucl égualité ou alors aucun allele présent --> je met des tiret
			if len(my_num)>1 or m==0:
				nucl="-"
				
			#J'inscris quelques stats dans le fichier de sortie
			tmp.write(contig+"\t"+str(num)+"\t"+str(nb_homo)+"\t"+str(nb_hetero)+"\t"+str(pourc_hetero)+"\t"+str(nb_missing)+ "\t"+str(nucl)+"\t"+str(nb_A)+"\t"+str(nb_C)+"\t"+str(nb_G)+"\t"+str(nb_T)+"\t"+str(genot_studied)+"\n")
		
		


print "\n\n\n---"
print("Fin du taff, Merci pour votre visite, votre fichier "+output+" est prêt, bonne lecture")
print "---\n\n\n"

	