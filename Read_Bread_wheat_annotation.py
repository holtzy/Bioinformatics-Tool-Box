#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT Read_Bread_Wheat annotation.py
#
#     				Les positions des UTR + des exons sont donnée dans les entetes d'un fichier fasta dans un format pas évident a déchiffrer. Ce script va permettre de lire ce fichier de ressortir 2 fichiers claires donnant les positions des CDS (1), et la position des exons (2)
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


parser = argparse.ArgumentParser(description= 'Les positions des UTR + des exons sont donnée dans les entetes d\'un fichier fasta dans un format pas évident a déchiffrer. Ce script va permettre de lire ce fichier de ressortir 2 fichiers claires donnant les positions des CDS (1), et la position des exons (2)')
parser.add_argument('-annot', required=True, help=' fichier d\'annotation de ensembl a traiter')
parser.add_argument('-output_pos_CDS', required=False, help=' fichier de sortie CDS : fichier texte : contig,début CDS, fin CDS, orientation')
parser.add_argument('-output_pos_exons', required=False, help=' fichier de sortie Exon : fichier texte : contig,exon1,exon2...')


args = parser.parse_args()
annot=args.annot
out1=args.output_pos_CDS
out2=args.output_pos_exons

#------------------------------------------------------------------------------------------------------









#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 1 : Collecte des positions des Introns + UTR dans le cas du blé tendre dans un dictionnaire
print "-----\n\nLecture du fichier\n\n----\n"

dico_intron=dict()
dico_UTR=dict()

nombre=0 


#A chaque nouvelle ligne je récupère toutes les infos : contig, positions des exons
for line in open(annot):
	line.strip()
	
	if line.startswith(">Traes"):
		nombre=nombre+1
		line=line.split("|")

		
		#Je vire les cas chelou :
		#ou je n ai pas d'utr ni de CDS (comment c est possible je ne sais pas)
		#Ou j'ai des UTR décomposé
		if len(line) < 7:
			continue
		if len(line[6].split(";"))>1 :
			continue
		if len(line)>10 and len(line[11].split(";"))>1 :
			continue
		

		#Je récupère les principales infos.
		contig=line[0].replace(">","")
		strand=line[5]
		deb_gene=line[3]
		fin_gene=line[4]


		#Quelle est la longueur du CDS ???
		longueur_CDS=0
		liste_deb=line[8].split(";")
		liste_fin=line[9].split(";")
		for i in range(0,len(liste_deb)) :
			longueur_CDS=longueur_CDS +  int(liste_fin[i]) - int(liste_deb[i]) + 1
		
		
		#Ou finit le dernier exon? Attention, c'est une position intra exon. En fait ca revient a donner la taille de l'exon
		end_dernier_exon= int(max(line[9].split(";"))) 
		
					
		
		#-----------Sens normal
		if strand == "1":
		
			_sens_="sens"
		
			#___En premier on coupe au dernier nucléotide de l'UTR5. Si pas d'UTR5, on met donc 0
			#Si il n'y a pas d'UTR5
			if line[7]=="":
				fin_UTR5=0
			#Sinon je calcule la taille de l'UTR5 !!! (fin - début +1)
			else:
				fin_UTR5 = int(line[7]) - int(line[6]) + 1
				
			#Je peux alors faire mon dictionnaire d'exons!
			for i in range(0,len(liste_deb)) :
				a = fin_UTR5 + int(liste_deb[i])
				b= fin_UTR5 + int(liste_fin[i])
				if contig not in dico_intron:
					dico_intron[contig]=str(a)+".."+str(b)
				else :
					dico_intron[contig]=dico_intron[contig]+","+str(a)+".."+str(b)
			
				
			#___En deuxième on coupe à la fin de tous les exons (+1 pour etre sur le premier nucléotide de l'UTR3)
			fin_des_exons=fin_UTR5 + longueur_CDS +1
			
			
			#___et en troisième on calcule ou est la fin de l'UTR3 ?! (ca devrait donner la taille de la séquence dans le fasta). on vérifie si il y a un UTR3
			#Si il n'y a pas d'UTR3
			if len(line) < 11:
				len_UTR3=0
			#Sinon je calcule la taille de l'UTR3 !!! (fin - début +1)
			else:
				len_UTR3 = int(line[11]) - int(line[10]) + 1
			#et donc la position sur le contig	
			fin_UTR3 = fin_des_exons + len_UTR3 -1		
					
	
	
			dico_UTR[contig]=str(fin_UTR5+1)+"\t"+str(fin_des_exons-1)+"\t"+_sens_





		#-----------reverse complement
		if strand == "-1":
	
			_sens_="antisens"
			
			#En premier on coupe au dernier nucléotide de l'UTR3. Si pas d'UTR3, on met donc 0
			#Si il n'y a pas d'UTR3
			if len(line) < 11:
				fin_UTR3=0
			#Sinon je calcule la taille de l'UTR3 !!! (fin - début +1)
			else:
				fin_UTR3 = int(line[11]) - int(line[10]) + 1

			#Je peux alors faire mon dictionnaire d'exons!
			for i in range(0,len(liste_deb)) :
				a = fin_UTR3 + int(liste_deb[i])
				b= fin_UTR3 + int(liste_fin[i])
				if contig not in dico_intron:
					dico_intron[contig]=str(a)+".."+str(b)
				else :
					dico_intron[contig]=dico_intron[contig]+","+str(a)+".."+str(b)

				
			#En deuxième on coupe à la fin de tous les exons (+1 pour etre sur le premier nucléotide de l'UTR5)
			fin_des_exons=fin_UTR3 + longueur_CDS +1		
	
			
			#___et en troisième on calcule ou est la fin de l'UTR5 ?! (ca devrait donner la taille de la séquence dans le fasta). on vérifie si il y a un UTR5
			#Si il n'y a pas d'UTR5
			if line[7]=="":
				len_UTR5=0
			#Sinon je calcule la taille de l'UTR5 !!! (fin - début +1)
			else:
				len_UTR5 = int(line[7]) - int(line[6]) + 1
			#Et donc la position sur le contig
			fin_UTR5 = fin_des_exons + len_UTR5 - 1

	
			#et je rempli mon dico UTR
			dico_UTR[contig]=str(fin_UTR3+1)+"\t"+str(fin_des_exons-1)+"\t"+_sens_


print "-----\n\nLe fichier a été lu correctement. Nombre de séquences annotées : "+str(nombre)+"\n\n----\n"

#------------------------------------------------------------------------------------------------------





#------------------------------------------------------------------------------------------------------

#Impression de ces 2 Dicos :

tmp1=open(out1,"w")
tmp2=open(out2,"w")

for contig in dico_UTR :
	tmp1.write(contig+"\t"+dico_UTR[contig]+"\n")

for contig in dico_intron :
	tmp2.write(contig+"\t"+dico_intron[contig]+"\n")




#------------------------------------------------------------------------------------------------------



