#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT Tranform_bayt_format
#     				Pour transformer les baytes comme les veux sylvain santoni
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


parser = argparse.ArgumentParser(description= 'Permet de créer les baytes en plusieurs format.')
parser.add_argument('-input1', required=True, help=' fichier contenant les baytes provenant de Dic2 * Silur')
parser.add_argument('-input2', required=False, help=' fichier contenant les baytes provenant de EPO')
parser.add_argument('-fasta', required=True, help=' fasta de EPO')
parser.add_argument('-output', required=True, help=' fichier de sortie')


args = parser.parse_args()
input1=args.input1
input2=args.input2
fasta=args.fasta
output=args.output






#ETAPE DE CE SCRIPT : 
	# Step 0 : Réalisation du dico des contigs de EPO
	
	# Step 1 : Mise en place de la fonction de traitement :
		# -- 0 -- : Sélection des séquences qui valent le coup de fabriquer un bayte 
		# -- 1 -- : type 1 : Le SNP est centré est la séquence fait 200 pb : Dans ce cas la je sors 4 séquences de 120pb : 2 au centres et 2 décalées.
		# -- 2 -- : type 2 : Le SNP n'est pas centré mais situé entre 80 et 120 pb + séquence fait 200 pb en tout, alors je créé 2 séquences décalées seulement
		# -- 3 -- : type 3 et 4 : Le SNP est en position extrême (<80 ou > 120) et séquence = 200pb , alors :
						#--> Si j'ai divergence entre le contig en question et son homéologue : je fais un type 4 : uniquement l'allèle de Dic2 répété 2 fois
						#--> Si je n'ai pas divergence : type 3, une séquence pour Dic2, une pour Silur
		# -- 4 -- : type 3 : Le SNP est centrée et la séquence fait 120 pb > alors je garde exactement cette séquences, j'obtiens 2 séquences alignées (une Dic2, une Silur)
		# -- 5 -- : type 3 : Le SNP n'est PAS centrée et la séquence fait 120 pb > alors je garde exactement cette séquences, j'obtiens 2 séquences alignées (une Dic2, une Silur)
	 
	# Step 2 : Utilisation des fonctions sur tous les fichiers d'entrées






#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------

# STEP 0 : dico du fasta EPO : il permettra juste de voir si il y  a divergence entre homéocontigs pour les types 4
dico_fasta=dict()

for line in open(fasta):
	line=line.strip()
	
	if line.startswith(">"):
		contig=line.replace(">","")
	
	if contig in dico_fasta :
		dico_fasta[contig]=dico_fasta[contig]+line
	else :
		dico_fasta[contig]=""

#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------











#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------

# STEP 0 : Fonction qui lit les baytes et les transforme en type 1,2,3 ou 4 en fonction de leur caractéristiques

def transform_bayte_format(fic , nb_seq_tot , nb_snp_initial , 	nb_snp_lost_car_redondant , nb_snp_final , tag):
	for line in open(fic) :
		
		line=line.strip()
		
		#---0/ Je commence par décider si la séquence lue vaut le coup de fabriquer un baytes ou non :
		
		#Pour chaque ligne d'en tete je récupère tous les infos du SNP
		if line.startswith(">"):
			SNP=line.split("|All")[0]
			SNP=SNP.replace(">","")
			allele=line.split("Allele_")[1]
			allele=allele.split("|")[0]
			position=line.split("pos_")[1]
			position=int(position.split("|")[0])
			contig=line.split("@")[0].replace(">","")
			entete=line
			cluster=line.split("|")[0]
			
			if allele=="Dic2":
				nb_snp_initial+=1
				
		#Je fais le taff que si ce cluster n'a pas déja été trouvé pour un AUTRE contig --> pour Silur
		if allele=="Silur" and not line.startswith(">") :
			if cluster in dico_cluster_sil and contig not in dico_contig_sil	:
				nb_snp_lost_car_redondant+=0.5
				continue
			else :
				dico_cluster_sil[cluster]="trouve"
				dico_contig_sil[contig]="--"
					
		#Je fais le taff que si ce cluster n'a pas déja été trouvé pour un AUTRE contig --> pour Dic2
		if allele=="Dic2" and not line.startswith(">") :
			if cluster in dico_cluster_dic and contig not in dico_contig_dic :
				nb_snp_lost_car_redondant+=0.5
				continue
			else :
				dico_cluster_dic[cluster]="trouve"	
				dico_contig_dic[contig]="--"
	
	
	
	
	
		#Je travail uniquement sur les séquences = lignes non en tete
		if not line.startswith(">") :
		
	
			#J'ai donc un SNP de plus ajouté a la liste
			nb_snp_final+=0.5
		
			
			#---1/ Si le SNP est centré et que le baytes fait 200 pb:
			if position == 101 and len(line)==200:
			
				#Je mettrai l'allele de dic2 sur la gauche et celui de Silur sur la droite
				if allele=="Dic2":
					side=line[0:120]
				else:
					side=line[80:200]
				
				#Je créé un bayte centré et un décalé :
				tmp.write(entete+"|type1@middle|"+tag+"\n")
				tmp.write(line[40:160]+"\n")
				tmp.write(entete+"|type1@side|"+tag+"\n")
				tmp.write(side+"\n")
				nb_seq_tot+=2
		
		
		
			#---2/ Si 80 < position SNP < 120 et longueur = 200:
			if position < 120 and position > 80 and position != 101 and len(line)==200 : 
				
				#Je mettrai l'allele de dic2 sur la gauche et celui de Silur sur la droite
				if allele=="Dic2":
					side=line[position-80-1:position-1+40]
				else:
					side=line[position-40-1:position-1+80]
		
				#Je créé un bayte  décalé seulement:
				tmp.write(entete+"|type2|"+tag+"\n")
				tmp.write(side+"\n")
				nb_seq_tot+=1
				
		
		
		
		
			#---3/ SNP en position extreme
			if (position <= 80 or position >= 120) and len(line)==200 : 
				
				#Le mode de base est "two", on va voir si on peux le passer en "one"
				#Le mode one correspond aux cas ou il y a divergence entre les 2 homéocontigs : dans ce cas on mettra uniquement un allèle pour le baytes
				type="two"
		
				#Si je suis dans un cas EPO
				if not entete.startswith(">Traes") :
					contig=SNP.split("@")[0]
					pos_SNP=int(SNP.split("@")[1])
					contig_short=contig.split("|")[0]+"|"+contig.split("|")[1]
					nature=contig.split("|")[-1]
					
					if nature=="complementarySeq":
						lautre=contig_short+"|"+"likelySeq"
					if nature=="likelySeq":
						lautre=contig_short+"|"+"complementarySeq"
					if nature=="original":
						lautre=""
						
					#Je récupère la séquence de l'autre, coupée au bon endroit
					if lautre in dico_fasta  and not line.startswith(">"):
						seq=dico_fasta[lautre]
						seq=seq[pos_SNP-position : 	pos_SNP-position+200]
						
						#Je compte le nombre de différence entre les 2 séquences?
						nb_diff=0
						for i in range(1,200):
							if line[i-1]!= seq[i-1] :
								nb_diff+=1
								
						#Si j'ai plus de 8 différences, alors je me mets en mode un seul baytes.
						if nb_diff > 8 :
							type="one"	
		
							
				#Maintenant j'imprime en fonction
				
				#Cas du "one"
				if type=="one" :
					if allele=="Dic2" :
						tmp.write(entete+"|type4|"+tag+"\n")
						if position <= 80:
							tmp.write(line[0:120]+"\n")
							tmp.write(entete+"|type4|"+tag+"\n")
							tmp.write(line[0:120]+"\n")
							nb_seq_tot+=2
						if position >= 120:
							tmp.write(line[80:200]+"\n")
							tmp.write(entete+"|type4|"+tag+"\n")
							tmp.write(line[80:200]+"\n")
							nb_seq_tot+=2
		
				#Cas du "two"
				if type=="two" :
					tmp.write(entete+"|type3|"+tag+"\n")
					if position <= 80:
						tmp.write(line[0:120]+"\n")
						nb_seq_tot+=1
					if position >= 120:
						tmp.write(line[80:200]+"\n")
						nb_seq_tot+=1
	
	
	
	
			#---4/ Si le SNP est centré et que le baytes fait 120 pb :
			if position == 61 and len(line)==120 : 
				
				#Je créé un bayte centré :
				tmp.write(entete+"|type3|"+tag+"\n")
				tmp.write(line+"\n")
				nb_seq_tot+=1
				continue
		
			
			#---5/ Si le SNP n'est PAS centré et que le baytes fait 120 pb:
			elif len(line)==120  : 
				#Je créé un bayte non centré :
				tmp.write(entete+"|type3|"+tag+"\n")
				tmp.write(line+"\n")
				nb_seq_tot+=1
		
		
			
	#Je renvoi le nbr de sequence totales récupérées
	return(nb_seq_tot , nb_snp_initial , 	nb_snp_lost_car_redondant , nb_snp_final )


#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------









#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------

# STEP 2 : Transformation des baytes !
tmp=open(output,"w")

#Nombre total de séquence ajoutée au fichier de sortie :
nb_seq_tot=0
nb_snp_initial=0
nb_snp_lost_car_redondant=0
nb_snp_final=0


#Je vais faire en sorte de ne jamais avoir 2 contigs différents provenant du même cluster...
dico_cluster_dic=dict() ; dico_cluster_sil=dict()
dico_contig_dic=dict()  ; dico_contig_sil=dict()


#J'envoie la fonction pour les baytes Dic2 * Silur pur
result=transform_bayte_format(input1 , nb_seq_tot , nb_snp_initial , 	nb_snp_lost_car_redondant , nb_snp_final  , "from_DicSil")
nb_seq_tot=result[0]
nb_snp_initial=result[1]
nb_snp_lost_car_redondant=result[2]
nb_snp_final=result[3]


#J'affiche le nombre de baytes total avant ajout des SNP provenant de EPO :
print "\n------------\nUtilisation des SNPs de Dic2 * Silur\n----------------"
print "\nle nombre de SNP présent dans le fichier de baits d'entrée est : "+str(nb_snp_initial)+"\n"					
print "\nle nombre de SNP perdus car redondant de cluster est : "+str(nb_snp_lost_car_redondant)+"\n"					
print "\nle nombre de SNP présents avant ajout des SNP en provenance de EPO est : "+str(nb_snp_final)+"\n"					
print "\nle nombre de sequence présentes avant ajout des SNP en provenance de EPO est : "+str(nb_seq_tot)+"\n"					
			

if input2 is not None:
	#J'envoie la fonction pour les baytes provenant des SNPs EPO pour compléter
	nb_snp_initial=0
	nb_snp_lost_car_redondant=0
	result=transform_bayte_format(input2 , nb_seq_tot , nb_snp_initial , 	nb_snp_lost_car_redondant , nb_snp_final , "from_EPO")
	nb_seq_tot=result[0]
	nb_snp_initial=result[1]
	nb_snp_lost_car_redondant=result[2]
	nb_snp_final=result[3]
	
	
	#J'affiche le nombre de baytes total avant ajout des SNP provenant de EPO :
	print "\n------------\nUtilisation des SNPs de EPO pour compléter \n----------------"
	print "\nle nombre de SNP présent dans le fichier de baits d'entrée est : "+str(nb_snp_initial)+"\n"					
	print "\nle nombre de SNP perdus car redondant de cluster est : "+str(nb_snp_lost_car_redondant)+"\n"					
	print "\nle nombre de SNP présents maintenant, avec l'ajout des SNP en provenance de EPO est : "+str(nb_snp_final)+"\n"					
	print "\nle nombre de sequence présentes maintenant, avec l'ajout des SNP en provenance de EPO est : "+str(nb_seq_tot)+"\n"					























