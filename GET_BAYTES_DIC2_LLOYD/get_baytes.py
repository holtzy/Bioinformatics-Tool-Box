#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON get_baytes
#     				Ce script permet de récupérer des baytes à partir de :
#						-un fichier de SNP
#						-un fichier donnant l'emplacement des introns dans les contigs
#						-un fichier donnant l'emplacement des CDS dans les contigs
#						-un fichier donnant les séquences de tous les contigs
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


parser = argparse.ArgumentParser(description= 'get the baits from SNPs, fasta and annotation. ')
parser.add_argument('-fasta', required=True, help='	fasta de EPO')
parser.add_argument('-position_intron', required=False, help='	position des introns dans les contigs. Format : Contig, debut exon1..fin exon1,debut exon2..finexon2, etc')
parser.add_argument('-position_CDS', required=False, help=' fichier indiquant la position des CDS dans les contigs. Format : contig,debut CDS, fin CDS')
parser.add_argument('-SNP', required=True, help='	fichier de SNP')
parser.add_argument('-out', required=True, help=' fichier de sortie')


args = parser.parse_args()
fasta=args.fasta
position_intron=args.position_intron
SNP=args.SNP
position_CDS=args.position_CDS
out=args.out




#Fichier d'information sur le déroulement du script.
tmpfichierinfo=open("fichier_info_baytes.txt","w")

#------------------------------------------------------------------------------------------------------










#------------------------------------------------------------------------------------------------------
### STEP 0 : Collecte de l'ensemble des contigs dans un dictionnaire
print "\n\nSTEP 0 : dictionnaire des contigs"
dico_fasta=dict() 
dico_name=dict()
nombre=0

for line in open(fasta):
	line=line.strip()
	if line.startswith(">"):
		nombre=nombre+1
		contig=line.replace(">","")
		
		#Attention si c'est un contig de blé tendre je dois raccourcir le nom
		if contig.startswith("Traes"):
			contig=contig.split("|")[0]
		
		dico_fasta[contig]=""
		
	else:
		dico_fasta[contig]=dico_fasta[contig]+line


	# Petite parenthèse : Pour EPO, l'annotation a été faite avant homéosplitter, donc avec des noms de contigs raccourcis. Je dois donc faire une table de lien...
	contig_short=contig.replace("|original" , "")
	contig_short=contig_short.replace("|likely" , "")
	contig_short=contig_short.replace("|complementary" , "")
	dico_name[contig_short]=contig


print "\n\n      ---> done !     "+str(nombre)+"     Contigs ont été récupérés dans le fichier fasta"
toprint= "nbr de contigs récupérés dans le fichier Fasta de référence: "+str(nombre)
tmpfichierinfo.write(toprint+"\n")

#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------







#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 1 : Collecte des positions des introns dans un dictionnaire.
dico_intron=dict()

if position_intron is not None:

	# Message d'avis + initialization des variables	
	print "\n\nSTEP 1 : dictionnaire des positions des introns"
	dico_blast=dict()
	nombre=0
	
	
	# A chaque nouvelle ligne je récupère toutes les infos : contig, scores de blast, et position des introns
	for line in open(position_intron):
		nombre=nombre+1
		line=line.strip()
		line=line.split()

		contig=line[0]
		positions=line[1]
		
		dico_intron[contig]=positions
		
	
	
	print "\n\n      ---> done  ! nbr de lignes dans le fichier position_intron : "+str(nombre)+"  cas" 
	toprint = "nbr de lignes dans le fichier position des introns : "+str(nombre)
	tmpfichierinfo.write(toprint+"\n")


#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------












#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 2 : Collecte des positions des UTR.
print "\n\nSTEP 2 : dictionnaire des positions des  UTR"

if position_CDS is not None:

	dico_UTR=dict()
	nombre=0 
	
	
	#A chaque nouvelle ligne je récupère toutes les infos : contig, positions des exons
	for line in open(position_CDS):
		line.strip()
		line=line.split()
		nombre+=1
		
		#Je récupère les infos : il s'agit des début et fin des CDS
		contig_short=line[0]
		deb=int(line[1])-1
		fin=int(line[2])+1
		
		#J'enregistre : attention, c'est les positions des UTR que j'inscrit, donc j'ajoute et enleve 1 ! Mais j'enregistre que si le contig est dans le fasta du début !
		if contig in dico_name :
			contig=dico_name[contig_short]
			longueur=len(dico_fasta[contig])
			dico_UTR[contig]="1.."+str(deb)+","+str(fin)+".."+str(longueur)

		
	
	print "\n\n      ---> done  ! nbr de lignes dans le fichier position_intron : "+str(nombre)+"  cas" 
	toprint = "nbr de lignes dans le fichier position des UTR : "+str(nombre)
	tmpfichierinfo.write(toprint+"\n")
	

#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------








#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 3 : Collecte des infos des SNPs
print "\n\nSTEP 3 : dictionnaire des positions des SNPs"

dico_SNP=dict() 
nombre=0
 
 
for line in open(SNP):
	
	# initialisation
	nombre=nombre+1
	line=line.strip()
	line=line.split()

	# Nom du contig et du SNP 
	contig=line[6].replace(">","")
	if contig.startswith("Trae"):
		contig=contig.split("|")[0]
	position=line[7] 
	SNP_name=contig+"@"+position
	
	
	# Je récupère les 2 allèles : je dois parser toute la ligne pour trouver les 2 allèles présent
	nucl1="-" ; nucl2="-"
	for i in range(10,len(line)+1 , 2) :
		nucl=line[i-1].split("|")[0][0]
		
		
		if nucl=="-":
			continue
		if nucl!=nucl1 and nucl1=="-":
			nucl1=nucl
		if nucl != nucl1 and nucl1!="-":
			nucl2=nucl
		
	# Enregistrement dans un dico
	dico_SNP[SNP_name]=str(position)+","+nucl1+","+nucl2		
			
			


print "\n\n      ---> done !    "+str(nombre)+"   SNPs ont été récupérés !!"
toprint = "nbr SNPs potentiellement utilisables pour faire des baytes : "+str(nombre)
tmpfichierinfo.write(toprint+"\n")


#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
















#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 4 : Synthèse des baits
print "\n\nSTEP 4 : Synthèse des baytes !! "
tmp=open(out,"w")
num_bayte=0
num_SNP_avec_annot=0


#Je procède SNP par SNP, pour voir a chaque fois si il est utilisable ou non. J'ai plein de critère à chécker :
for SNP_name in dico_SNP:
	
	
	#Récupération de toutes les infos du SNP et du contig :
	position_SNP=int(dico_SNP[SNP_name].split(",")[0])
	contig=SNP_name.split("@")[0]
	sequence=dico_fasta[contig]
	allele1=dico_SNP[SNP_name].split(",")[1]
	allele2=dico_SNP[SNP_name].split(",")[2]

	
	# Si le contig a des introns connus ou des UTR connus , alors je peux essayer de faire un bait, sinon je laisse tomber
	if contig in dico_intron or contig in dico_UTR:

		#Je marque un SNP qui a une annotation.
		num_SNP_avec_annot+=1
		
		#Je vais faire une grande liste de toutes les tranches ou on peut faire un bait = les exons + les UTRs
		if contig in dico_intron : 	
			position=dico_intron[contig]
		if contig in dico_UTR :
			if position is None :
				position=position+","+dico_UTR[contig]
			else :
				position=dico_UTR[contig]
		
		#Je vais parcourir ces tranches une par une et voir si il y en a une qui correspond au SNP
		position=position.split(",")
		for exon in position :
		
		
			#print contig
			#print exon
			#print dico_intron[contig]
			#print dico_UTR[contig]
			
			#Je récupère les borners de l'exon ou de l'UTR
			fin=int(exon.split("..")[1])
			debut=int(exon.split("..")[0])
			

			#J'ai récupéré toutes les infos, je peux poser mes conditions pour voir si mon SNP permet de faire un baytes!!!
			#Petite sortie pour voir les caractéristiques du bayte
			#Condition : au moins 10 nucléotides avant et après le SNP sur l'exon, et une taille d'exon supérieur ou égale à 120pb
			
			
			# Puis je faire un bait de 120pb de long ? C'est a dire que l'exon est >=120 et que le SNP est écarté d'au moins 10 bases des extrémités
			if position_SNP>(debut+10) and position_SNP<(fin-10) and (fin-debut)>=121 :
			
				#Alors je compte un bait en plus
				num_bayte+=1

				#Puis je mettre le SNP au milieu (au moins 60pb d'exon de chaque coté)
				if  position_SNP>=(debut+60) and position_SNP<=(fin-60) :
					position_on_baytes=60
					debut_bayte=int(position_SNP-position_on_baytes) 
					cas=1
					
					
				# Si je n'ai pas pu, je peux alors faire un bait avec le SNP décentré, avec le SNP vers le début
				elif position_SNP>(debut+10) and position_SNP<(fin-10) and (fin-debut)>=121 and position_SNP<(debut+60):
					debut_bayte=debut
					position_on_baytes=position_SNP-debut_bayte
					cas=2

				# Si je n'ai pas pu, je peux alors faire un bait avec le SNP décentré, avec le SNP vers la fin
				elif position_SNP>(debut+10) and position_SNP<(fin-10) and (fin-debut)>=121 and position_SNP>(fin-60):
					debut_bayte=fin-120
					position_on_baytes=position_SNP-debut_bayte
					cas=3
	
				
				# Maintenant que j'ai la position de mon baits et celle de mon SNP dans le bait, je peux l'imprimer
				bayte=sequence[(debut_bayte-1):(debut_bayte+120-1)]
				bayte1=list(bayte) ; bayte2=list(bayte)
				bayte1[position_on_baytes]=allele1
				bayte1="".join(bayte1)
				bayte2[position_on_baytes]=allele2 
				bayte2="".join(bayte2)
				toprint=">"+contig+"@"+str(position_SNP)+"|Allele_Dic2|pos_"+str(position_on_baytes+1)+"\n"+str(bayte1)+"\n"
				tmp.write(toprint)
				toprint=">"+contig+"@"+str(position_SNP)+"|Allele_Silur|pos_"+str(position_on_baytes+1)+"\n"+str(bayte2)+"\n"
				tmp.write(toprint)

	
	
	
	
	

print "\n\n      ---> done !    "+str(num_SNP_avec_annot)+"  SNPs ont des données d'annotation sur UTR ou exons !!"
print "\n\n"+str(num_bayte)+"  paires de Baits ont été formés !!"

toprint = "nbr de SNP avec annotation : "+str(num_SNP_avec_annot)
tmpfichierinfo.write(toprint+"\n")
toprint = "nbr de paires de baytes obtenues : "+str(num_bayte)
tmpfichierinfo.write(toprint+"\n")


#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------












































































