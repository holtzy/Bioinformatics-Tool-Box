#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON Get_Contigs_Infos.py
#     				Va résumer toutes les infos de divers fichiers pour chaque contig EPO.
#					Les fichiers d'entrées sont :
#						- Fichier fasta des contigs
#						- Fichier des SNPs détectés
#						- Blast des contigs (* 6)
#						- 
#						- autre...
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


parser = argparse.ArgumentParser(description= 'Va résumer toutes les infos de divers fichiers pour chaque contig EPO.')
parser.add_argument('-fasta', required=True, help='	séquences fasta des contigs d\'étude')
parser.add_argument('-fasta2', required=False, help='	séquences fasta des contigs de blé dur lorsque je met 2 fastas')
parser.add_argument('-psql', required=False, help=' fichier psql donnant les début et fins des CDS et UTR')
parser.add_argument('-SNP', required=True, help=' fichier de SNP')
parser.add_argument('-Orge_fasta', required=False, help=' fasta de l\'orge qui contient la position des contigs')
parser.add_argument('-Orge_liaison', required=False, help=' fichier de liaison contig blé / contig orge')
parser.add_argument('-Orge_blast', required=False, help=' fichier de blast contig EPO / contig orge')
parser.add_argument('-UTR5', required=False, help='	données de diversité génétique pour l\'UTR5')
parser.add_argument('-CDS', required=False, help='	données de diversité génétique pour l\'CDS')
parser.add_argument('-UTR3', required=False, help='	données de diversité génétique pour l\'UTR3')
parser.add_argument('-RPKM', required=False, help='	données concernant le RPKM')
parser.add_argument('-fic_3B', required=False, help='	position sur le 3B')
parser.add_argument('-out', required=True, help='	output')


args = parser.parse_args()
fic_fasta=args.fasta
fic_fasta2=args.fasta2
fic_snp=args.SNP
fic_psql=args.psql
orge_fasta=args.Orge_fasta
orge_liaison =args.Orge_liaison
orge_blast =args.Orge_blast
fic_UTR5=args.UTR5
fic_CDS=args.CDS
fic_UTR3=args.UTR3
out=args.out
fic_3B=args.fic_3B
fic_RPKM=args.RPKM




# Résumé du travail de ce script en étape :
	# Step 0 : Lecture des fichier fasta, récupération du nom de tous les contigs et de leurs longueurs
	# Step 1 : Récupération de la taille des UTR5 et UTR3, à partir du fichier .psql (blé dur) et du fasta avec les noms complets (blé tendre)
	# Step 2 : Intégration des infos relatives au SNP : nombre, position, FIS etc...
	# Step 3 : Données relatives aux positions des contigs sur l'orge
	# Step 4 : Ajout des données de DIversité génétiques
	# Step 5 : Utilisation des RPKMs
	# Step 6 : Ajout des positions des contigs sur le chromosome 3B
	# Step FINAL : Impression de tout ceci !






#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 0 : Collecte de l'ensemble des contigs dans un dictionnaire, j'en profite pour récupérer la longueur des contigs
dico_des_contigs=dict() ; tot=0

for line in open(fic_fasta):
	line=line.strip()
	if line.startswith(">"):
		tot=tot+1
		contig_name=line.replace(">","")

		#Je raccourci le nom si c'est un blé tendre, car sinon c'est inbuvable
		if contig_name.startswith("Traes"):
			contig_name=contig_name.split("|")[0]
		
		longueur_du_contig=0
		
	else:
		longueur_du_contig=longueur_du_contig + len(line)
		dico_des_contigs[contig_name]=str(longueur_du_contig)
		
#Fabrication de l'en tete
entete="contig_name"+"\t"+"longueur_en_pb"+"\t"
		
	
		
##### Meme travail si il y a un deuxième Fasta

if fic_fasta2 is not None :
	for line in open(fic_fasta2):
		line=line.strip()
		if line.startswith(">"):
			tot=tot+1
			contig_name=line.replace(">","")
	
			#Je raccourci le nom si c'est un blé tendre, car sinon c'est inbuvable
			if contig_name.startswith("Traes"):
				contig_name=contig_name.split("|")[0]
			
			longueur_du_contig=0
			
		else:
			longueur_du_contig=longueur_du_contig + len(line)
			dico_des_contigs[contig_name]=str(longueur_du_contig)
	
	
	
#Affichage final de l'étape		
print "\n-------\n\nTous les contigs de blé dur ont été répertorié, il y en a "+str(tot)+" \n\n------\n"
print "\n-------\n\nLa longueur des contigs a été calculé \n\n------\n"
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------











#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 1 : A l'aide de l'annotation, on va déterminer pour chaque contig, les bornes de ses CDS.
# Attention cette étape va etre faite en 2 fois : en effet : les contigs ont 2 origines possibles:
#	- blé dur --> dans ce cas on ira lire le fichier .psql
#	- blé tendre --> dans ce cas les infos sont dans l'entete des contigs.
# on va créer un grand dico qui va regrouper toutes ces infos pour les 2 types de contigs. Il contiendra : début CDS, fin CDS, orientation.

#---1/ Cas 1 : contig origine EPO :
if fic_psql is not None :
	dico_cds=dict() 
	
	#---1/ Contig en provenance du blé dur
	for line in open(fic_psql):
		line=line.strip()
		line=line.split(",")
		
		
		#Petite partie : au cas ou j'ai des , dans les noms, je dois recoller les morceaux a posteriori
		if len(line)==12:
			line[0]=line[0]+","+line[1]
			line[1]=line[2]+","+line[3]
			line.pop(2) ; line.pop(2)
	
		if len(line)==14:
			line[0]=line[0]+","+line[1]+","+line[2]
			line[1]=line[3]+","+line[4]+","+line[5]
			line.pop(2) ; line.pop(2) ; line.pop(2) ; line.pop(2)
	
		
		
		#Quel est le contig de la ligne?
		contig_name=line[0]
		
		#Quel est le début du CDS?
		if line[4] == "":
			start=line[5]
		else:
			start=line[4]
		
		#Quel est la fin?
		if line[8] == "":
			end=line[7]
		else:
			end=line[8]
	
	
		#Quel est l'orientation?
		if int(line[6])<0:
			orientation= "antisens"
		else:
			orientation= "sens"
			
		#Y a t'il un frameshift? Je n'enregistre le contig que si il n'y a pas de frameshift
		if int(line[6])==int(line[9]):
			dico_cds[contig_name]=str(start)+"\t"+str(end)+"\t"+str(orientation)
		else:
			dico_cds[contig_name]="-"+"\t"+"-"+"\t"+"-"
			
	



#---2/ Contig en provenance du blé tendre ---> Attention c'est vachement plus compliqué .....
dico_size=dict()

for line in open(fic_fasta) :
	line.strip()
	
	if line.startswith(">Traes"):
		line=line.split("|")
		
		#### - Je vire les cas chelou :
		#ou je n ai pas d'utr ni de CDS --> l'annotation du contig n'a pas été possible
		#Ou j'ai des UTR décomposé.. Comment c'est possible ?
		if len(line) < 7 or len(line[6].split(";"))>1 or len(line)>10 and len(line[11].split(";"))>1:
			continue
			
		### - Je récupère les principales infos: 
		#--Nom du contig
		contig=line[0].replace(">","")
		
		#--Taille de l'UTR5 : si j'ai juste les pipe sans valeurs, c'est que pas d'UTR5, donc taille=0
		if line[7] == "" and line[6] == "" :
			size_UTR5=0
		else:
			size_UTR5=abs(int(line[7])-int(line[6]))+1
		
		#--Taille de l'UTR3 : attention, si je n'ai pas d'UTR3, alors je n'ai meme pas les pipes à la fin...
		if len(line)<11 :
			size_UTR3=0
		else:
			size_UTR3=abs(int(line[10])-int(line[11]))+1


		#j'enregistre dans le dico
		dico_size[contig]=str(size_UTR5)+"\t"+str(size_UTR3)


#Maintenant que j'ai finit de faire mes 2 nouveaux dicos, je vais mettre toutes les infos dans mon dico bilan dico_des_contigs
for contig_name in dico_des_contigs:
	if contig_name in dico_size:
		dico_des_contigs[contig_name]=dico_des_contigs[contig_name]+"\t"+dico_size[contig_name]
	elif contig_name in dico_cds :
		dico_des_contigs[contig_name]=dico_des_contigs[contig_name]+"\t"+dico_cds[contig_name].split("\t")[0]+"\t"+dico_cds[contig_name].split("\t")[1]
	else:
		dico_des_contigs[contig_name]=dico_des_contigs[contig_name]+"\t"+"-"+"\t"+"-"


#Fabrication de l'en tete
entete=entete+"UTR5_size"+"\t"+"UTR3_size"+"\t"


print "\n-------\n\nLes infos relatives a la position des CDS ont été récupérées\n\n------\n"


#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------














#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 2 : Récupération des infos relatives au SNP,a partir du fichier de SNP. Je fais un nouveau Dico. A partir de ce nouveau dico je remplirai l'ancien!
memoire_contig=""
dico_info_snp=dict()

for line in open(fic_snp):

	#Récupération des infos de la ligne :
	line=line.split()
	he=line[3] ; FIS=line[4] ; nbr_indiv_genotypes=line[5] ; 
	contig_name=line[6].replace(">","") 
	if contig_name.startswith("Trae"):
		contig_name=contig_name.split("|")[0]
	snp_position=line[7] 


	###Si je commence un nouveau contig
	if memoire_contig=="" or memoire_contig!=contig_name:
		
		##Si ce n'est pas le premier, 
		if memoire_contig != "":
		
			#je calcule mes résultats pour ce contig,
			if nb_de_he > 0:
				he_mean=he_tot/nbr_snp
				FIS_mean=FIS_tot/nbr_snp
				pourcentage_FIS_faible=float(nbr_snp_FIS_faible)/float(nbr_snp)*100 
				pourcentage_FIS_faible=round(pourcentage_FIS_faible,1)
			else:
				he_mean="tri-allelic"
				FIS_mean="tri-allelic"
			
			
			nbr_indiv_genotypes_mean=nbr_indiv_genotypes_tot/nbr_snp
			
			#a chaque fois j'ajoute au dictionnaire
			dico_info_snp[memoire_contig]=str(nbr_snp)
			dico_info_snp[memoire_contig]=dico_info_snp[memoire_contig]+"\t"+str(snp_position_tot)
			dico_info_snp[memoire_contig]=dico_info_snp[memoire_contig]+"\t"+str(he_mean)
			dico_info_snp[memoire_contig]=dico_info_snp[memoire_contig]+"\t"+str(FIS_mean)
			dico_info_snp[memoire_contig]=dico_info_snp[memoire_contig]+"\t"+str(pourcentage_FIS_faible)
			dico_info_snp[memoire_contig]=dico_info_snp[memoire_contig]+"\t"+str(nbr_indiv_genotypes_mean)
		
		
		##Et dans tous les cas je remets mes compteurs à 0
		nbr_snp=0 ; he_tot=0 ; nb_de_he=0 ; FIS_tot=0 ; nbr_indiv_genotypes_tot=0 ; snp_position_tot=[] ; nbr_snp_FIS_faible=0
		
		##Et je remet ma mémoire à jour pour ce nouveau contig
		memoire_contig=contig_name


	###Dans tous les cas j'incrémente mes compteurs
	nbr_snp += 1
	if he != "-":
		nb_de_he+=1
		he_tot += float(he)
		FIS_tot+= float(FIS)
		if float(FIS) < 0.3:
			nbr_snp_FIS_faible+=1
		
	nbr_indiv_genotypes_tot += int(nbr_indiv_genotypes)
	snp_position_tot.append(snp_position)
	



#Maintenant que j'ai finit de faire mon nouveau dico, je vais mettre toutes les infos dans mon dico bilan dico_des_contigs
for contig_name in dico_des_contigs:
	if contig_name in dico_info_snp:
		#Juste, pour pas que les position snp soient réparties dans plusieurs colonnes, mais bien dans une seule
		dico_info_snp[contig_name]=dico_info_snp[contig_name].replace(", ","@@")
		
		dico_des_contigs[contig_name]=dico_des_contigs[contig_name]+"\t"+dico_info_snp[contig_name]
	else:
		dico_des_contigs[contig_name]=dico_des_contigs[contig_name]+"\t"+"0"+"\t"+"-"+"\t"+"-"+"\t"+"-"+"\t"+"-"+"\t"+"-"


#Fabrication de l'en tete
entete=entete+"nbr_snp"+"\t"+"snp_positions"+"\t"+"he_moyen"+"\t"+"FIS_moyen"+"\t"+"%_de_FIS_faible"+"\t"+"nbr_moyen_indiv_par_snp"+"\t"


print "\n-------\n\nLes infos relatives aux SNPs ont été récupérées\n\n------\n"
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
















#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 3 : Mise en relation avec les données de l'orge !
# il me faut le fasta de l'orge qui contient les positions des contigs
# Et le fichier d'ensembl qui fait le lien entre contig de blé tendre / contig d'orge :
# Et éventuellement un fichier de blast qui fait le lien entre contig de blé dur et contig d'orge.

if orge_fasta is not None :

	# -- Je commence par faire un dico des contigs d'orges et de leur position.
	dico_orge=dict()
	for line in open(orge_fasta):
		line=line.strip()
		if line.startswith(">") and re.search("morex",line) is None:
			line=line.split()
			contig_name_orge=line[0].replace(">","")
			contig_name_orge=contig_name_orge.split(".")[0]
			chromosome=line[2].split(":")[2]
			blast_start=line[2].split(":")[3]
			blast_end=line[2].split(":")[4]
			dico_orge[contig_name_orge]=str(chromosome)+"\t"+str(blast_start)

	
	# -- Ensuite j'intègre la liaison entre orge et blé tendre si elle est dispo
	if orge_liaison is not None :
		dico_liaison=dict()
		for line in open(orge_liaison):
			line=line.strip()
			line=line.split()
			if len(line)>=3:
				dico_liaison[line[0]] = line[2]


 	# -- Ensuite j'intégre la liaison entre orge et blé dur EPO si elle est dispo
 	if orge_blast is not None :
 		for line in open(orge_blast):
			line=line.strip()
			line=line.split()
			gene_orge=line[1].split(".")[0]
			dico_liaison[line[0]] = gene_orge


	#Maintenant que j'ai finit de faire mon nouveau dico, je vais mettre toutes les infos dans mon dico bilan dico_des_contigs
	for contig_name in dico_des_contigs:
		if contig_name in dico_liaison :
			orge_corres=dico_liaison[contig_name]
			if orge_corres in dico_orge :
				dico_des_contigs[contig_name]=dico_des_contigs[contig_name]+"\t"+dico_orge[orge_corres]
			else:
				dico_des_contigs[contig_name]=dico_des_contigs[contig_name]+"\t"+"-"+"\t"+"-"
		else:
			dico_des_contigs[contig_name]=dico_des_contigs[contig_name]+"\t"+"-"+"\t"+"-"
	
	

	#Fabrication de l'en tete
	entete=entete+"chromo_hordeum"+"\t"+"position_hordeum"





	print "\n-------\n\nLes infos concernant les positions sur l'orge ont été récupérés\n\n------\n"
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------









	
	


#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 4 : Ajout des données de diversité calculées par le script de Benoit Nabholz : 3 fichiers d'entré : UTR5, CDS, UTR3
def add_diversity_info(fic):
	
	dico_diversite=dict()
	for line in open(fic):
		line=line.split()
	
		#Récupération du nom du contig
		contig_name=line[0].replace("_Contig","|Contig")
		contig_name=contig_name.replace("_less","|less")
		contig_name=contig_name.replace("_orig","|orig")	
		contig_name=contig_name.replace("_like","|like")
		contig_name=contig_name.replace("_compl","|compl")
	
		#Récupération des variables calculées par egglib
		size=line[1]
		S=line[2]
		Pi=line[3]
		W=line[4]
		D=line[5]
		Ts=line[6]
		PS=line[7]
		PN=line[8]
		NSS=line[9]
		GC3=line[10]
	
	
		#Ajout au dictionnaire
		dico_diversite[contig_name]=str(size)+"\t"+str(S)+"\t"+str(Pi)+"\t"+str(W)+"\t"+str(D)+"\t"+str(Ts)+"\t"+str(PS)+"\t"+str(PN)+"\t"+str(NSS)+"\t"+str(GC3)
	

	#Et liaison avec le dico total ! :
	for contig_name in dico_des_contigs:
		if contig_name in dico_diversite:
			dico_des_contigs[contig_name]=dico_des_contigs[contig_name]+"\t"+str(dico_diversite[contig_name])
		else:
			dico_des_contigs[contig_name]=dico_des_contigs[contig_name]+"\t"+"-"+"\t"+"-"+"\t"+"-"+"\t"+"-"+"\t"+"-"+"\t"+"-"+"\t"+"-"+"\t"+"-"+"\t"+"-"+"\t"+"-"


if fic_UTR5 is not None and fic_CDS is not None and fic_UTR3 is not None:
	add_diversity_info(fic_UTR5)
	add_diversity_info(fic_CDS)
	add_diversity_info(fic_UTR3)
	
	#Fabrication de l'en tete
	entete=entete+"UTR5_size"+"\t"+"UTR5_S"+"\t"+"UTR5_Pi"+"\t"+"UTR5_W"+"\t"+"UTR5_D"+"\t"+"UTR5_Ts"+"\t"+"UTR5_PS"+"\t"+"UTR5_PN"+"\t"+"UTR5_NSS"+"\t"+"UTR5_GC3"+"\t"
	entete=entete+"CDS_size"+"\t"+"CDS_S"+"\t"+"CDS_Pi"+"\t"+"CDS_W"+"\t"+"CDS_D"+"\t"+"CDS_Ts"+"\t"+"CDS_PS"+"\t"+"CDS_PN"+"\t"+"CDS_NSS"+"\t"+"CDS_GC3"+"\t"
	entete=entete+"UTR3_size"+"\t"+"UTR3_S"+"\t"+"UTR3_Pi"+"\t"+"UTR3_W"+"\t"+"UTR3_D"+"\t"+"UTR3_Ts"+"\t"+"UTR3_PS"+"\t"+"UTR3_PN"+"\t"+"UTR3_NSS"+"\t"+"UTR3_GC3"+"\t"
	
	
	print "\n-------\n\nLes infos de diversite génétique ont été récupérées\n\n------\n"
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------








#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 5 : Ajout des données de RPKM
if fic_RPKM is not None :
	dico_RPKM=dict()
	line_number=0
	
	
	for line in open(fic_RPKM):
		line=line.split()
		nb_indiv=len(line)-1
		line_number+=1
		
		if line_number>2:
			#Récupération du nom du contig
			contig_name=line[0]
			if contig_name.startswith("Trae"):
				contig_name=contig_name.split("|")[0]
			
			#ajout des vleurs au dico rpkm, et Calcul de la moyenne des RPKM
			tot=0 ; _nb_=0
			for i in range(1,nb_indiv):
				#dico_RPKM[contig_name]=dico_RPKM[contig_name]+"\t"+str(line[i])
				if float(line[i]) != 0:
					_nb_+=1
					tot=tot+float(line[i])	
			
			if _nb_ != 0:
				mean_RPKM=tot/float(_nb_)	
				#remplissage du dico	
				dico_RPKM[contig_name]=str(mean_RPKM)	
	
	#Ajout des données au dico général
	for contig_name in dico_des_contigs:
		if contig_name in dico_RPKM:
			dico_des_contigs[contig_name]=dico_des_contigs[contig_name]+"\t"+str(dico_RPKM[contig_name])
		else:
			dico_des_contigs[contig_name]=dico_des_contigs[contig_name]+"\t"+"-"
	
	#Fabrication de l'en tete
	entete=entete+"\t"+"RPKM_mean"
	
			
	
	print "\n-------\n\nLes infos de RPKM ont été récupérées\n\n------\n"
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------









#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 6 : Ajout des données de position sur le 3B
if fic_3B is not None :
	dico_3B=dict()
	
	for line in open(fic_3B):
		line=line.strip()
		line=line.split()
		contig=line[0]
		dico_3B[contig]= line[2]
	
		
	#Ajout des données au dico général
	for contig_name in dico_des_contigs:
		if contig_name in dico_3B:
			dico_des_contigs[contig_name]=dico_des_contigs[contig_name]+"\t"+str(dico_3B[contig_name])
		else:
			dico_des_contigs[contig_name]=dico_des_contigs[contig_name]+"\t"+"-"
	
	
	#Fabrication de l'en tete
	entete=entete+"\t"+"position_3B"

	
	

	
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
	
	
	








#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP FINAL : Je transforme mon dico en fichier de sortie!

#Enregistrement
tmp=open(out,"w")
tmp.write(entete+"\n")
for contig_name in dico_des_contigs:
	tmp.write(contig_name+"\t"+dico_des_contigs[contig_name]+"\n")
tmp.close()

print "\n-------\n\nLe fichier de sortie --- "+out+" --- est prêt!\n\n------\n"
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------





















