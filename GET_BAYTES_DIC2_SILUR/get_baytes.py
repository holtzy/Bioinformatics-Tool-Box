#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON get_baytes
#     				Ce script permet de récupérer des baytes à partir de :
#						-un fichier de SNP
#						-un fichier donnant l'emplacement des introns dans les contigs
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


parser = argparse.ArgumentParser(description= 'get the baytes from durum wheat SNPs')
parser.add_argument('-fasta', required=True, help='	fasta de EPO')
parser.add_argument('-position_intron', required=False, help='	position des introns dans les contigs')
parser.add_argument('-SNP', required=True, help='	fichier de SNP')
parser.add_argument('-UTR', required=False, help=' fichier de SNP avec l info UTR ou non')
parser.add_argument('-out', required=True, help=' fichier de sortie')


args = parser.parse_args()
fasta=args.fasta
position_intron=args.position_intron
SNP=args.SNP
UTR=args.UTR
out=args.out




tmpfichierinfo=open("fichier_info_baytes.txt","w")




#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 0 : Collecte de l'ensemble des contigs dans un dictionnaire
print "\n\nSTEP 0 : dictionnaire des contigs"
dico_EPO=dict() 
nombre=0

for line in open(fasta):
	line=line.strip()
	if line.startswith(">"):
		nombre=nombre+1
		contig=line.replace(">","")
		
		#Attention si c'est un contig de blé tendre je dois raccourcir le nom
		if contig.startswith("Traes"):
			contig=contig.split("|")[0]
		
		dico_EPO[contig]=""
		
	else:
		dico_EPO[contig]=dico_EPO[contig]+line


print "\n\n      ---> done !     "+str(nombre)+"     Contigs ont été récupérés dans le fichier fasta"
toprint= "nbr de contigs récupérés dans le fichier Fasta : "+str(nombre)
tmpfichierinfo.write(toprint+"\n")

#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------







#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 1 : Collecte des positions des introns dan un dictionnaire dans le cas EPO
dico_intron=dict()

if position_intron is not None:
	print "\n\nSTEP 1 : dictionnaire des positions des introns"
	dico_blast=dict()
	nombre=0 ; nbr_uniq=0
	
	
	#A chaque nouvelle ligne je récupère toutes les infos : contig, scores de blast, et position des introns
	for line in open(position_intron):
		nombre=nombre+1
		line=line.strip()
		line=line.split()
		contig=line[1]
		contig=contig.replace("_Con","|Con")
		contig=contig.replace("_like","|like")
		contig=contig.replace("_compl","|compl")
		contig=contig.replace("_orig","|orig")
		contig=contig.replace("_less","|less")
		identity=line[2]
		coverage=line[3]
		positions=line[4]
		
		#Je stocke les valeurs de blast dans un dico à part
		dico_blast[contig]=identity+","+coverage
		
		#Si on recontre ce contig pour la premiere fois, on récupère la position de ses introns.
		if contig not in dico_intron:
			nbr_uniq+=1
			dico_intron[contig]=positions
		
		#Si le contig a déja été repertorié, il faut alors choisir son MEILLEUR blast et les positions d'introns correspondantes.
		else:
		
			#Si le nouveau contig a un meilleur % d'identiter, alors on le garde d'office, il remplace alors l'ancien dans les dicos
			if identity > dico_blast[contig].split(",")[0]:
				dico_intron[contig]=positions
				dico_blast[contig]=identity+","+coverage
				
				
				
			#Si c est pile égal, alors on regarde la distance de couverture.
			if identity == dico_blast[contig].split(",")[0]:
				if coverage > dico_blast[contig].split(",")[1]:
					dico_intron[contig]=positions
					dico_blast[contig]=identity+","+coverage
									
	
	
	print "\n\n      ---> done  ! nbr de lignes dans le fichier position_intron : "+str(nombre)+"  cas" 
	toprint = "nbr de lignes dans le fichier position des introns : "+str(nombre)
	tmpfichierinfo.write(toprint+"\n")
	toprint = "nbr de contig pour lesquels on a au final l'info des introns : "+str(nbr_uniq)
	tmpfichierinfo.write(toprint+"\n")


#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------












#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 1 : Collecte des positions des introns dan un dictionnaire
print "\n\nSTEP 1 : dictionnaire des positions des introns"




dico_UTR=dict()
nombre=0 
affichage_message="no"


#A chaque nouvelle ligne je récupère toutes les infos : contig, positions des exons
for line in open(fasta):
	line.strip()
	
	if line.startswith(">Traes"):
		affichage_message="yes"
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
					
	
	
			dico_UTR[contig]="1..."+str(fin_UTR5)+","+str(fin_des_exons)+"..."+str(fin_UTR3)





		#-----------reverse complement
		if strand == "-1":
	
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
			dico_UTR[contig]="1..."+str(fin_UTR3)+","+str(fin_des_exons)+"..."+str(fin_UTR5)



if affichage_message=="yes" :	
	print "\n\n      ---> done  ! nbr de lignes dans le fichier position_intron : "+str(nombre)+"  cas" 
	toprint = "nbr de lignes dans le fichier position des introns : "+str(nombre)
	tmpfichierinfo.write(toprint+"\n")
	toprint = "nbr de contig pour lesquels on a au final l'info des introns : "+str(len(dico_intron))
	tmpfichierinfo.write(toprint+"\n")


#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------








#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 2 : Collecte des infos des SNPs
print "\n\nSTEP 2 : dictionnaire des positions des SNPs"
dico_SNP=dict() 
nombre=0 ; nombre_filtre=0

for line in open(SNP):
	nombre=nombre+1
	line=line.strip()
	line=line.split()
	He=float(line[3])
	nb_indiv=float(line[5])
	nb_indiv_potentiel= float( (len(line)-8)/2 )
	contig=line[6].replace(">","")
	
	#Je raccourci le nom pour le blé tendre
	if contig.startswith("Trae"):
		contig=contig.split("|")[0]
		
	position=line[7] 
	SNP_name=contig+"@"+position
	
	#Je sélectionne uniquement les SNPs ayant plus de x% d'individus génotypés
	if nb_indiv > (nb_indiv_potentiel * 0.14) and He > 0.2:
		nombre_filtre+=1
		nucl1=line[9].split("|")[0][0]
		nucl2=line[11].split("|")[0][0]
	
		dico_SNP[SNP_name]=str(He)+","+str(nb_indiv)+","+str(position)+","+nucl1+","+nucl2		
			
			


print "\n\n      ---> done !    "+str(nombre)+"   SNPs ont été testés !!"
print "           Parmis eux    "+str(nombre_filtre)+"   SNPs ont passés les filtres de nbr d'individus genotypes et de He !!"
toprint = "nbr SNPs potentiellement utilisables pour faire des baytes : "+str(nombre)
tmpfichierinfo.write(toprint+"\n")


#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------














#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 3 : Au final quels SNPs peut on utiliser pour fabriquer les baytes???
print "\n\nSTEP 3 : récupération des baytes !! "
tmp=open(out,"w")
num_bayte=0
inventaire_des_baytes=dict()


#Je procède SNP par SNP, pour voir a chaque fois si il est utilisable ou non. J'ai plein de critère à chécker :
for SNP_name in dico_SNP:
	position_SNP=int(dico_SNP[SNP_name].split(",")[2])
	contig=SNP_name.split("@")[0]
	
	#Du coup il faut que j'aille récupérer la séquence du contig
	sequence=dico_EPO[contig]
	allele1=dico_SNP[SNP_name].split(",")[3]
	allele2=dico_SNP[SNP_name].split(",")[4]
	
	# *** SI je connais la position des introns pour le contig en question
	if contig in dico_intron:

		#Je vais voir si le SNP peut faire l'affaire exon par exon
		position=dico_intron[contig].split(",")
		for exon in position:

			fin=int(exon.split("..")[1])
			debut=int(exon.split("..")[0])
			

			#J'ai récupéré toutes les infos, je peux poser mes conditions pour voir si mon SNP permet de faire un baytes!!!
			#Petite sortie pour voir les caractéristiques du bayte
			#Condition : au moins 10 nucléotides avant et après le SNP sur l'exon, et une taille d'exon supérieur ou égale à 200pb
			#Si je peux faire un baytes de 200 pb avec le SNP au milieu je le fais.
			#Sinon j'en fait de 120 pb
			
			
			# *** SI je peux faire un bayte de 200pb avec le SNP au milieu
			if (fin-debut)>=201 and position_SNP>=(debut+100) and position_SNP<=(fin-100):
				num_bayte+=1
				inventaire_des_baytes[SNP_name]="OK"
				debut_bayte=int(position_SNP-100) 
				bayte=sequence[(debut_bayte-1):(debut_bayte+200-1)]
				bayte1=list(bayte) ; bayte2=list(bayte)
				bayte1[101-1]=allele1
				bayte1="".join(bayte1)
				bayte2[101-1]=allele2 
				bayte2="".join(bayte2)
				toprint=">"+contig+"@"+str(position_SNP)+"|Allele_Dic2|pos_101"+"\n"+str(bayte1)+"\n"
				tmp.write(toprint)
				toprint=">"+contig+"@"+str(position_SNP)+"|Allele_Silur|pos_101"+"\n"+str(bayte2)+"\n"
				tmp.write(toprint)


			# *** SI je peux faire un bayte de 200pb avec le SNP décalé vers le début
			elif position_SNP>(debut+10) and position_SNP<(fin-10) and (fin-debut)>=201 and position_SNP<(debut+100):
				num_bayte+=1
				inventaire_des_baytes[SNP_name]="OK"
				debut_bayte=debut
				position_on_baytes=position_SNP-debut_bayte+1
				bayte1=sequence[(debut_bayte-1):debut_bayte+200-1]
				bayte2=sequence[(debut_bayte-1):debut_bayte+200-1]
				bayte1=list(bayte1) ; bayte2=list(bayte2)
				bayte1[position_on_baytes-1]=allele1
				bayte2[position_on_baytes-1]=allele2
				bayte1="".join(bayte1)
				bayte2="".join(bayte2)
				toprint=">"+contig+"@"+str(position_SNP)+"|Allele_Dic2|pos_"+str(position_on_baytes)+"\n"+bayte1+"\n"
				tmp.write(toprint)
				toprint=">"+contig+"@"+str(position_SNP)+"|Allele_Silur|pos_"+str(position_on_baytes)+"\n"+bayte2+"\n"
				tmp.write(toprint)


			# *** SI je peux faire un bayte de 200pb avec le SNP décalé vers la fin
			elif position_SNP>(debut+10) and position_SNP<(fin-10) and (fin-debut)>=201 and position_SNP>(fin-100):
				num_bayte+=1
				inventaire_des_baytes[SNP_name]="OK"
				debut_bayte=fin-199
				position_on_baytes=position_SNP-debut_bayte+1
				bayte1=sequence[(debut_bayte-1):debut_bayte+200-1]
				bayte2=sequence[(debut_bayte-1):debut_bayte+200-1]
				bayte1=list(bayte1) ; bayte2=list(bayte2)
				bayte1[position_on_baytes-1]=allele1
				bayte2[position_on_baytes-1]=allele2
				bayte1="".join(bayte1)
				bayte2="".join(bayte2)
				toprint=">"+contig+"@"+str(position_SNP)+"|Allele_Dic2|pos_"+str(position_on_baytes)+"\n"+bayte1+"\n"
				tmp.write(toprint)
				toprint=">"+contig+"@"+str(position_SNP)+"|Allele_Silur|pos_"+str(position_on_baytes)+"\n"+bayte2+"\n"
				tmp.write(toprint)


			# *** SI je ne peux pas faire un bayte de 200pb avec le SNP au milieu, mais que je peux en faire un de 120 avec le SNP au milieu
			elif position_SNP>(debut+10) and position_SNP<(fin-10) and (fin-debut)>=121 and position_SNP>=(debut+60) and position_SNP<=(fin-60) :
				num_bayte+=1
				inventaire_des_baytes[SNP_name]="OK"
				debut_bayte=int(position_SNP-60) 
				bayte=sequence[(debut_bayte-1):(debut_bayte+120-1)]
				bayte1=list(bayte) ; bayte2=list(bayte)
				bayte1[61-1]=allele1
				bayte1="".join(bayte1)
				bayte2[61-1]=allele2 
				bayte2="".join(bayte2)
				toprint=">"+contig+"@"+str(position_SNP)+"|Allele_Dic2|pos_61"+"\n"+str(bayte1)+"\n"
				tmp.write(toprint)
				toprint=">"+contig+"@"+str(position_SNP)+"|Allele_Silur|pos_61"+"\n"+str(bayte2)+"\n"
				tmp.write(toprint)

	
			# *** SI je ne peux pas faire un bayte de 200pb avec le SNP au milieu, mais que je peux en faire un de 120 avec le SNP vers le début
			elif position_SNP>(debut+10) and position_SNP<(fin-10) and (fin-debut)>=121 and position_SNP<(debut+60):
				num_bayte+=1
				inventaire_des_baytes[SNP_name]="OK"
				debut_bayte=debut
				position_on_baytes=position_SNP-debut_bayte+1
				bayte1=sequence[(debut_bayte-1):debut_bayte+120-1]
				bayte2=sequence[(debut_bayte-1):debut_bayte+120-1]
				bayte1=list(bayte1) ; bayte2=list(bayte2)
				bayte1[position_on_baytes-1]=allele1
				bayte2[position_on_baytes-1]=allele2
				bayte1="".join(bayte1)
				bayte2="".join(bayte2)
				toprint=">"+contig+"@"+str(position_SNP)+"|Allele_Dic2|pos_"+str(position_on_baytes)+"\n"+bayte1+"\n"
				tmp.write(toprint)
				toprint=">"+contig+"@"+str(position_SNP)+"|Allele_Silur|pos_"+str(position_on_baytes)+"\n"+bayte2+"\n"
				tmp.write(toprint)


			# *** SI je ne peux pas faire un bayte de 200pb avec le SNP au milieu, mais que je peux en faire un de 120 avec le SNP vers la fin
			elif position_SNP>(debut+10) and position_SNP<(fin-10) and (fin-debut)>=121 and position_SNP>(fin-60):
				num_bayte+=1
				inventaire_des_baytes[SNP_name]="OK"
				debut_bayte=fin-119
				position_on_baytes=position_SNP-debut_bayte+1
				bayte1=sequence[(debut_bayte-1):debut_bayte+120-1]
				bayte2=sequence[(debut_bayte-1):debut_bayte+120-1]
				bayte1=list(bayte1) ; bayte2=list(bayte2)
				bayte1[position_on_baytes-1]=allele1
				bayte2[position_on_baytes-1]=allele2
				bayte1="".join(bayte1)
				bayte2="".join(bayte2)
				toprint=">"+contig+"@"+str(position_SNP)+"|Allele_Dic2|pos_"+str(position_on_baytes)+"\n"+bayte1+"\n"
				tmp.write(toprint)
				toprint=">"+contig+"@"+str(position_SNP)+"|Allele_Silur|pos_"+str(position_on_baytes)+"\n"+bayte2+"\n"
				tmp.write(toprint)
					
						
						
						

print "\n\n      ---> done !    "+str(num_bayte)+"  SNPs ont donnés naissance a des baytes !!"
toprint = "nbr de paires de baytes obtenues : "+str(num_bayte)
tmpfichierinfo.write(toprint+"\n")

#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------


















###################################
### PARTIE 2 : BAYTES FROM UTR ####
####################################

#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------

#STEP1 : dico des positions des UTR
if UTR is not None:
	num=0
	
	for line in open(UTR) :
		num+=1
		
		#Pour pas prendre en compte la ligne d'en tete
		if num>1:
			line.strip()
			line=line.split()
			
					
			contig=line[0]
			contig=contig.split("@")[0]
			deb_CDS=line[2]
			fin_CDS=line[3]
			fin_contig=line[1]
		
			dico_UTR[contig]="1..."+str(deb_CDS)+","+str(fin_CDS)+"..."+str(fin_contig)








#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------

#STEP2 : récupération des baytes from UTR
print "\n\nSTEP 4 : récupération des baytes from UTR !! "
nbr_SNP_in_UTR_and_exon=0


#Je procède SNP par SNP, pour voir a chaque fois si il est utilisable ou non.
for SNP_name in dico_SNP:
	position_SNP=int(dico_SNP[SNP_name].split(",")[2])
	contig=SNP_name.split("@")[0]
	
	#Du coup il faut que j'aille récupérer mon baytes dans la séquence du contig
	sequence=dico_EPO[contig]
	allele1=dico_SNP[SNP_name].split(",")[3]
	allele2=dico_SNP[SNP_name].split(",")[4]

	
	if contig in dico_UTR:

		#Je vais voir si le SNP peut faire l'affaire UTR par UTR
		pos = dico_UTR[contig].split(",")
		for i in range(1,3) :
			UTR=pos[i-1]
			fin=UTR.split("...")[1]
			debut=UTR.split("...")[0]

			#Juste pour bien avoir l'info des UTR
			if fin != "-" and debut != "-" :
				fin=int(fin)
				debut=int(debut)

				#Et je vérifie que le SNP n a pas déja été trouvé dans un exon 
				if SNP_name not in inventaire_des_baytes :


					# *** SI je peux faire un bayte de 200pb avec le SNP au milieu
					if (fin-debut)>=201 and position_SNP>=(debut+100) and position_SNP<=(fin-100):
						num_bayte+=1
						inventaire_des_baytes[SNP_name]="OK"
						debut_bayte=int(position_SNP-100) 
						bayte=sequence[(debut_bayte-1):(debut_bayte+200-1)]
						bayte1=list(bayte) ; bayte2=list(bayte)
						bayte1[101-1]=allele1
						bayte1="".join(bayte1)
						bayte2[101-1]=allele2 
						bayte2="".join(bayte2)
						toprint=">"+contig+"@"+str(position_SNP)+"|Allele_Dic2|pos_101"+"\n"+str(bayte1)+"\n"
						tmp.write(toprint)
						toprint=">"+contig+"@"+str(position_SNP)+"|Allele_Silur|pos_101"+"\n"+str(bayte2)+"\n"
						tmp.write(toprint)
		
		
					# *** SI je peux faire un bayte de 200pb avec le SNP décalé vers le début
					elif position_SNP>(debut+10) and position_SNP<(fin-10) and (fin-debut)>=201 and position_SNP<(debut+100):
						num_bayte+=1
						debut_bayte=debut
						position_on_baytes=position_SNP-debut_bayte+1
						bayte1=sequence[(debut_bayte-1):debut_bayte+200-1]
						bayte2=sequence[(debut_bayte-1):debut_bayte+200-1]
						bayte1=list(bayte1) ; bayte2=list(bayte2)
						bayte1[position_on_baytes-1]=allele1
						bayte2[position_on_baytes-1]=allele2
						bayte1="".join(bayte1)
						bayte2="".join(bayte2)
						toprint=">"+contig+"@"+str(position_SNP)+"|Allele_Dic2|pos_"+str(position_on_baytes)+"\n"+bayte1+"\n"
						tmp.write(toprint)
						toprint=">"+contig+"@"+str(position_SNP)+"|Allele_Silur|pos_"+str(position_on_baytes)+"\n"+bayte2+"\n"
						tmp.write(toprint)
		
		
					# *** SI je peux faire un bayte de 200pb avec le SNP décalé vers la fin
					elif position_SNP>(debut+10) and position_SNP<(fin-10) and (fin-debut)>=201 and position_SNP>(fin-100):
						num_bayte+=1
						debut_bayte=fin-199
						position_on_baytes=position_SNP-debut_bayte+1
						bayte1=sequence[(debut_bayte-1):debut_bayte+200-1]
						bayte2=sequence[(debut_bayte-1):debut_bayte+200-1]
						bayte1=list(bayte1) ; bayte2=list(bayte2)
						bayte1[position_on_baytes-1]=allele1
						bayte2[position_on_baytes-1]=allele2
						bayte1="".join(bayte1)
						bayte2="".join(bayte2)
						toprint=">"+contig+"@"+str(position_SNP)+"|Allele_Dic2|pos_"+str(position_on_baytes)+"\n"+bayte1+"\n"
						tmp.write(toprint)
						toprint=">"+contig+"@"+str(position_SNP)+"|Allele_Silur|pos_"+str(position_on_baytes)+"\n"+bayte2+"\n"
						tmp.write(toprint)
		
		
					# *** SI je ne peux pas faire un bayte de 200pb avec le SNP au milieu, mais que je peux en faire un de 120 avec le SNP au milieu
					elif position_SNP>(debut+10) and position_SNP<(fin-10) and (fin-debut)>=121 and position_SNP>=(debut+60) and position_SNP<=(fin-60) :
						num_bayte+=1
						debut_bayte=int(position_SNP-60) 
						bayte=sequence[(debut_bayte-1):(debut_bayte+120-1)]
						bayte1=list(bayte) ; bayte2=list(bayte)
						bayte1[61-1]=allele1
						bayte1="".join(bayte1)
						bayte2[61-1]=allele2 
						bayte2="".join(bayte2)
						toprint=">"+contig+"@"+str(position_SNP)+"|Allele_Dic2|pos_61"+"\n"+str(bayte1)+"\n"
						tmp.write(toprint)
						toprint=">"+contig+"@"+str(position_SNP)+"|Allele_Silur|pos_61"+"\n"+str(bayte2)+"\n"
						tmp.write(toprint)
		
			
					# *** SI je ne peux pas faire un bayte de 200pb avec le SNP au milieu, mais que je peux en faire un de 120 avec le SNP vers le début
					elif position_SNP>(debut+10) and position_SNP<(fin-10) and (fin-debut)>=121 and position_SNP<(debut+60):
						num_bayte+=1
						debut_bayte=debut
						position_on_baytes=position_SNP-debut_bayte+1
						bayte1=sequence[(debut_bayte-1):debut_bayte+120-1]
						bayte2=sequence[(debut_bayte-1):debut_bayte+120-1]
						bayte1=list(bayte1) ; bayte2=list(bayte2)
						bayte1[position_on_baytes-1]=allele1
						bayte2[position_on_baytes-1]=allele2
						bayte1="".join(bayte1)
						bayte2="".join(bayte2)
						toprint=">"+contig+"@"+str(position_SNP)+"|Allele_Dic2|pos_"+str(position_on_baytes)+"\n"+bayte1+"\n"
						tmp.write(toprint)
						toprint=">"+contig+"@"+str(position_SNP)+"|Allele_Silur|pos_"+str(position_on_baytes)+"\n"+bayte2+"\n"
						tmp.write(toprint)
		
		
					# *** SI je ne peux pas faire un bayte de 200pb avec le SNP au milieu, mais que je peux en faire un de 120 avec le SNP vers la fin
					elif position_SNP>(debut+10) and position_SNP<(fin-10) and (fin-debut)>=121 and position_SNP>(fin-60):
						num_bayte+=1
						debut_bayte=fin-119
						position_on_baytes=position_SNP-debut_bayte+1
						bayte1=sequence[(debut_bayte-1):debut_bayte+120-1]
						bayte2=sequence[(debut_bayte-1):debut_bayte+120-1]
						bayte1=list(bayte1) ; bayte2=list(bayte2)
						bayte1[position_on_baytes-1]=allele1
						bayte2[position_on_baytes-1]=allele2
						bayte1="".join(bayte1)
						bayte2="".join(bayte2)
						toprint=">"+contig+"@"+str(position_SNP)+"|Allele_Dic2|pos_"+str(position_on_baytes)+"\n"+bayte1+"\n"
						tmp.write(toprint)
						toprint=">"+contig+"@"+str(position_SNP)+"|Allele_Silur|pos_"+str(position_on_baytes)+"\n"+bayte2+"\n"
						tmp.write(toprint)
					
						
							
						


print "\n\n      ---> done !    "+str(num_bayte)+"  SNPs ont donné naissance a des baytes maintenant que ceux des UTR ont été ajoutés!!"


toprint = "nbr de paires de baytes obtenues avec l'ajout des baytes des UTR : "+str(num_bayte)
tmpfichierinfo.write(toprint+"\n")
toprint = "nbr de SNP présent dans un UTR ET dans un exon : "+str(nbr_SNP_in_UTR_and_exon)
tmpfichierinfo.write(toprint+"\n")



#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------




































































































