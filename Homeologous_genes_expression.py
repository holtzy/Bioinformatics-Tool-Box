#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON Homeologous Genes expression
#     				Ce script permet de calculer l'expression des sites divergent entre homéogénomes
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


parser = argparse.ArgumentParser(description= 'permet de calculer l\'expression des sites divergent entre homéogénomes')
parser.add_argument('-fasta', required=True, help='	fasta contenant les séquences des homéogènes')
parser.add_argument('-blast', required=True, help='	résultat de l\'autoblast de ces séquences')
parser.add_argument('-homeogenes', required=True, help=' liste des couples d\'homéocontigs')
parser.add_argument('-alr', required=True, help=' fichier .alr de résultat du mapping')
parser.add_argument('-nbr_reads', required=True, help=' nbr de reads total par indiv')
parser.add_argument('-out', required=True, help=' fichier de sortie')


args = parser.parse_args()
fasta=args.fasta
blast=args.blast
homeogenes=args.homeogenes
alr=args.alr
output=args.out
reads_tot=args.nbr_reads





dico_info_expression=dict()
tmp2=open("fichier_info_expression.txt","w")






#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP -1 : liste du nbr total de reads par individus + réalisation en tête

#Dictionnaire des homéocontigs (orienté A - B)
nbr_reads_tot=list()
for line in open(reads_tot):
	nb=line.split()[1]
	nbr_reads_tot.append(nb)










#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 0 : Réalisation des dictionnaires à partir de ces fichiers.

#Dictionnaire des homéocontigs (orienté A - B)
dico_homeo=dict()
nb=0
for line in open(homeogenes):
	#Je saute la ligne d'en tete
	nb+=1
	if nb>1:
		line=line.strip()
		line=line.split()
		
		#Je ne garde que les cas ou j'ai des couples d'homéologues, pas les cas ou un contig ressemble à plusieurs autres
		if line[3] == "homoeolog_one2one" :
			contig1=line[0]
			contig2=line[2]
	
			#Dictionnaire des homéologues:
			if contig1[7]=="A" and contig1 not in dico_homeo:
				dico_homeo[contig1]=contig2	
			if contig1[7]=="B" and contig2 not in dico_homeo:
				dico_homeo[contig2]=contig1
		
			#Fichier d'information
			dico_info_expression[contig1]=contig2
			dico_info_expression[contig2]=contig1		
	
			

#Je commence a remplir mon fichier d'information
toprint="nombre de couple d'homéologue one2one total = "+str(len(dico_homeo))+"\n"
tmp2.write(toprint)




print "step0 done, bravo"




#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 1 : Ajout des infos de blast
dico_position_align=dict()


for line in open(blast):
	line=line.strip()
	line=line.split()
	
	#Je récupère toutes les infos de ma ligne de blast : nom du contig , début de l'alignement dans le contig, fin de l'alignement dans le contig.
	contig1=line[0] ; contig1=contig1.split("|")[0] 
	contig2=line[1] ; contig2=contig2.split("|")[0] 
	deb1=int(line[2]) ; end1=int(line[3]) ; seq_aligned1=line[11]
	deb2=int(line[5]) ; end2=int(line[6]) ; seq_aligned2=line[12]
	similarity=float(line[9])
	
		
	#Attention, Je vais garder uniquement les cas ou le contig a blasté sur le bon homéologue, avec le A a gauche et le B a droite
	#Est ce que le premier contig de la ligne est un contig A de la liste des homéologues? Dans ce cas la je récupére le contig B théorique
	#Plus filtre de similarité
	if contig1 in dico_homeo:
		contigB= dico_homeo[contig1]
		
		#Le contig 2 est il bien le contig théorique? Dans ce cas la uniquement je récupère les infos
		if contig2 == contigB and similarity > 90 :
			dico_position_align[contig1]=str(deb1)+"\t"+str(end1)+"\t"+seq_aligned1
			dico_position_align[contig2]=str(deb2)+"\t"+str(end2)+"\t"+seq_aligned2


#OK ! maintenant, pour chaque couple d'homéocontigs, je connais le début et la fin du chevauchement des 2.


#Je rempli mon fichier d'information
toprint="nombre de contigs pour lesquels j'ai l'info de blast = "+str(len(dico_position_align))+"\n"
tmp2.write(toprint)
nbr=float(len(dico_position_align)/2)
toprint="soit en nombre de couple = "+str(nbr)+"\n"
tmp2.write(toprint)


print "step1done mais c'est pas fini "







#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 2 : Dico des séquences fasta
dico_fasta=dict()

for line in open(fasta):
	line=line.strip()
	
	if line.startswith(">"):
		contig=line.replace(">","")
		contig=contig.split("|")[0]

	
	if contig in dico_fasta :
		dico_fasta[contig]=dico_fasta[contig]+line
	else :
		dico_fasta[contig]=""
	
	


print "step2 done, ça fait combien sur 12 ?"
	
	
	
	
	
	
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 3  ; je trouve la position des sites divergents ! et je les mets dans un dico
dico_site_divergent=dict()
nbr_total_site_divergent=0
nbr_site_analysed=0
nb_couple_contig_analysed=0


#Pour chaque couple d'homéocontigs
for contig in dico_homeo:
	nb_couple_contig_analysed+=1
	
	#Je récupère le nom des contigs A et B
	contigA=contig
	contigB=dico_homeo[contig]


	#Je récupère les séquences des contigs A et B ENTIERES
	sequenceA=dico_fasta[contigA]
	sequenceB=dico_fasta[contigB]


	#Je coupe ces 2 séquences de manière à avoir uniquement les parties chevauchantes
	#Attention il faut vérifier qu'il y a bien eu blast etre les 2 homéo contigs...
	#En fait je récupère directement la séquence alignée donnée par blastn
	if contigA in dico_position_align  :	
		debA=dico_position_align[contigA].split()[0]  ;   endA=dico_position_align[contigA].split()[1]
		sequenceA_aligned=dico_position_align[contigA].split()[2]

		if contigB in dico_position_align  :
			debB=dico_position_align[contigB].split()[0]  ;   endB=dico_position_align[contigB].split()[1]
			sequenceB_aligned=dico_position_align[contigB].split()[2]
	
	

			#Je lis alors nucléotide par nucléotide pour voir si c'est divergent ou pas. Attention, un gap = divergence
			nbr_gapA=0 ; nbr_gapB=0
			for nucl in range(1, len(sequenceA_aligned)):
				nbr_site_analysed+=1
				
				A=sequenceA_aligned[nucl-1]
				B=sequenceB_aligned[nucl-1]
					
				#Je relève les positions divergentes. Attention dans le cas d'un gap, il ne faut pas compter la position dans le nombre de position total.
				if A!=B:
					nbr_total_site_divergent+=1
				
					if A=="-" :
						nbr_gapA+=1
					if B=="-" :
						nbr_gapB+=1
						
					
					posA=int(debA)-1+nucl-nbr_gapA
					posB=int(debB)-1+nucl-nbr_gapB
					
					
					#J'enregistre chaque positions divergentes dans un dictionnaire.	
					if contigA not in dico_site_divergent :
						dico_site_divergent[contigA] = [ posA]
					else :
						dico_site_divergent[contigA].append(posA)
					
					if contigB not in dico_site_divergent :
						dico_site_divergent[contigB] = [ posB]
					else :
						dico_site_divergent[contigB].append(posB)
					
						
						
			#Je note ces infos dans mon fichier d'information.		
			if contigA in dico_site_divergent :
				dico_info_expression[contigA]=dico_info_expression[contigA]+"\t"+str(len(sequenceA))+"\t"+str(len(dico_site_divergent[contigA])) +"\t@@\t"+str(dico_site_divergent[contigA])
			else :
				dico_info_expression[contigA]=dico_info_expression[contigA]+"\t"+str(len(sequenceA))+"\t"+"0"+"\t@@\t"+"-"
			
			if contigB in dico_site_divergent :
				dico_info_expression[contigB]=dico_info_expression[contigB]+"\t"+str(len(sequenceB))+"\t"+str(len(dico_site_divergent[contigB]))+"\t@@\t"+str(dico_site_divergent[contigB])
			else :
				dico_info_expression[contigB]=dico_info_expression[contigB]+"\t"+str(len(sequenceB))+"\t"+"0"+"\t@@\t"+"-"
		
		

#Inscription dans le fichier d'informations des principales caractéristiques des calculs				
toprint1="nombre de couple de contig dans lesquels on a cherché des sites divergents =  "+str(nb_couple_contig_analysed)+"\n"
toprint1bis="nbr de contigs dans lesquels on a trouvé au moins 1 site divergent = "+str(len(dico_site_divergent))+"\n"
toprint2="nombre de sites divergents total = "+str(nbr_total_site_divergent)+"\n"
toprint3="nombre de sites analysés au total = "+str(nbr_site_analysed)+"\n"
proportion_de_site_divergent=float(nbr_total_site_divergent)/float(nbr_site_analysed)*100
toprint4="proportion moyenne de site divergents (sur cent) = "+str(proportion_de_site_divergent)+"\n"
nbr_moyen_site_divergent_par_couple=float(nbr_total_site_divergent)/float(nb_couple_contig_analysed)
toprint5="nombre moyen de sites divergent par paire de contig = "+str(nbr_moyen_site_divergent_par_couple)+"\n"
longueur_moyenne_alignment=float(nbr_site_analysed)/float(nb_couple_contig_analysed)
toprint5="longueur moyenne de l'alignement des contigs = "+str(longueur_moyenne_alignment)+"\n"
toprint=toprint1+toprint1bis+toprint2+toprint3+toprint4+toprint5+"\n\nliste des positions divergentes de chaque contig : \n"
tmp2.write(toprint)
					


for i in dico_info_expression :
	toprint=i+"\t"+str(dico_info_expression[i])+"\n"					
	tmp2.write(toprint)
	


print "step3 done, youhouhou"











#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 4  ; Je vais chercher le nombre de reads ayant mappé à chaque position !

beginning="yes"
tmp=open(output,"w")
nb_contig_avec_RPDS=0



#Je parcours mon fichier ALR
for line in open(alr) :
	line.strip()
	
	#Si je change de contig :
	if line.startswith(">") :
	

		#------1/ Je finalise le précédent, sauf dans le cas du début
		if beginning != "yes" and contig in dico_site_divergent :
				
				nb_contig_avec_RPDS += 1
		
		
				#je moyenne le nombre de reads par le nombre de site
				liste_des_comptes2=  [float(z)/len(liste) for z in liste_des_comptes] 

				#Je pondère chaque élément de la liste par le nbr de reads total et je multiplie par 1 million
				for i in range(0,177):
					liste_des_comptes2[i]=float(liste_des_comptes2[i])*1000000 / float(nbr_reads_tot[i])
				
				toprint=str(contig)
				for element in liste_des_comptes2:
					toprint=toprint+"\t"+str("{0:.2f}".format(element))
				toprint=toprint+"\n"
				tmp.write(toprint)

	
	
	
		#------2/ Je prépare le suivant je récupère la nouvelle liste des locus à récupérer si le contig est un homéo
		beginning="not_anymore"
		
		contig=line.replace(">","")
		contig=contig.split("|")[0]
		print contig
		
		#Je réinitialise
		num_ligne=-1
		liste_des_comptes=list()
		for i in range(0,177):
			liste_des_comptes.append(0)
		
		#Je récupère la liste des sites à récupérer
		if contig in dico_site_divergent :
			liste=dico_site_divergent[contig]	
			
			
	
	
	#Si je ne change pas de contig :
	else:
		line=line.split()
	
		#Si j'ai récupérer une liste, alors je vais chercher les lignes de la liste
		if contig in dico_site_divergent :
			num_ligne+=1
	
			#Si je suis dans un site divergent, j'ajoute le nombre de reads pour chaque individu.
			if num_ligne in liste :
				for col in range(2,len(line)):
					if line[1]=="P":
						nbr_reads=line[col].split("[")[0]
					if line[1]=="M":
						nbr_reads=line[col]
					
					liste_des_comptes[col-2]=liste_des_comptes[col-2]+int(nbr_reads)
				
				
		
tmp.close()


tmp2.write("nbr de contig avec un RPDS calculé = "+str(nb_contig_avec_RPDS))


print "step4 done ptit mel contente"













