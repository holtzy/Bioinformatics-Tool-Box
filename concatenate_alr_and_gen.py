#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON CONCATENATE_ALR_AND_GEN
#                               permet d'utiliser les informations de fichier alr et alr.gen pour détecter les snps finaux.
#                                       La version 2 permet de ne pas utiliser le scripts de Julien qui aligne le nombre de colonne, car avec le novueau reads2snp on a tout le temps le meme nombre de colonne.
#                                       Yan Holtz, yan1166@hotmail.com
#-------------------------------------

import os
import sys
import re
try:
  import argparse
except ImportError:
  print"oops, the import /argparse/ didn't work"




parser = argparse.ArgumentParser(description= "permet d\'utiliser les infos d\'un fichier .alr et d\'un .gen pour détecter des snps.\n\nOUTPUT (.tsv):\n- Nombre d individus homozygotes\n-Nombre d\'individus hétérozygotes\n-Ho\n-He\n-FIS\n-Nombre d\'individus génotypés pour le SNP\n-Nom du contig\n-Position du SNP dans le contig\n-Couverture ind i\n-génotype ind i et probabilité associée\n")
parser.add_argument('-alr', required=True, help=' fichier .alr d\'entrée')
parser.add_argument('-gen', required=True, help=' fichier .gen d\'entrée')
parser.add_argument('-out', required=True, help='Nom du fichier de sortie contenant les snp attendu')
parser.add_argument('-cov', required=True, help='couverture minimum pour valider le génotype d\'un individu')
parser.add_argument('-pvalue', required=True, help='pvalue minimum pour valider le génotype d\'un individu')

#parser.add_argument('-toto', choices=['yes', 'no'], default='no', help='test branch neighbor ? yes or no (no by default)')

args = parser.parse_args()

fic_alr=args.alr
fic_gen=args.gen
output=args.out
cov=int(args.cov)
pvalue=float(args.pvalue)








#------------------------------------------------------------------------------------------------------

print "\n\n--------\n\nCollecte des individus\n\n------------\n\n"              


### STEP 0 : collecte de l'ensemble des individus.
for line in open(fic_alr):
                if line.startswith('maj'):
                        info=line.split()
                        
                        nbr_indiv=len(info)-2
                        break
                
                        

print "\n\nvotre jeux de données concerne "+str(nbr_indiv)+" individus !!!"
        
#------------------------------------------------------------------------------------------------------









#------------------------------------------------------------------------------------------------------

print "\n\n--------\n\nRécupération des SNPs\n\n------------\n\n"		




tmp=open(output,"w")
nb_contig=0

fichier=open(fic_gen)
for line_alr in open(fic_alr):																	#Pour chaque ligne du fichier .alr	
	nb_indiv_genotype=0
	line_gen=fichier.readline()																#Et chaque ligne du .gen correspondante
	line_alr=line_alr.strip() ; line_gen=line_gen.strip()

	if line_alr.startswith('>'):															#Si la ligne commence par ">"
		contig_name=line_alr
		nb_contig+=1
		print str(nb_contig)+"\t--\t"+contig_name																#Je récupère le nom du contig
		locus=0	
	else:
		if not line_alr.startswith('m'):													#Pour toutes les lignes de données (non en tête)
			locus=locus+1
			if line_alr.split()[1] == "P":													#Et en particulier celles qui sont polymorphes
				info_alr=line_alr.split() ; info_gen=line_gen.split()						
				nom_locus_a_conserver=contig_name							#Nom du locus = nom du contig + position dans le contig
				nucleotide=info_alr[1]														#Récupération du nucléotide majoritaire
				ligne_final=nom_locus_a_conserver+"\t"+str(locus)							#Je crée l'objet ligne final que je fais grandir petit a petit
				
				id_indiv=0
				while id_indiv < nbr_indiv:													#Pour tous mes indivs,
					nbr_reads=info_alr[id_indiv+2].split("[")[0]							#Le nombre de reads total est
					p_value=info_gen[id_indiv+1].split("|")[1]								#La pvalue du génotype est			
					if float(nbr_reads) > cov and float(p_value) > pvalue:						#Si cela respecte les seuils, j'affiche la valeur, sinon je mets un tiret
						ligne_final=ligne_final+"\t"+info_alr[id_indiv+2]+"\t"+info_gen[id_indiv+1]
						nb_indiv_genotype=nb_indiv_genotype+1								#Le nombre d'individu génotypé dans le fichier de snip augmente de 1
					else:
						ligne_final=ligne_final+"\t-\t-"
					id_indiv=id_indiv+1
					
				#A Partir de la, j'ai obtenu une ligne_final par locus, contenant seulement les infos intéressantes (> seuil)
				#Il faut a présent sélectionner les lignes intéressantes =  au moins 2 homozygotes différents
				#Au passage, on calculera le FIS !!
				AA=0 ; CC=0 ; GG=0 ; TT=0 ; AC=0 ; AG=0 ; AT=0 ; CG=0 ; CT=0 ; GT=0 ; A=0 ; C=0 ; G=0 ; T=0
				#Travail sur chaque individu de la ligne
				for i in ligne_final.split():												#A travers cette ligne finale
					geno=i.split("|")[0]
					
					#Les homozygotes
					if geno == "AA":
						AA=AA+1 ; A=A+2
					if geno == "CC":
						CC=CC+1 ; C=C+2
					if geno == "GG":
						GG=GG+1 ; G=G+2
					if geno == "TT":
						TT=TT+1 ; T=T+2
					
					#Les hétérozygotes
					if geno == "AC":
						AC=AC+1 ; A=A+1 ; C=C+1
					if geno == "AG":
						AG=AG+1 ; A=A+1 ; G=G+1
					if geno == "AT":
						AT=AT+1 ; A=A+1 ; T=T+1
					if geno == "CG":
						CG=CG+1 ; C=C+1 ; G=G+1
					if geno == "CT":
						CT=CT+1 ; C=C+1 ; T=T+1
					if geno == "GT":
						GT=GT+1 ; G=G+1 ; T=T+1
				
				#Obtention des chiffres pour le locus entier :
				#nombre d'allèle différent présent (max=4)
				nb_allele_diff=0
				if A>0:
					nb_allele_diff=nb_allele_diff+1
				if C>0:
					nb_allele_diff=nb_allele_diff+1
				if G>0:
					nb_allele_diff=nb_allele_diff+1
				if T>0:
					nb_allele_diff=nb_allele_diff+1

				#nombre d'allèle différent présent (max=4)
				nb_ind_homo_diff=0
				if AA>0:
					nb_ind_homo_diff=nb_ind_homo_diff+1
				if CC>0:
					nb_ind_homo_diff=nb_ind_homo_diff+1
				if GG>0:
					nb_ind_homo_diff=nb_ind_homo_diff+1
				if TT>0:
					nb_ind_homo_diff=nb_ind_homo_diff+1

								
				nb_hetero=AC+AG+AT+CG+CT+GT 	; nb_hetero=int(nb_hetero)
				nb_homo=AA+CC+GG+TT 			; nb_homo=int(nb_homo)
				nb_total_allele=A+C+G+T# 		; nb_total_allele=float(nb_total_allele)
				nb_indiv=nb_hetero + nb_homo	; nb_indiv=float(nb_indiv)


				#Conditions à remplir pour être un SNP : avoir au moins 2 indiv homo différents
				if int(nb_ind_homo_diff) >= 2 :													#On s'intéresse au cas ou on a au moins 2 individus homozygotes différents différents présent !

					#Calcul de He et du FIS
					ho=(nb_hetero)/(nb_indiv)																	#On calcule FIS si et seulement si on a 2 allèles différents. Si on a 3 allèles, ont affiche la ligne, mais avec des tirets pour le FIS
					if float(nb_allele_diff) ==2:       #Que pour les BI allèliques
						he=1-pow(float(A)/float(nb_total_allele),2)-pow(float(C)/float(nb_total_allele),2)-pow(float(G)/float(nb_total_allele),2)-pow(float(T)/float(nb_total_allele),2)
						FIS=(he-ho)/he
					else:
						he="-" ; FIS="-"
					
					#Remplissage du fichier de SNP
					ligne_final=str(nb_homo)+"\t"+str(nb_hetero)+"\t"+str(ho)+"\t"+str(he)+"\t"+str(FIS)+"\t"+str(nb_indiv_genotype)+"\t"+str(ligne_final)     #J'ajoute des infos a la ligne finale : nbr homo, nbr hetero et FIS
					tmp.write(ligne_final+'\n')
					
print "\n \n ----------------- \n \n Vos snips ont été obtenus avec succès !!!! \n Bonne journée et que la force soit avec vous ! \n \n -----------------"		
tmp.close()				


					
					
					
					
					
							

