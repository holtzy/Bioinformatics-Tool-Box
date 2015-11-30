#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON Get_Cds_From_Pfas.py
#     				Permet d'obtenir les CDS et UTR d'un fichier .pfas (sortie du programme geno2pfas) avec l'aide d'un fichier d'annotation (.psql, sortie de prot4est)
#  					Yan Holtz, yan1166@hotmail.com
#-------------------------------------

import os
import sys
import re
try:
  import argparse
except ImportError:
  print"oops, the import /argparse/ didn't work"


parser = argparse.ArgumentParser(description= 'permet d\'utiliser les infos d\'un fichier pfas et d\'un psql pour soustirer les CDS et UTR')
parser.add_argument('-pfas', required=True, help=' fichier .pfas d\'entrée')
parser.add_argument('-psql', required=True, help=' fichier d\'annotation .psql d\'entrée')


args = parser.parse_args()

fic_pfas=args.pfas
fic_psql=args.psql







#------------------------------------------------------------------------------------------------------
print "-----\n\nLe dictionnaire des annotations est en cours de réalisation ... "


# STEP_1 : crééer mon dictionnaire du fichier .psql. Pour chaque contig j'inscrirai une liste de 3 éléments : début, fin et orientation.
info_contig=dict() ; nb_contig_annote=0 ; nb_frameshift=0 ; nb_contig_utilisable=0

for line in open(fic_psql):
	nb_contig_annote+=1
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
	Contig=line[0]
	
	#Quel est le début?
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
		info_contig[Contig]=(start,end,orientation)
		nb_contig_utilisable+=1
	else:
		nb_frameshift +=1		
		

print "-----\n\nLe dictionnaire des annotations a été réalisé avec succès\n\tNbr de contigs annotés : "+str(nb_contig_annote)+"\n\tNbr de Contig ayant un frameshift : "+str(nb_frameshift)+"\n\tNbr de Contigs utilisables = "+str(nb_contig_utilisable)+"\n\n----\n"



#------------------------------------------------------------------------------------------------------

#STEP_2 : Je crée une fonction qui lit la ligne suivante et retourne le nouveau nom de contig si il y en a un
def getContigName(fic,contig):
	line=fic.readline()
	line=line.strip()
	if line.startswith(">"):
		contig=line.replace(">","")
		contig=re.sub(r"\|original.*",r"",contig)
		contig=re.sub(r"\|like.*",r"",contig)
		contig=re.sub(r"\|compl.*",r"",contig)
		contig=contig.replace("|","_")
		return contig,line
	if line=="":										#Pour la fin quand il n'y a plus de ligne
		contig=""
		return contig,line
	else:
		return contig,line



#Je créé une fonction qui me rogne une séquence à partir des données d'annotation
def cutSequence(sequence,contig):
	start_CDS=info_contig[contig][0] ; 	end_CDS=info_contig[contig][1] ; 	sens=info_contig[contig][2]
	
	#print start_CDS,end_CDS,sens
	#print "@"+sequence
	
	#Cas ou la séquence est dans le bon sens
	#Attention en python on commence en 0. Donc début = deb-1
	#Attention, avec les crochet, donc on dit a quel nucl on s'arrete, donc fin=fin -1 +1
	if sens=="sens":
	
		#Je connais le début de mon CDS. Pour le début de l'UTR5, il faut avoir le même cadre de lecture, donc il faut que ce soit un multiple de trois. Le -1, c'est pour couper la séquence en python : premier nucléotide = position 0
		start_UTR5=(int(start_CDS)-1)%3 
	
		sequence_UTR5=sequence[start_UTR5:int(start_CDS)-1]
		sequence_CDS=sequence[int(start_CDS)-1:int(end_CDS)]
		sequence_UTR3=sequence[int(end_CDS):]
			
	#Cas ou la séquence est dans le mauvais sens : je reverse complement simplement ma séquence.
	if sens=="antisens":
	
		#Attention, comme on commence par couper la séquence, on tombe sur le UTR3 d abord. Et le UTR5 ne devra pas aller jusqu'au bout!
		nbr_nucl_a_exclure=((len(sequence)-int(end_CDS))%3)
		start_UTR5=len(sequence)-nbr_nucl_a_exclure-1
		
		sequence_UTR5=sequence[int(end_CDS):start_UTR5+1]
		sequence_CDS=sequence[int(start_CDS)-1:int(end_CDS)]
		sequence_UTR3=sequence[0:int(start_CDS)-1]
		
		#Je dois reverse complémenté mes séquences :
		sequence_UTR5=multiple_replace(sequence_UTR5, {'A':'T', 'T':'A', 'C':'G', 'G':'C'})
		sequence_CDS=multiple_replace(sequence_CDS, {'A':'T', 'T':'A', 'C':'G', 'G':'C'})
		sequence_UTR3=multiple_replace(sequence_UTR3, {'A':'T', 'T':'A', 'C':'G', 'G':'C'})
		sequence_UTR5=sequence_UTR5[::-1]
		sequence_CDS=sequence_CDS[::-1]
		sequence_UTR3=sequence_UTR3[::-1]
	
		print 0
		print start_CDS
		print end_CDS
		print start_UTR5
	
	return sequence_UTR5,sequence_CDS,sequence_UTR3
	
	

def multiple_replace(string, rep_dict):
    pattern = re.compile("|".join([re.escape(k) for k in rep_dict.keys()]), re.M)
    return pattern.sub(lambda x: rep_dict[x.group(0)], string)
    		

print "-----\n\nLes fonctions sont prêtes ... ... "







#------------------------------------------------------------------------------------------------------
print "-----\n\nFabrication des CDS et des UTR ... ... "


# STEP_3 : Je lis mon fichier et je rogne les séquences

#Je créé mon fichier de sortie
tmp_UTR5=open("EPO_106_UTR5.pfas","w")
tmp_CDS=open("EPO_106_CDS.pfas","w")
tmp_UTR3=open("EPO_106_UTR3.pfas","w")

#J'ouvre mon fichier .pfas
fic_pfas=open(fic_pfas)


#J'initialise en lisant la première ligne à la main.
deb="" ; contig="" ; save=""
info=getContigName(fic_pfas,contig)
contig=info[0] ; line=info[1]




#Ensuite je vais lire mon fichier jusqu'à la dernière ligne, c'est à dire tant que contig != à rien
while contig != "" :
	
	#Si je commence un nouveu contig
	if line.startswith(">"):
	
		#Je commence par finaliser l'ancien a l'aide de la fonction cutSequence
		#Cas particulier pour le début = on ne le fait pas au premier passage
		if deb!= "" :

			if save in info_contig:
				in_info_contig="yes"
				sequence_UTR5=cutSequence(sequence,save)[0] ; sequence_CDS=cutSequence(sequence,save)[1] ; sequence_UTR3=cutSequence(sequence,save)[2]
				tmp_UTR5.write(sequence_UTR5+"\n") ; tmp_CDS.write(sequence_CDS+"\n") ; tmp_UTR3.write(sequence_UTR3+"\n")
			else:
				in_info_contig="no"
						
		
		
		#Puis je prépare le nouveau
		deb="non"								#On a passé l'étape du premier nucléotide
		sequence=""								#Je remets la séquence à 0
		if contig in info_contig:
			tmp_UTR5.write(line+"\n") ; tmp_CDS.write(line+"\n") ; tmp_UTR3.write(line+"\n")
			save=contig
			in_info_contig="yes"
		else:
			in_info_contig="no"
			save=""

	#Si je suis toujours dans le même, et que ce contig est dans le dictionnaire
	elif in_info_contig=="yes":
		sequence=sequence+line


	#Je récupère les infos de la ligne suivante
	info=getContigName(fic_pfas,contig)
	contig=info[0] ; line=info[1]
	
	#Cas ou on est a la derniere ligne
	if contig=="":
		if save in info_contig:
			sequence_UTR5=cutSequence(sequence,save)[0] ; sequence_CDS=cutSequence(sequence,save)[1] ; sequence_UTR3=cutSequence(sequence,save)[2]
			tmp_UTR5.write(sequence_UTR5+"\n") ; tmp_CDS.write(sequence_CDS+"\n") ; tmp_UTR3.write(sequence_UTR3+"\n")
			
			
tmp_UTR5.close()
tmp_CDS.close()
tmp_UTR3.close()


print "\n\n-----\n\nC'est dans la poche !! \n\n------- "









