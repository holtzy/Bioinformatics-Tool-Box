#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON Get_SNP_Infos.py
#     				Va résumer toutes les infos de divers fichiers pour chaque SNIPs
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
parser.add_argument('-SNP', required=True, help=' fichier contenant tous les SNPs | obligatoire')
parser.add_argument('-psql', required=False, help=' fichier .psql (annotation) d\'entrée')
parser.add_argument('-fasta1', required=True, help=' fichier fasta des contigs | obligatoire')
parser.add_argument('-fasta2', required=False, help=' fichier fasta des contigs (si tous les contigs ne sont pas dans le fichier 1')
parser.add_argument('-bayte', required=False, help=' fichier des baytes')
parser.add_argument('-scaffold', required=False, help=' dans quel scaffold est situé mon contig?')
parser.add_argument('-wang', required=False, help=' données de SNP de wang, elle donne une position des SNP en cM dans les scaffoldds')
parser.add_argument('-SNP_EPO', required=False, help=' fichier des SNPs EPO ')
parser.add_argument('-rep_element', required=False, help=' résultat de blast sur la base MIPS des éléments répétés ')
parser.add_argument('-out', required=True, help=' fichier de sortie | obligatoire')


args = parser.parse_args()
fic_snp=args.SNP
fic_psql=args.psql
fasta1=args.fasta1
fasta2=args.fasta2
bayte=args.bayte
scaffold=args.scaffold
wang=args.wang
SNP_EPO=args.SNP_EPO
rep_element=args.rep_element
out=args.out








#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
#LE CODE GENETIQUE
map = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
       "UCU":"S", "UCC":"s", "UCA":"S", "UCG":"S",
       "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
       "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
       "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
       "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
       "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
       "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
       "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
       "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
       "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
       "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
       "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
       "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
       "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
       "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------






#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
#FONCTIONS PRELIMINAIRES


#Fonction qui permet de faire plusieurs remplacement de nucléotide.
def multiple_replace(string, rep_dict):
    pattern = re.compile("|".join([re.escape(k) for k in rep_dict.keys()]), re.M)
    return pattern.sub(lambda x: rep_dict[x.group(0)], string)



#Petite fonction qui permet de tester si des caracteres sont des chiffres ou pas
def is_number(s):
    try:
        float(s) # for int, long and float
    except ValueError:
        try:
            complex(s) # for complex
        except ValueError:
            return False
    return True
  
  
  
#Fonction qui dira synonime ou non pour un SNP donné.
def is_synonimus(sequence,debut_CDS,fin_CDS,position,allele1,allele2,sens) : 
	sortie="-"
	longueur=len(sequence)
	position=int(position)
	debut_CDS=int(debut_CDS)
	fin_CDS=int(fin_CDS)
	sens=str(sens)
	
	#Je cree 2 sequences avec 2 alleles
	sequence1=list(sequence) ; sequence1[position-1]=allele1 ; sequence1="".join(sequence1)
	sequence2=list(sequence) ; sequence2[position-1]=allele2 ; sequence2="".join(sequence2)
	
	#SENS NORMAL
	if sens=="sens":
	
		#Je chope le bon codon
		for i in range(debut_CDS,longueur,3):
			if position>=i and position <i+3:

				codon1=sequence1[i-1:i-1+3] ; codon1=codon1.replace("T","U")
				codon2=sequence2[i-1:i-1+3] ; codon2=codon2.replace("T","U")
				
				if str(codon1) in map and str(codon2) in map:
					if map[str(codon1)] == map[str(codon2)]:
						sortie="S"
					else:
						sortie="NS"
	


	
	#ANTISENS --> besoin de faire le reverse complement, on met tout a l'envers.
	if sens=="antisens":
		position=longueur-position+1
		sequence1=multiple_replace(sequence1,{'A':'T', 'T':'A', 'C':'G', 'G':'C'}) ; sequence1=sequence1[::-1]
		sequence2=multiple_replace(sequence2,{'A':'T', 'T':'A', 'C':'G', 'G':'C'}) ; sequence2=sequence2[::-1]
		nouveau_debut=longueur-fin_CDS+1

		#Puis la meme chose
		for i in range(nouveau_debut,longueur,3):
			if position>=i and position <i+3:
				codon1=sequence1[i-1:i-1+3] ; codon1=codon1.replace("T","U")
				codon2=sequence2[i-1:i-1+3] ; codon2=codon2.replace("T","U")
			
				if str(codon1) in map and str(codon2) in map:
					if map[str(codon1)] == map[str(codon2)]:
						sortie="S"
					else:
						sortie="NS"

	return sortie
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------















#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 0 : Collecte de l'ensemble des SNP dans un dictionnaire, j'en profite pour récupérer quelques les variables les concernants
dico_des_snp=dict() ; tot=0

for line in open(fic_snp):
	tot+=1
	line=line.strip()
	line=line.split()
	contig=line[6].replace(">","")
	
	#Je raccourci le nom si c'est un blé tendre, car sinon c'est inbuvable
	if contig.startswith("Traes"):
		contig=contig.split("|")[0]
		
	snp_name=contig+"@"+line[7]
	
	#Récupération des valeurs des infos
	nb_homo=line[0] ; nb_hetero=line[1] ; ho=round(float(line[2]),2) ; 
	if line[3] != "-" :
		he=round(float(line[3]),2)
	else :
		he="-"
	if line[4] != "-" :
		FIS=round(float(line[4]),2)
	else:
		FIS="-"	
	nb_indiv=line[5] ; position=line[7]

	
	
	#Je récupère aussi les 2 alleles du SNP
	allele1="" ; allele2=""
	for i in range(9,len(line),2):
		nucl=line[i][0]
		if nucl in ("A","C","G","T"):
			if allele1 == "":
				allele1=nucl
			if nucl != allele1 :
				allele2=nucl
	
	dico_des_snp[snp_name]=contig+"\t"+position+"\t"+nb_homo+"\t"+nb_hetero+"\t"+str(ho)+"\t"+str(he)+"\t"+str(FIS)+"\t"+nb_indiv+"\t["+allele1+"/"+allele2+"]"+"\t"+"-"

entete="contig_name"+"\t"+"SNP_position"+"\t"+"nb_homo"+"\t"+"nb_hetero"+"\t"+"Ho"+"\t"+"He"+"\t"+"FIS"+"\t"+"#_ind"+"\t"+"Alleles"+"\t"+"Alleles"	
		
print "\n-------\n\nTous les SNPs ont été répertorié, il y en a "+str(tot)+" \n\n------\n"
print "\n-------\n\nleurs infos basiques ont été récupérées \n\n------\n"
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------










#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
# LE SNP est il présent chez EPO (optionnel). Si oui quel est son FIS EPO ?

if SNP_EPO is not None:
	dico_EPO=dict()
	for line in open(SNP_EPO):
		line=line.strip()
		line=line.split()
		contig=line[6].replace(">","")
		FIS=line[4]
	
		#Je raccourci le nom si c'est un blé tendre, car sinon c'est inbuvable
		if contig.startswith("Traes"):
			contig=contig.split("|")[0]
		
		snp_name=contig+"@"+line[7]
		dico_EPO[snp_name]=FIS


	#Ce SNP est il présent chez EPO?
	for snp_name in dico_des_snp:
		if snp_name in dico_EPO:
			to_print="yes"+"\t"+str(dico_EPO[snp_name])
		else:
			to_print="no"+"\t"+"-"
		dico_des_snp[snp_name]=dico_des_snp[snp_name]+"\t"+to_print
		
	entete=entete+"\t"+"presence_EPO"+"\t"+"FIS_EPO"

#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------






#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
dico_sequence=dict()
#STEP XX : je récupère le fasta de tous les contigs dans un dico

def read_fasta(fasta) :
	for line in open(fasta):
		line=line.strip()
		if line.startswith(">"):
			contig=line.replace(">","")
			#Je raccourci le nom si c'est un blé tendre, car sinon c'est inbuvable
			if contig.startswith("Traes"):
				contig=contig.split("|")[0]
			#et pour le blé dur je remplace les pipe par des _ a cause du fichier prot_main.psql.
			else :
				contig=contig.split("|")[0]+"_"+contig.split("|")[1]

			#J'initialise le dico
			dico_sequence[contig]=""
		else:
			dico_sequence[contig]=dico_sequence[contig]+line
	return dico_sequence
	
read_fasta(fasta1)


if fasta2 is not None:
	read_fasta(fasta2)

#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------





#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 1 : A l'aide de l'annotation, on va déterminer pour chaque contig, les bornes de ses CDS.
#Attention cette étape va etre faite en 2 fois : en effet : les SNPS ont 2 origines possibles:
#	- blé dur --> dans ce cas on ira lire le fichier .psql
#	- blé tendre --> dans ce cas les infos sont dans l'entete des contigs.
#on va créer un grand dico qui va regrouper toutes ces infos pour les 2 types de contigs. Il contiendra : début CDS, fin CDS, orientation.
dico_cds=dict() 

if fic_psql is not None :
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
dico_intron=dict()
dico_UTR=dict()



#A chaque nouvelle ligne je récupère toutes les infos : contig, positions des exons
for line in open(fic_snp):
	line.strip()
	line=line.split()[6]
	line=line.replace(">","")


	#Suis je dans le blé tendre?
	if line.startswith("Traes"):
		line=line.split("|")		
		contig=line[0]

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
					
	
			dico_cds[contig]=str(fin_UTR5+1)+"\t"+str(fin_des_exons)+"\t"+"sens"





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
			dico_cds[contig]=str(fin_UTR3+1)+"\t"+str(fin_des_exons)+"\t"+"antisens"

print "\n-------\n\nrécupération de toutes les données d'annotation = OK\n-------\n\n"
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------










#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP 2 : Je sais ou sont les CDS --> je vais voir pour chaque SNP ou il est situé !

for snp_name in dico_des_snp:
	contig_name_tmp=snp_name.split("@")[0]
	
	#Il va falloir modifier le nom du contig pour qu'il corresponde au nom de contig du dictionnaire des CDS.
	if contig_name_tmp.startswith("Cluster"):
		contig_name=contig_name_tmp.split("|")[0]+"_"+contig_name_tmp.split("|")[1]
	else:
		contig_name=contig_name_tmp
		

	
	if contig_name in dico_cds:
		deb=dico_cds[contig_name].split()[0] ; fin=dico_cds[contig_name].split()[1] ; sens=dico_cds[contig_name].split()[2]
		position=dico_des_snp[snp_name].split()[1]

		if deb != "-" :	
			issyn="-"	
			
			#Cas du CDS
			if int(position) >= int(deb) and int(position) <= int(fin):
				type="CDS"
				
				#Dans ce cas il faut voir si le SNP est synonyme ou non synonyme !
				sequence=dico_sequence[contig_name]
				issyn=is_synonimus(sequence,deb,fin,position,dico_des_snp[snp_name].split()[8],dico_des_snp[snp_name].split()[9],sens)
		
			#Cas de l'UTR		
			elif int(position) < int(deb):
				if sens=="sens":
					type="UTR5"
				else:
					type="UTR3"
			elif int(position) > int(fin):
				if sens=="sens":
					type="UTR3"
				else:
					type="UTR5"
			
			dico_des_snp[snp_name]=dico_des_snp[snp_name]+"\t"+type+"\t"+str(issyn)
		else:
			dico_des_snp[snp_name]=dico_des_snp[snp_name]+"\t"+"-"+"\t"+"-"


	else:
		dico_des_snp[snp_name]=dico_des_snp[snp_name]+"\t"+"-"+"\t"+"-"


		
entete=entete+"\t"+"SNP_annotation"+"\t"+"SNP_Status"

print "\n-------\n\nLes infos relatives aux CDS ont été récupérées\n\n------\n"
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------








#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### RECUPERATION  Données de Wang

if wang is not None:
	
	dico_wang_tmp=dict()
	dico_wang=dict()
	
	#Je fais un premier dico temporaire qui va lister pour tous les scaffolds les positions des SNPs en cM
	for line in open(wang):
		line=line.strip()
		line=line.split()
		__scaffold=line[3]
		position=line[6]
		if __scaffold not in dico_wang_tmp:
			dico_wang_tmp[__scaffold]=""
		dico_wang_tmp[__scaffold]=dico_wang_tmp[__scaffold]+"\t"+position
		
		
	#Je vais lire le dico temporaire, et calculer la moyenne des positions des SNPs de chaque scaffold : ce sera la position de mon scaffold !
	for __scaffold in dico_wang_tmp:
		list=dico_wang_tmp[__scaffold].split()
		
		#Calcul moyenne
		num=0 ; tot=0
		for i in list:
			
			if is_number(i) == True :		
				num=num+1
				tot=tot+float(i)
		if num==0 :
			continue
			
		moyenne=float(tot)/float(num)
		moyenne= float("{0:.2f}".format(moyenne))
		dico_wang[__scaffold]=moyenne

#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------






















#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### POSITIONNEMENT : Dans quel scaffold de blé tendre est on situé ?
if scaffold is not None :
	dico_blast=dict()
	
	#Récupération des scaffolds pour EPO
	for line in open(scaffold):
		line=line.strip()
		line=line.split()
		__scaffold=line[0]
		contig=line[1]
		contig=contig.replace("_Con","|Con")
		contig=contig.replace("_like","|like")
		contig=contig.replace("_compl","|compl")
		contig=contig.replace("_orig","|orig")
		contig=contig.replace("_less","|less")
		identity=line[2]
		coverage=line[3]
		
		#Je stocke les valeurs de blast dans un dico à part
		dico_blast[contig]=identity+","+coverage+","+__scaffold
		
		#Si on a déja rencontré ce contig, il faut déterminer si ce blast est meilleur ou moins bien que son précédent.
		if contig in dico_blast:
		
			#Si le nouveau contig a un meilleur % d'identité, alors on le garde d'office, il remplace alors l'ancien dans les dicos
			if identity > dico_blast[contig].split(",")[0]:
				dico_blast[contig]=identity+","+coverage+","+__scaffold
							
			#Si c est pile égal, alors on regarde la distance de couverture.
			if identity == dico_blast[contig].split(",")[0]:
				if coverage > dico_blast[contig].split(",")[1]:
					dico_blast[contig]=identity+","+coverage+","+__scaffold
	
	
	
	
	#Récupération des scaffolds pour le blé tendre (juste regarder dans le nom de la séquence).
	for line in open(fic_snp):
		line=line.strip()
		line=line.split()
		contig=line[6].replace(">","")
		
		#lorsque je suis sur un SNP blé tendre
		if contig.startswith("Traes"):
			contig=contig.split("|")[0]
		
			#récupération du scaffold
			__scaffold=line[6].split("|")[2]
			dico_blast[contig]="-"+","+"-"+","+__scaffold
	
	
	
									
	#Et ajout au fichier bilan
	for snp_name in dico_des_snp:
		contig=snp_name.split("@")[0]
		if contig in dico_blast:
			dico_des_snp[snp_name]=dico_des_snp[snp_name]+"\t"+dico_blast[contig].split(",")[2]
		else:
			dico_des_snp[snp_name]=dico_des_snp[snp_name]+"\t"+"-"
			
	
	entete=entete+"\t"+"scaffold"
	
			
	#Si j'ai les positions des scaffolds données par Wang, il faut les ajouter !
	if wang is not None :
		for snp_name in dico_des_snp:
			contig=snp_name.split("@")[0]
			if contig in dico_blast:
				scaffold=dico_blast[contig].split(",")[2]
				chromo=scaffold.split("_")[0]
				
				if scaffold in dico_wang:
					dico_des_snp[snp_name]=dico_des_snp[snp_name]+"\t"+str(chromo)+"\t"+str(dico_wang[scaffold])
				else:
					dico_des_snp[snp_name]=dico_des_snp[snp_name]+"\t"+"-"+"\t"+"-"
			else:
				dico_des_snp[snp_name]=dico_des_snp[snp_name]+"\t"+"-"+"\t"+"-"
		
		entete=entete+"\t"+"K"+"\t"+"pos_wang_cM"
		

			












#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### OPTIONNEL : Le SNP a t'il donné naissance à un baytes?


if bayte is not None:
	dico_bayte=dict()

	for line in open(bayte):
		line=line.strip()
		if line.startswith(">"):
			snp_name=line.replace(">","")
			snp_name=snp_name.split("|All")[0]
	
		dico_bayte[snp_name]="-"

	for snp_name in dico_des_snp:
		if snp_name in dico_bayte:
			dico_des_snp[snp_name]=dico_des_snp[snp_name]+"\t"+"yes"
		else:
			dico_des_snp[snp_name]=dico_des_snp[snp_name]+"\t"+"no"
			
	
	entete=entete+"\t"+"bayte"
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
	










#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### OPTIONNEL : element repetes ?

if rep_element is not None :
	
	dico_rep=dict()
	
	for line in open(rep_element) :
		line=line.strip()
		line=line.split()
		
		if line[9] > 80 and line[8] > 100 :
			contig_EPO=line[0]
			contig_MIPS=line[1].split("|")[0]
			dico_rep[contig_EPO]=contig_MIPS
			
	for snp_name in dico_des_snp:
		contig=snp_name.split("@")[0]
		if contig in dico_rep :
			dico_des_snp[snp_name]=dico_des_snp[snp_name]+"\t"+dico_rep[contig]
		else:
			dico_des_snp[snp_name]=dico_des_snp[snp_name]+"\t"+"-"

	entete=entete+"\t"+"REB"

















#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
### STEP FINAL : Je transforme mon dico en fichier de sortie!
tmp=open(out,"w")
tmp.write(entete+"\n")

for snp_name in dico_des_snp:
	tmp.write(dico_des_snp[snp_name]+"\n")
tmp.close()

print "\n-------\n\nLe fichier de sortie est prêt!\n\n------\n"
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------






