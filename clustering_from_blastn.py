#!/usr/bin/python
# -*- coding: utf-8 -*-


#-------------------------------------
#
#    SCRIPT PYTHON Clustering_from_blasn
#     				Take A fasta file. Make an auto-blastn of this sequences using Blastn. This script permits to build cluster of sequences that have blasted together.
#					Example : contig1 blast on contig2. Contig2 blast on contig5. This script will build a new fasta file with sequences of contig1, 2 and 5.
#  					Yan Holtz, yan1166@hotmail.com
#-------------------------------------

import os
import sys
import re
try:
  import argparse
except ImportError:
  print"oops, the import /argparse/ didn't work"


parser = argparse.ArgumentParser(description= '\n\nTake A fasta file. Make an auto-blastn of this sequences using Blastn. This script permits to build cluster of sequences that have blasted together.        ---------           Example : contig1 blast on contig2. Contig2 blast on contig5. This script will build a new fasta file with sequences of contig1, 2 and 5.\n\n')
parser.add_argument('-blast', required=True, help=' Blastn Output')
parser.add_argument('-fasta', required=True, help=' Fasta file containing the sequences used for the Blastn')
parser.add_argument('-marge', required=True, help='nb de nucleotide max avant et ares la zone de blast')
parser.add_argument('-similarity', required=False, help='seuil de similarite de blast')
parser.add_argument('-overlap', required=True, help='seuil de couverture du blast')


args = parser.parse_args()

fic_blast=args.blast
fic_fasta=args.fasta
marge=args.marge
similarity=args.similarity
overlap=args.overlap






#----------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION 0 : BLASTN OUTPUT SELECTION
# Cette fonction permet de filtrer le blastn : on considèrera uniquement les couples qui nous intéresse (longueur overlap, % identity, longueur overhang)
def blastn_output_selection(input,marge,p,overlap):
	selectionned_contigs_couple_list=[]
	id1=0 ; id2=1 ; deb1=2 ; end1=3 ; len1=4 ; deb2=5 ; end2=6 ; len2=7 ; len_align=8 ; pc_ident=9
	for line in open(input):
		infos=line.split()
		id1= infos[0]; id2=infos[1]; pc_ident=float(infos[9]); len_align=int(infos[8]);
		deb1=int(infos[2]) ; end1=int(infos[3]) ; len1=int(infos[4]) ;
		deb2=int(infos[5]) ; end2=int(infos[6]) ; len2=int(infos[7]) ;

		if pc_ident >= float(p) and len_align > int(overlap):								# On sélectionne sur % identity, et longueur overlap

			if id1 !=id2:																# On veux pas les blast sur soi même
				if deb2 < end2:															# On va tester les 6 cas de figure possible.
					chevauche= ( (len1-end1)<marge and deb2<marge )
					est_inclu_id1= ( (len1-end1)<marge and deb1<marge )
					est_inclu_id2= ( (len2-end2)<marge and deb2<marge )
					if( (chevauche or est_inclu_id1 or est_inclu_id2) == True):
						selectionned_contigs_couple_list.append([id1,id2])
				if deb2 > end2:
					chevauche= ( (len1-end1)<marge and (len2-end2)<marge )
					est_inclu_id1= ( (len1-end1)<marge and deb1<marge )
					est_inclu_id2= ( (len2-end2)<marge and deb2<marge )
					if( (chevauche or est_inclu_id1 or est_inclu_id2) == True):			# Si un des cas est OK
						selectionned_contigs_couple_list.append([id1,id2])				# alors on garde le couple
	#print selectionned_contigs_couple_list
	return selectionned_contigs_couple_list	
#----------------------------------------------------------------------------------------------------------------------------------------------
	
	





#----------------------------------------------------------------------------------------------------------------------------------------------
	
# FUNCTION 1 : BUILDING THE HASHTABLE
def building_hash_table(table):
	"""	A partir de la sortie du blastn filtre, on fabrique un dictionnaire"""
	id=1
	Code_Contig=dict()
	for line in table:
		for contig in line:
			if contig not in Code_Contig:
			 	Code_Contig[contig]=id
			 	id=id+1
	#print Code_Contig
	return Code_Contig,id
#----------------------------------------------------------------------------------------------------------------------------------------------








	
	
#--------------------------------------------------------------------------------------------

# FUNCTION 2 : LINKAGE TABLE
def building_linkage_table(Code_Contig,id,table):
	"""Creation de la table de lien"""
	Link_Table=[]
	i=1
	while i<id:
		Link_Table.append([])
		i=i+1
	for line in table:
		C1=line[0]
		C2=line[1]
		id_C1=Code_Contig[C1]
		id_C2=Code_Contig[C2]
		Link_Table[id_C1-1].append(id_C2)
		Link_Table[id_C2-1].append(id_C1)
	#print Link_Table
	return Link_Table

#----------------------------------------------------------------------------------------------------------------------------------------------






	
#--------------------------------------------------------------------------------------------	
#STEP 3 : BUILDING CLUSTERS !

def formation_des_clusters(Link_Table,id):
	###Creation de la table d'exploration : 0 si la ligne n'a pas encore ete traite, 1 sinon
	Explore_Table=[]
	i=0
	while i<id:
		Explore_Table.append(0)
		i=i+1
	
	#Obtention des Clusters.
	Num_Cluster=0
	i=0
	Bilan_Cluster=dict()
		
	while i < len (Link_Table):
		if Explore_Table[i] == 0 :
			Num_Cluster=Num_Cluster+1
			Contigs_Du_Cluster=[]
			Contigs_Du_Cluster.append(i+1)
			Contigs_A_Traiter=[]
			Contigs_A_Traiter.append(i+1)
			while len(Contigs_A_Traiter)>0:
				Contig_Courant=Contigs_A_Traiter[0]
				Bilan_Cluster[Contig_Courant]="Cluster_"+str(Num_Cluster)
				del Contigs_A_Traiter[0]
				if Explore_Table[Contig_Courant-1]==0:
					Explore_Table[Contig_Courant-1]=1
					Contigs_Du_Cluster=Contigs_Du_Cluster+Link_Table[Contig_Courant-1]
					Contigs_A_Traiter=Contigs_A_Traiter+Link_Table[Contig_Courant-1]
			#print("Cluster_"+str(i+1)+"   =   "+str(Contigs_Du_Cluster))
		i=i+1

	#print(Bilan_Cluster)
	return(Bilan_Cluster)
#----------------------------------------------------------------------------------------------------------------------------------------------
	
	
	
	
#----------------------------------------------------------------------------------------------------------------------------------------------
	
#STEP 4 : RECUPERATION DES SEQUENCES

def recuperation_sequences(Bilan_Cluster,sequences_fasta,Code_Contig):
	#Je fabrique tous les fichiers de Cluster.
	for cluster in Bilan_Cluster.values():
		open(cluster,"w")
	open("cluster_singlet","w")

	id_cluster=""

	fichier=open(sequences_fasta,"r")
	for line in fichier:
		if line[0:1]==">":
			seq_name=line[1:]
			seq_name=seq_name.replace("\n","")
			if seq_name in Code_Contig:
				id_seq_name=Code_Contig[seq_name]
				id_cluster=Bilan_Cluster[id_seq_name]
			else:
				id_cluster="cluster_singlet"
		#print(id_cluster)
		#print(line)
		tmp=open(id_cluster,"a")
		tmp.write(line)
		tmp.close
	
#----------------------------------------------------------------------------------------------------------------------------------------------








#----------------------------------------------------------------------------------------------------------------------------------------------

#LET'S USE FUNCTIONS !

#Etape 0: selectionner les couples de contigs a regrouper
result_etape0=blastn_output_selection(fic_blast, marge, similarity, overlap)
print "OK :  contigs a regrouper ensemble sélectionnés"

#Etape 1: Faire un dictionnaire de tous les contigs selectionnes en etape 0
result_etape1=building_hash_table(result_etape0)
[a,b]=result_etape1
#Donc a=Code_Contig et b=id
print "

#Etape 2: je cree ma linkage table a partir de mon dictionnaire
result_etape2=building_linkage_table( a , b, result_etape0)

#Etape 3: formation des clusters !
resultat_etape3=formation_des_clusters(result_etape2,b)

#Etape 4: recuperation des sequences correspondantes
recuperation_sequences(resultat_etape3, fic_fasta, a)

#----------------------------------------------------------------------------------------------------------------------------------------------









