#!/bin/bash 

#OBJECTIF DU SCRIPT
	# --> Préparation du fichier "adapter.fasta" en vue d'utiliser le programme Trimomatic
	# --> Prend en entrée la séquence du Tag de l'individu (argument 1) + La séquence de l'index (argument 2).
	# --> Attention, bien vérifier que l'adaptateur universel Ilumina utilisé est bien celui inscrit ci dessous.


#FONCTION PRELIMINAIRE : permet de faire le reverse complement d'une sequence
function rev_comp { echo $1 | sed 's/A/mm/g' | sed 's/C/ww/g' | sed 's/G/C/g' | sed 's/T/A/g' | sed 's/mm/T/g' | sed 's/ww/G/g' | awk 'BEGIN{FS=""}{ ans="" ; for(i=NF+1 ; i=i-1 ; i>0){ans=ans$i}}END{print ans}' ; }


#RECUPERATION DES TAGS ET INDEXS
seqtag=$1
seqtagrev=$(rev_comp $seqtag)
seqindex=$2
seqindexrev=$(rev_comp $2 )

#CREATION DU FICHIER

#--Les lignes Prefix_PE1 et Prefix_PE2 sont les 2 lignes qui vont servir a faire l'alignement. Il va aligner ces 2 séquences ensembles. Ensuite, si il y a Read-trough, il y aura chevauchement,et donc coupure --> on aura supprimer les adaptateurs !
echo -e ">Prefix_PE/1\nAATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"${seqtag} > adapter.fasta
echo -e ">Prefix_PE/2\nCAAGCAGAAGACGGCATACGAGAT"${seqindexrev}"GTGACTGGAGTTCAGACGTGTCGTCTTCCGATCT" >> adapter.fasta

#--Les lignes suivantes servent vont juste être "blastées"
echo -e ">PE1\nAATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"${seqtag} >> adapter.fasta
echo -e ">PE2\nCAAGCAGAAGACGGCATACGAGAT"${seqindexrev}"GTGACTGGAGTTCAGACGTGTCGTCTTCCGATCT" >> adapter.fasta
echo -e ">PE1rc\n"${seqtagrev}"AGATCGGAAGAGCACACGTCTGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTA" >> adapter.fasta
echo -e ">PE2rc\nAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"${seqindex}"ATCTCGTATGCCGTCTTCTCGTTGA" >> adapter.fasta


