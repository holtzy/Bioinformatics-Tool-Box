
#------------------------------------------------------------------------------------------------------------#
#
# 		SCRIPT : CALCULATE_LD_ADJUSTED.R
#					Permet de faire les calculs de déséquilibre de liaison corrigés par l'apparentement !
#
#-------------------------------------------------------------------------------------------------------------#


# OBJECTIF DU PROGRAMME
#Calcul la DL entre chaque paire de marqueur fournit. Calcul basique + corrigé par l'apparentement.

#FICHIER D'ENTREE
	# -- Matrice de génotype de cette forme : les allèles sont codés 0,1,2 (1pour les hétérozygotes). Les NA sont des "-"; Chaque colonne est un individu.
		#	snp_name	Tm0310_B02	Tm0406_B02	Tm0502_B02	Tm0598_B02
		#	cfn1938419	0	0	0	0
		#	cfn1958045	0	2	0	2
		#	cfn2000481	0	2	2	2
	# -- Matrice de kinship = similarité de similarté, donc la diagonal vaut 1. Demi matrice supérieure.  Matrice carré.
		#	versus	Tm0310_B02	Tm0406_B02	Tm0502_B02	Tm0598_B02
		#	Tm0310_B02	1.0	0.842880523732	0.773646801531	0.816625241646
		#	Tm0406_B02	-	1.0	0.801966675772	0.791942604857
		#	Tm0502_B02	-	-	1.0	0.795907079646
		#	Tm0598_B02	-	-	-	1.0

#Arguments
args <- commandArgs(trailingOnly = TRUE)
fic_geno=args[1]
fic_distance=args[2]
output=args[3]

#Chargement des mes données
library(LDcorSV) 

	
	#Matrice de génotype
	geno=read.table(fic_geno , sep="\t" , na.strings="-" , header=T)
	names=geno[,1] ; geno=geno[,-1]
	geno=apply(geno,2,as.numeric)
	rownames(geno)=names
	geno=t(geno)
	
	#Matrice de kinship
	kinship=read.table(fic_distance , sep="\t"  , na.strings="-" , header=T)
	names=kinship[,1]
	kinship=kinship[,-1]
	rownames(kinship)=names
	kinship=(1-kinship)
	kinship[is.na(kinship)]=0
	kinship=kinship+t(kinship)
	kinship=kinship / 2


#Calcul du DL corrigé par la matrice de Kinship ! (et observation de la matrice de kinship)
DL=LD.Measures(geno,V=kinship , S=NA )

#Impression du fichier de sortie
write.table(DL, file = output , quote = FALSE , row.names = FALSE , sep="\t")


