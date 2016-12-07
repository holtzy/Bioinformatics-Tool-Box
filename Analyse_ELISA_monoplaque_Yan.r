###############################
###		SCRIPT ELISA    #######
###############################
#
#

#OBJECTIF DU SCRIPT
	#Transformer les données d une manière lisible
	#Tester si la plaque est OK.
	#Sortir un bilan de la plaque facilement lisible avec eventuellement les échantillons à refaire.
	
#DEROULEMENT DU SCRIPT
	#-STEP 0 : sur quelles donnees je veux travailler : rentrer les infos de sa plaque.
	#-STEP 1 : je cree mon tableau bilan de plaque, avec les moyennes, les moyennes pondérées etc...
	#-STEP 2 : je teste la validité de ma plaque
	#-STEP 3 : je compare mes plaques pour décider quelle est la meilleure méthode de pondération.

#FICHIERS D'ENTREES DU SCRIPT
	#-->le fichier plaque1lecture1.csv
	#-->le fichier liste qui donne le contenu de chaque puits (liste des échantillons)




####################################
#   STEP 0 : POUR QUELLE PLAQUE JE VEUX TRAVAILLER
####################################



#---------------------------------
# VERO A REMPLIR !!!!!!!!!!!!!
#---------------------------------
		#Je charge mon répertoire de travail:
		#setwd("mon répertoire a véronique")

#setwd("D:/Documents and Settings/Viader/Mes documents/AAELISA2014/dataELISA2014/data_modif")
#getwd()

#setwd("H:/travail20140822/AAELISA2014/dataELISA2014/data_modif")
#getwd()
rm(list=ls())
		#J'écris mes fichiers d'entrées et de sortie:
		#fic_liste= "listeplV18bisWSSMV.csv"
		#fic_lecture= "plV18bis.csv"

args <- commandArgs(trailingOnly = TRUE)
fic_liste=args[1]
fic_lecture=args[2]



#---------------------------------

####################################
#	STEP 1 : FAIRE MON TABLEAU BILAN
####################################

#Nom du fichier de sortie.
name=paste(gsub(".csv","",fic_lecture) , "_concat.csv", sep="" )

#Je charge les données
liste=read.table(fic_liste , sep=";" , header=TRUE)
lecture=read.table(fic_lecture , dec= "." , sep= ";" , header=TRUE , na.strings="NA" )


#Je commence par faire un vecteur avec toutes les valeurs de DO.
vecteur_DO=c()

for(col in c(4,6,8,10,12)){
for (line in c(1:8)){

vecteur_DO=c(vecteur_DO,lecture[line,col])
vecteur_DO=c(vecteur_DO,lecture[line,col+1])
}}

#DO moyenne des temois + coeff variation ?
moy_blancs = round( mean(lecture[,2] , na.rm=TRUE ) ,2) 					;  coeff.var_blancs=round( (sd(lecture[,2] , na.rm=TRUE) *100  / mean(lecture[,2] , na.rm=TRUE)) , 2)					
DO_temoins_negatifs=round(mean(lecture[c(1:4),3] , na.rm=TRUE )	,2)	;  coeff.var_negatifs=round( (sd(lecture[c(1:4),3] , na.rm=TRUE ) *100  / mean(lecture[c(1:4),3] , na.rm=TRUE)) , 2)
DO_temoins_positifs=round(mean(lecture[c(5:8),3] , na.rm=TRUE )	,2)	;  coeff.var_positifs=round( (sd(lecture[c(5:8),3] , na.rm=TRUE) *100  / mean(lecture[c(5:8),3] , na.rm=TRUE)) , 2)

#Moyenne pondérée des DO
moy_ponder_neg_blanc=round( (mean(lecture[,2], na.rm=TRUE) - moy_blancs) / (DO_temoins_negatifs - moy_blancs) , 2)	;	moy_ponder_pos_blanc=round( (mean(lecture[,2]) - moy_blancs) / (DO_temoins_positifs - moy_blancs) , 2)
moy_ponder_neg_temoins_negatifs=round( (mean(lecture[c(1:4),3], na.rm=TRUE) - moy_blancs) / (DO_temoins_negatifs - moy_blancs) , 2)	;	moy_ponder_pos_temoins_negatifs=round( (mean(lecture[c(1:4),3], na.rm=TRUE) - moy_blancs) / (DO_temoins_positifs - moy_blancs) , 2)
moy_ponder_neg_temoins_positifs=round( (mean(lecture[c(5:8),3], na.rm=TRUE) - moy_blancs) / (DO_temoins_negatifs - moy_blancs) , 2)	;	moy_ponder_pos_temoins_positifs=round( (mean(lecture[c(5:8),3], na.rm=TRUE) - moy_blancs) / (DO_temoins_positifs - moy_blancs) , 2)


#Je créé et rempli ma matrice BILAN avec les moyennes de chaque indiv
bilan=matrix(0,nrow(liste),12)
colnames(bilan)=c("plaque","lecture","virus","anticorps","blanc","negatif","positif","plante","moyenne","coef_var","moy_pondere_negatif","moy_pondere_positif")

#Colonnes 1 a 5 = généralités
bilan[,1]=as.character(lecture[1,14])
bilan[,2]=as.character(lecture[2,14])
bilan[,3]=as.character(lecture[3,14])
bilan[,4]=as.character(lecture[4,14])

bilan[,5]=moy_blancs
bilan[,6]=DO_temoins_negatifs
bilan[,7]=DO_temoins_positifs

bilan[,8]=c(as.character(liste[,2]))

#Colonnes 9 a 12 = Calcul
	#Pour les plantes
num=0
for (i in c(1:nrow(liste))){ 
	num=num+1
	rep1=vecteur_DO[num]
	num=num+1
	rep2=vecteur_DO[num]
	bilan[i,9]=mean(c(rep1,rep2))
	bilan[i,10]=round( sd(c(rep1,rep2)) *100  / mean(c(rep1,rep2)) , 2)
	bilan[i,11]=round( (mean(c(rep1,rep2)) - moy_blancs) / (DO_temoins_negatifs - moy_blancs) , 2)
	bilan[i,12]=round( (mean(c(rep1,rep2)) - moy_blancs) / (DO_temoins_positifs - moy_blancs) , 2)
	}


#Quel est l'individu avec la plus grande DO (pourra servir de temoin positif si celui n'est pas pr�sent)
tmp=bilan[as.numeric(as.character(bilan[,10]))<30 & !is.na(bilan[,10]) , ]
tmp=tmp[order(as.numeric(as.character(tmp[,9])) , decreasing=T) , ]
new_max_value=as.numeric(as.character(tmp[1,9]))
new_max_value

#Si le temoin positif n'est pa valable, alors il faut recalculer toutes les moyennes ponderees avec cette nouvelle valeur.
if (is.na(DO_temoins_positifs)){ bilan[,12] = round( (as.numeric(as.character(bilan[,9])) - moy_blancs) / (new_max_value - moy_blancs) , 2)}




write.table(bilan, file = name , quote = FALSE , row.names = FALSE , sep=";" , na = "-")







####################################
#	STEP 2 : CETTE PLAQUE EST ELLE BONNE?
####################################

#DO moyenne des témois positifs >= 2* DO moyenne témoins négatifs?
rapport=DO_temoins_positifs / DO_temoins_negatifs


#Quel est la meilleur DO max de remplacement?
bilan_OK=data.frame(bilan)
bilan_OK=bilan[bilan[,7]<5,]
bilan_OK=bilan_OK[order(bilan_OK[,6] , decreasing = TRUE),]
remplacement=bilan_OK[1,]

#Quels sont les individus ayant un coeff de variation trop fort = A REFAIRE:
nb_a_refaire=length(which(as.numeric(bilan[,7])>10))
a_refaire=bilan[which(as.numeric(bilan[,7])>10),]

#petite manip si il n y a que un seul a refaire
if(nb_a_refaire==1){a_refaire=t(data.frame(a_refaire))}

name=paste("a_refaire_",fic_lecture,sep="")
write.table(a_refaire , file = name ,sep=";"   , row.names = FALSE )



#Affichage des conclusions
cat("\n\n\n\n\n\n\n\n")
print("Je travaille sur la plaque :")
print(fic_lecture)
print("le rapport DO témoins positif / DO témoins negatifs vaut :")
print(rapport)
print("Si ce rapport n'est pas bon (<3), alors je peux utiliser la DO de remplacement suivante:")
print(remplacement)
print("les indivs ayant trop de différences entre 2 mesures sont :")
print(a_refaire)
cat("\n\n\n\n\n\n\n\n")































































###USELESS : TABLEAU VERTICALE


#Combien y'a t'il de données en tout
nb_de_donnees_total=8+8+80
#Combien de colonnes?
nb_colonnes=8
#Je créé une matrice rempli de 0 avec ces dimensions
matrice_verticale=matrix(0,nb_de_donnees_total , nb_colonnes)


#Maintenant je rempli ma matrice petit à petit.
#colonne 1
matrice_verticale[,1]=lecture[1,14]
#pour les BLANCS
for (i in c(1:8)){matrice_verticale[i,2]="blanc" ; matrice_verticale[i,3]= "blanc" ; matrice_verticale[i,4]=i ; matrice_verticale[i,5]=lecture[i,2]}
#Pour les témoins
num=0
for (i in c(9:12)){num=num+1 ; matrice_verticale[i,2]="tnegatif" ; matrice_verticale[i,3]="tnegatif" ; matrice_verticale[i,4]=num ; matrice_verticale[i,5]=lecture[num,3]}
num=0
for (i in c(13:16)){num=num+1 ; matrice_verticale[i,2]="tpositif" ; matrice_verticale[i,3]="tpositif"  ; matrice_verticale[i,4]=num ; matrice_verticale[i,5]=lecture[num+4,3]}
#Pour les échantillons
for (i in c(17:nb_de_donnees_total)){matrice_verticale[i,3]="echantillon"  }
num=0
for (i in c(1:nrow(liste))){ 
	num=num+1
	matrice_verticale[ 2*i-1+16 ,2] = as.character(liste[i,2])
	matrice_verticale[ 2*i-1+16 ,4] = "1"
	matrice_verticale[ 2*i-1+16 ,5]  = vecteur_DO[num]
	num=num+1
	matrice_verticale[ 2*i+16 ,2] = as.character(liste[i,2]) ;  
	matrice_verticale[ 2*i+16 ,4] = "2"
	matrice_verticale[ 2*i+16 ,5] = vecteur_DO[num]
	}
	
#col6
matrice_verticale[,6]=as.character(lecture[8,14])
#col7
matrice_verticale[,7]=as.character(lecture[3,14])
#col8
matrice_verticale[,8]=as.character(lecture[4,14])


#------> J'ai obtenu ma matrice verticale.

