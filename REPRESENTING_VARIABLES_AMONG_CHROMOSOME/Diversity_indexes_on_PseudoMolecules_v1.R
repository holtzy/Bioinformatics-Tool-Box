


########################################
#------ REPRESENTATION SUR LE 3B------##
########################################


diversity_calculation = function(fichier, taille_fenetre){
data=read.table(fichier , sep="\t" , header=TRUE , na.strings="-" , dec =".")
attach(data)
numeric=c(2,3,5,6,7,8,11:42)
factor=c(1,17)
for (i in numeric){data[,i]=as.numeric(as.character(data[,i]))}
for (i in factor){data[,i]=as.factor(as.character(data[,i]))}

#--3/ Représentation
#Déterminer la taille des segments de chromosome (K) d'orge à considérer.
taille_fenetre=as.numeric(taille_fenetre)

	#Récupération des données pour lequelles on a une position sur le 3B
	don=data[!is.na(data[,42]),]
	max_position=max(don[,42])
	
	#Calcul des Pi, PiN et PiS pondéré par la taille des CDS
	don$Pi_par_base=don$CDS_Pi 
	don$PN_par_base=don$CDS_PN 		
	don$PS_par_base=don$CDS_PS 

	#Création du tableau bilan de ce K. Combien de ligne faut il?
	nb_ligne=(max_position-(max_position%%taille_fenetre))/taille_fenetre+1
	AA=matrix(0,nb_ligne,15)
	
	#Remplissage du tableau du chromosome d'étude (AA)
	line=0
	for (fenetre in seq(0,(max_position-1),taille_fenetre)){										#On va faire des fenetres glissantes : on calcule la moyenne des variables pour les tronçons 1 et 2, puis les tronçcons 2 et 3, etc...
		line=line+1

		#Calcul des bornes
		min=fenetre
		max=fenetre+20*taille_fenetre

		#Taille total des CDS pour le tronçon		
		CDS_size_troncon=sum(don$CDS_size[don$position_orge< max & don$position_orge>=min] , na.rm=TRUE)

		#Je fais un sous ensemble de don concernant ce troncon
		don_troncon=don[don$position_3B< max & don$position_3B>=min , ]
		
		#Calcul des variable pour le tronçon et remplissage du tableau AA de résultat
		AA[line,1]=min/1000000
		AA[line,2]=max/1000000
		AA[line,3]=length(don_troncon$longueur_en_pb)
		AA[line,4]=sum(don_troncon$nbr_snp , na.rm=TRUE)
		AA[line,5]=sum(don_troncon$longueur_en_pb , na.rm=TRUE )
		AA[line,7]=mean(don_troncon$CDS_Pi , na.rm=TRUE) 
		AA[line,8]=mean(don_troncon$FIS_moyen , na.rm=TRUE)
		AA[line,9]=mean(don_troncon$RPKM_mean , na.rm=TRUE)
		AA[line,10]=mean(don_troncon$CDS_PS , na.rm=TRUE) 
		AA[line,11]=mean(don_troncon$CDS_PN , na.rm=TRUE) 
		AA[line,12]=mean(don_troncon$CDS_D , na.rm=TRUE)
		AA[line,13]=mean(don_troncon$he_moyen , na.rm=TRUE)
		AA[line,15]=mean(don_troncon$CDS_GC3 , na.rm=TRUE)
		}
		
	#Calcul du nbr de snp par nucleotide blasté
	AA[,6]=AA[,4]/AA[,5]*1000
	#Calcul PN/PS sur CDS 
	AA[,14]=AA[,11]/AA[,10]	

	#Nom des colonnes
	colnames(AA)=c("position sur l'orge","fin de la fenetre","nbr de contigs presents","nbr de snp detectes","nbr de nucléotides EPO correspondant","densite SNP (/1000pb)","Pi moyen","FIS moyen","RPKM moyen","PISmoy","PIN moy","Tajima","He moy","PiN / PiS","GC3")
	#Je supprime la dernière ligne qui est incomplète.
	AA=AA[ -c((nrow(AA)-20) : nrow(AA)), ]
	
return(AA)
}	
	






------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
