
 
########################################
#------ REPRESENTATION SUR L'ORGE ------##
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

#On va taffer K par K pour les K de 1 à 7. Pour chaque K, on créé un tableau bilan. on place ce tableau dans une liste qui s'appelle "Bilan_K_par_K"
Bilan_K_par_K=list()
for (chromosome in c("1","2","3","4","5","6","7")){
	
	#Récupération des données du K d'étude:
	don=data[which(chromosome_orge==chromosome)  ,]													#sous ensemble des données pour un K en particulier
	max_position=max(don$position_orge)															#Position du dernier contig blasté = taille du K de l'orge.

	#Je peux séparer chromosome A et B ?
	don$K_ble=substr(don[,1],8,8)
	don=don[don$K_ble == "B" , ]
	print(nrow(don))
	
	#Calcul des Pi, PiN et PiS pondéré par la taille des CDS
	don$Pi_par_base=don$CDS_Pi #* don$CDS_size
	don$PN_par_base=don$CDS_PN #* don$CDS_size			
	don$PS_par_base=don$CDS_PS #* don$CDS_size

	
	#Création du tableau bilan de ce K. Combien de ligne faut il?
	nb_ligne=(max_position-(max_position%%taille_fenetre))/taille_fenetre+1
	AA=matrix(0,nb_ligne,16)
	
	#Remplissage du tableau du chromosome d'étude (AA)
	line=0
	for (fenetre in seq(0,(max_position-1),taille_fenetre)){										#On va faire des fenetres glissantes : on calcule la moyenne des variables pour les tronçons 1 et 2, puis les tronçcons 2 et 3, etc...
		line=line+1

		#Calcul des bornes
		min=fenetre
		max=fenetre+2*taille_fenetre

		#Taille total des CDS pour le tronçon		
		CDS_size_troncon=sum(don$CDS_size[don$position_orge< max & don$position_orge>=min] , na.rm=TRUE)

		#Calcul des variable pour le tronçon
		AA[line,1]=min/1000000
		AA[line,2]=max/1000000
		AA[line,3]=length(don$longueur_en_pb[don$position_orge< max & don$position_orge>=min])
		AA[line,4]=sum(don$nbr_snp[don$position_orge< max & don$position_orge>=min] , na.rm=TRUE)
		AA[line,5]=sum(don$longueur_en_pb[don$position_orge< max & don$position_orge>=min] , na.rm=TRUE )
		AA[line,7]=mean(don$CDS_Pi[don$position_orge< max & don$position_orge>=min] , na.rm=TRUE) 
		AA[line,8]=mean(don$FIS_moyen[don$position_orge< max & don$position_orge>=min] , na.rm=TRUE)
		AA[line,9]=mean(don$RPKM_mean[don$position_orge< max & don$position_orge>=min] , na.rm=TRUE)
		AA[line,10]=mean(don$CDS_PS[don$position_orge< max & don$position_orge>=min] , na.rm=TRUE) 
		AA[line,11]=mean(don$CDS_PN[don$position_orge< max & don$position_orge>=min] , na.rm=TRUE) 
		AA[line,12]=mean(don$CDS_D[don$position_orge< max & don$position_orge>=min] , na.rm=TRUE)
		AA[line,13]=mean(don$he_moyen[don$position_orge< max & don$position_orge>=min] , na.rm=TRUE)
		AA[line,15]=mean(don$CDS_GC3[don$position_orge< max & don$position_orge>=min] , na.rm=TRUE)
		AA[line,16]=mean(don$divergence[don$position_orge< max & don$position_orge>=min] , na.rm=TRUE)
		}
		
	#Calcul du nbr de snp par nucleotide blasté
	AA[,6]=AA[,4]/AA[,5]*1000
	#Calcul PN/PS sur CDS 
	AA[,14]=AA[,11]/AA[,10]	

	#Nom des colonnes
	colnames(AA)=c("position sur l'orge","fin de la fenetre","nbr de contigs presents","nbr de snp detectes","nbr de nucléotides EPO correspondant","densite SNP (/1000pb)","Pi moyen","FIS moyen","RPKM moyen","PISmoy","PIN moy","Tajima","He moy","PiN / PiS","GC3","divergence")
	#Je supprime la dernière ligne qui est incomplète.
	AA=AA[-nrow(AA),]
	
	#Stockage des résultats
	Bilan_K_par_K=c(Bilan_K_par_K,list(AA))
	}
return(Bilan_K_par_K)
}	
	










