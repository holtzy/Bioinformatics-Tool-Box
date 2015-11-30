	
#---------------------------------------------------------------------------------------------------
	
# Chargement des Arguments Arguments
args <- commandArgs(trailingOnly = TRUE)
diversity_file=args[1]
position_file=args[2]
taille_fenetre=args[3] ; taille_fenetre=as.numeric(taille_fenetre)
glissement=args[4] ; glissement=as.numeric(glissement)
output=args[5]


# Chargement des données
data1=read.table(diversity_file , sep=" " , header=TRUE , na.strings="-" , dec =".")
data2=read.table(position_file , sep=" " , header=TRUE , na.strings="-" , dec =".")

#raccourcissment du chromosome (plus de Short ou Long)
data2[,2]=as.factor(substr(data2[,2] , 1 , 2))

#Je merge les 2 fichiers :
data=merge(data2 , data1 , by.x=1 , by.y=1 )
data[,2]=droplevels(data[,2])

# On va taffer K par K pour tous les chromosomes présents dans la colonne 2 du fichier position. 
# Pour chaque K, on créé un tableau bilan. on place ce tableau dans une liste qui s'appelle output
Bilan=""
for (chromosome in levels(data[,2])){
	
	# Message d'avis :
	print(paste(" On attaque le chromosome ",chromosome , sep=""))
	
	#Récupération des données du K d'étude:
	don=data[which(data[,2]==chromosome)  ,]													#sous ensemble des données pour un K en particulier
	max_position=max(don[,3])															#Position du dernier contig blasté = taille du K de l'orge.
	
	#Création du tableau bilan de ce K. Combien de ligne faut il?
	nb_ligne=ceiling((max_position-(taille_fenetre-1)) / glissement)
	AA=as.data.frame(matrix(0,nb_ligne,13))
	
	#Remplissage du tableau du chromosome d'étude (AA)
	line=0
	for (deb in seq(0,(max_position-taille_fenetre+1),glissement)){										
		line=line+1

		#Calcul des bornes et insertion dans le tableau bilan
		min=deb
		max=deb+taille_fenetre
		AA[line,1]=chromosome
		AA[line,2]=min
		AA[line,3]=max

		#Sous ensemble de données pour ces bornes
		fic_tmp = don[don[,3]< max & don[,3] >= min , ]
		
		#Calcul des variable pour le tronçon
		AA[line,4]=nrow(fic_tmp)
		AA[line,5]=sum(fic_tmp$S_Pop1 , na.rm=TRUE)
		AA[line,6]=sum(fic_tmp$size , na.rm=TRUE )
		AA[line,7]=mean(fic_tmp$Pi_Pop1 , na.rm=TRUE) 
		AA[line,8]=mean(fic_tmp$D_Pop1[fic_tmp$S_Pop1 != 0] , na.rm=TRUE) 
		AA[line,9]=mean(fic_tmp$PS_Pop1 , na.rm=TRUE) 
		AA[line,10]=mean(fic_tmp$PN_Pop1 , na.rm=TRUE) 
		AA[line,11]=mean(fic_tmp$gc3 , na.rm=TRUE) 
		}
		
	#Calcul du nbr de snp total par le nombre de base total (densité en SNP)
	AA[,12]=AA[,5]/AA[,6]*1000
	
	#Calcul PN/PS sur CDS 
	AA[,13]=AA[,10]/AA[,9]	

	#Ajout au bilan
	Bilan=rbind(Bilan,AA)

	}

#Nom des colonnes
colnames(Bilan)=c("chromosome","start_window","end_window","number_of_contigs","number_of_snp","number_of_nucleotide","meanPi","mean_Tajima","PISmoy","PIN moy" , "SNP_density(/1000pb)","PiN / PiS")



#Etape finale, enregistrement du fichier obtenu
write.table(Bilan, file = output , quote = FALSE , row.names = FALSE , sep="\t")


#Message de sortie
print("ok, travail terminé ! ")


