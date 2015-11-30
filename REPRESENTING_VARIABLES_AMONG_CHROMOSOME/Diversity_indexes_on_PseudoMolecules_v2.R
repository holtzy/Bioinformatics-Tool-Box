


	#--------------------------------------------------------------------------------------------------
	#   
	#		RSCRIPT : CALCULER LES VARIABLES PAR FENETRE GLISSANTE SUR LES K DE L'ORGE
	#
	#					Script réalisé par Yan Holtz (yan1166@hotmail.com / holtz@supagro.inra.fr)
	#
	#---------------------------------------------------------------------------------------------------


#Attention, le fichier d'entrée doit avoir les colonnes dans un ordre bien particulier !


# Chargement des Arguments Arguments
args <- commandArgs(trailingOnly = TRUE)
fichier=args[1]
taille_fenetre=args[2]
glissement=args[3]
output=args[4]



#Chargement des données
data=read.table(fichier , sep="\t" , header=TRUE , na.strings="-" , dec =".")


#Déterminer la taille des segments de chromosome (K) d'orge à considérer , et la taille du glissement.
taille_fenetre=as.numeric(taille_fenetre)
glissement=as.numeric(glissement)

#On va taffer K par K pour les K de 1 à 7. Pour chaque K, on créé un tableau bilan. on place ce tableau dans une liste qui s'appelle output
Bilan=""
for (chromosome in c("1","2","3","4","5","6","7")){
	
	# Message d'avis :
	print(paste(" On attaque le chromosome ",chromosome , sep=""))
	
	#Récupération des données du K d'étude:
	don=data[which(data$chromo_hordeum==chromosome)  ,]													#sous ensemble des données pour un K en particulier
	max_position=max(don$position_hordeum)															#Position du dernier contig blasté = taille du K de l'orge.
	
	#Création du tableau bilan de ce K. Combien de ligne faut il?
	nb_ligne=ceiling((max_position-(taille_fenetre-1)) / glissement)
	AA=matrix(0,nb_ligne,12)
	
	#Remplissage du tableau du chromosome d'étude (AA)
	line=0
	for (deb in seq(0,(max_position-taille_fenetre+1),glissement)){										
		line=line+1

		#Calcul des bornes et insertion dans le tableau bilan
		min=deb
		max=deb+taille_fenetre
		AA[line,1]=as.numeric(chromosome)
		AA[line,2]=min/1000000
		AA[line,3]=max/1000000

		#Sous ensemble de données pour ces bornes
		fic_tmp =don[don$position_hordeum< max & don$position_hordeum >= min , ]
		
		#Calcul des variable pour le tronçon
		AA[line,4]=nrow(fic_tmp)
		AA[line,5]=sum(fic_tmp$nbr_snp , na.rm=TRUE)
		AA[line,6]=sum(fic_tmp$longueur_en_pb , na.rm=TRUE )
		AA[line,8]=mean(fic_tmp$CDS_Pi , na.rm=TRUE) 
		AA[line,9]=mean(fic_tmp$FIS_moyen , na.rm=TRUE)
		AA[line,10]=mean(fic_tmp$CDS_PS , na.rm=TRUE) 
		AA[line,11]=mean(fic_tmp$CDS_PN , na.rm=TRUE) 
		}
		
	#Calcul du nbr de snp par nucleotide blasté
	AA[,7]=AA[,5]/AA[,6]*1000
	#Calcul PN/PS sur CDS 
	AA[,12]=AA[,11]/AA[,10]	

	#Ajout au bilan
	Bilan=rbind(Bilan,AA)

	}

#Nom des colonnes
colnames(Bilan)=c("chromosome","position_on_barley","end_window","number_of_contigs","number_of_snp","number_of_nucleotide","SNP_density(/1000pb)","meanPi","mean_FIS","PISmoy","PIN moy","PiN / PiS")



#Etape finale, enregistrement d'un objet R
write.table(Bilan, file = output , quote = FALSE , row.names = FALSE , sep="\t")


#Message de sortie
print("ok, travail terminé ! ")



