
			#--------------------------------------------------------------------------------------------------
			#   
			#		SCRIPT R : FROM GENO PHENO MAP TO BILAN_SIMPLE_MARKER
			#
			#					Holtz Yan
			#					Besnard Alban
			#---------------------------------------------------------------------------------------------------



#-- OBJECTIF 
	# Ce script permet de calculer le LOD de chaque marqueurs pour chaque variables phéno données.
	# Il sort un fichier nommé bilan_simple_marker qui servira d'entrée à l'application shiny de visualisation des QTL.

#-- FICHIER INPUT DANS L'ORDRE : 
	# Fichier de génotypage :
	# Fichier de phénotypage :
	# Carte génétique
	



# -- Récupération des Arguments
args <- commandArgs(trailingOnly = TRUE)
fic_geno=args[1]
fic_pheno=args[2] 
fic_map=args[3]





#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

# -------------------------------------------------
# PARTIE 1 : RECUPERATION DES 3 TABLEAUX D'ENTREE
# -------------------------------------------------


# --- Génotypage
# Le format du génotpage doit respecter les règles suivantes : 1/ Nom du fichier = fic_genotypage.csv 2/ en tete des colonnes = nom des geno 3/ données manquantes = "-" 4/ Séparateur=";" 
geno <- read.table(fic_geno, sep = ";" , header = F, na.strings = "-")
geno=as.matrix(geno)
colnames(geno)=geno[1,]
geno=as.data.frame(geno[-1 , ])


# --- Carte 
# Format de la carte : 3 colonnes : LG, nom du marqueur, position dans le LG
map <- read.table(fic_map , header=T , dec = ".", na.strings = "-" , check.names=F)
colnames(map) <- c("LG", "marqueur", "Distance","group_physique","Posi_physique")
rownames(map) <- map$marqueur
map$LG <- as.factor(map$LG)


# --- Phénotypage
Y=read.table(file = fic_pheno, header = TRUE, sep = ";", dec = ".", na.strings = "NA")
colnames(Y)[1]="geno"






#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

# -------------------------------------------------
# PARTIE 2 : VERIFICATION DES COMPATIBILITES
# -------------------------------------------------


# Vérification des données pour voir si tout va bien
verif=as.data.frame(matrix(0,6,2)) ; verif[,1]=c("nbr de geno dans la map" , "nbr geno dans fichier de génot" , "nbr indiv dans fichier de génot" , "nbr d'indiv dans fic phénot" , "nbr geno communs carte / genotypage","nbr indiv communs genot / phenot")
verif[1,2]=nrow(map) ; verif[2,2]=ncol(geno)-1 ; verif[3,2]=nrow(geno) ; verif[4,2]=nrow(Y) ; verif[5,2]=length(colnames(geno)[colnames(geno)%in%map$marqueur==TRUE])  ; verif[6,2]=length(Y[,1][ Y[,1]%in%geno[,1]==TRUE])
print(verif)






#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# -------------------------------------------------
# PARTIE 3 : CREATION DE FONCTION QUI CALCULENT LES LODS
# -------------------------------------------------

#Fonction SIM : pour un marqueur et un caractere : retourne pvalue, R2 moyenne etc..
sim <- function(x, y) {
  modele <- lm(y ~ x)
  variance <- anova(modele)
  moy <- aggregate(y, by = list(marqueur = x), mean, na.rm = TRUE)
  resu <- c(variance$"Pr(>F)"[1], variance$"Sum Sq"[1]/(variance$"Sum Sq"[1] + variance$"Sum Sq"[2]), moy$x[moy$marqueur == "A"], moy$x[moy$marqueur =="B"], abs(moy[1, 2] - moy[2, 2])/2)
  names(resu) <- c("pvalue", "R2", "moy.A", "moy.B", "a")
  return(resu)
}

#Fonction simple marker : pour un caractere donné (num de colonne), va retourner 2 objet : "final" = un tableau avec tous les marqueurs significatifs et leur infos (LOD, moyenne, pos...)	; "res"=idem mais pour TOUS les marqueurs
simple_marker_analysis=function(colonne){
  pheno <- Y[, c(1,colonne)]
  pheno[,2]=as.numeric(as.character(pheno[,2]))
  pheno=na.omit(pheno)
  
  donnees <- merge(geno , pheno , by.x=1 , by.y=1)
  colnames(donnees)[ncol(donnees)] <- "pheno"
  indiv <- donnees$geno
  
  #J'enleve les locus pourris
  vec=c() ; for(i in c(2:(ncol(donnees)-1))){a=nlevels(droplevels(donnees[,i])) ;  if(a!=2){vec=c(vec,i)}} ; if(length(vec)>0){donnees=donnees[ , -vec]}
  
  #On applique la fonction a tous les marqueurs (toutes les colonnes)
  res2 <- apply(donnees[, 2:(ncol(donnees) - 1)], MARGIN = 2, FUN = function(x,y) sim(x, y), y = donnees$pheno)
  res2 <- t(res2)
  res <- as.data.frame(res2)
  res$LOD=-log10(res$pvalue)
  return(res)
}

#Fonction treeshold_results pour obtenir la liste des marqueurs supérieurs au treeshold précisé
treeshold_results=function(res,LOD_seuil=4){	
  #on récupère les marqueurs ayant un effet significatif : LOD > 4
  msignif <- res[res$LOD >= LOD_seuil, ]
  msignif[,7]=row.names(msignif)
  map_signif=map[map[,2] %in% rownames(msignif), ]
  final=merge(map_signif, msignif , by.x=2 , by.y=7 , all=T)
  return(final)
}





#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# -------------------------------------------------
# PARTIE 4 : CALCUL DU BILAN_SIMPLE_MARKER
# -------------------------------------------------
to_check=names(Y)[sapply(Y,is.numeric)==TRUE]

bilan_simple_marker=data.frame()
for (i in to_check){
  a=simple_marker_analysis(which(colnames(Y) == i))
  a$marker=rownames(a)
  a=merge(map,a , by.x=2 , by.y=c("marker") , all.x=T )
  a=a[ , c(2,1,3:ncol(a))]
  if(nrow(a)>0){
    a$variable=i}
  bilan_simple_marker=rbind(bilan_simple_marker,a)
}

# écriture du bilan
write.csv(x=bilan_simple_marker,file="bilan_simple_marker",row.names=F)











