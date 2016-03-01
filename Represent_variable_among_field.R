
#---------------------------------------------------------------------#
# -- REPRESENTATION D'UNE VARIABLE SUR UNE PARCELLE
#---------------------------------------------------------------------#

#---------------------------------------------------------------------#
#OBJECTIF
#Cette fonction permet de représenter une variable sur une parcelle. Chaque plante est représenté par un carré. La couleur du carré dépend de la valeur de la variable

# FICHIER INPUT
	# Le fichier d'entrée est un tableau. La première colonne donne le nom de l'échantillon. Le tableau peut contenir autant de colonne que nécessaire mais 2 colonnes sont indispensables : 
		# - colonne "cadrillage" : elle donne la position géométrique de l'échantillon. On considère un quadrillage A,B,C,D.... sur 1,2,3,4..... La coordonnée est donnée sous le format "B-6" par exemple (respectez la syntaxe)
		# - colonne "analyse" : c'est la variable que l'on étudie et pour laquelle on calculera la moyenne des voisins, exemple :  note de symptôme mosaïque.

# UTILISATION : 
    # Pour charger la fonction : source("Represent_variable_among_field.R")
	# Si on considère un tableau nommé "data", avec la colonne 3 donnant le cadrillage, la variable à étudier en colonne 5 et le nom des échantillons en colonne 1, alors on peut faire : 
	# new_data=calculate_neighbour( data , 3 , 5, 1)
	# Un plot apparait
#---------------------------------------------------------------------#


#---------------------------------------------------------------------#
# REPRESENT FIELD FUNCTION
library(lattice)
represent_field=function(input, col_cadrillage, col_analyse , col_text){

	#column col_text is optionnal, check if it's filled or not
	if(missing(col_text)==T){col_text=NA}


	#J'impose la nature de mes colonnes d'intérêts
	input[,col_cadrillage]=as.character(input[,col_cadrillage])
	input[,col_analyse]=as.numeric(as.character(input[,col_analyse]))	

	#Fabrication de 2 petites fonctions qui récupèrent les coordonées x et y de la colonne de cadrillage
	fun_x=function(x){strsplit(x,"-")[[1]][1]}
	fun_y=function(x){strsplit(x,"-")[[1]][2] }

	# Récupération des informations concernant le cadrillage. (Je transforme le format A,B,C... en 1,2,3... + je retourne le vecteur pour que le 1 soit en haut (par défault la fonction levelplot me met le A en bas...)) :
	X=as.character(lapply(input[,col_cadrillage] , fun_x))
	X=match(X, sort(unique(X)))
	X=max(X)-X+1
	Y=as.numeric(lapply(input[,col_cadrillage] , fun_y))
	
	#choose features of the plot
	myPanel <- function(x, y, z, ...) {
    	panel.fill(col = "grey")
    	panel.levelplot(x,y,z, ... )
    	if (!is.na(col_text)){panel.text(x, y, input[,col_text] , cex=0.5  ) }	
    	}
	
	#Make the plot !!!
	levelplot(input[,col_analyse] ~ Y * X , main=paste("representation of ",colnames(input)[col_analyse], sep="") , panel = myPanel , col.regions = heat.colors(100)[length(heat.colors(100)):1] , scales=list( y=list(at=seq(1:max(X,na.rm=T)), labels=letters[max(X,na.rm=T):1]) , x=list(at=seq(1:max(Y,na.rm=T)), labels=seq(1:max(Y,na.rm=T))) ) )
	}
#---------------------------------------------------------------------#



#---------------------------------------------------------------------#
#TEST
# Utilisation de la fonction sur un exemple : 
#par(mfrow=c(1,2))
#a=paste(rep("sample",100) , seq(1,100) , sep="_")
#b=paste(sort(rep(c("A","B","C","D","E","F","G","H","I","J"),10)) , rep(seq(1,10) , 10) , sep="-")
#c=sample(0:5 , 100 , replace=T)
#data=data.frame(a,b,c)
#represent_field(data , 2 , 3 , 1)

#And the same remresenting the mean of the neighbourhood:
#source("Calculate_mean_neighbour.R")
#data_new=calculate_neighbour(data , 2 , 3)
#represent_field(data_new , 2 , 4)

#---------------------------------------------------------------------#


