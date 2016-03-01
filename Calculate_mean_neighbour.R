
#---------------------------------------------------------------------#
# -- CALCUL DES VARIABLES DE VOISINAGES
#---------------------------------------------------------------------#


#---------------------------------------------------------------------#
# OBJECTIF
	# Cette fonction permet de calculer la moyenne des plus proches voisins de chaque individus d'une parcelle pour une variable donnée.
	# Exemple, si j'ai un individus en position C-8, je vais calculer la moyenne de ma variable pour les positions : B-7, B-8, B-9, C-7, C-9, D-7, D-8, D-9, 
	
# FICHIER INPUT
	# Le fichier d'entrée est un tableau. La première colonne donne le nom de l'échantillon. Le tableau peut contenir autant de colonne que nécessaire mais 2 colonnes sont indispensables : 
		# - colonne "cadrillage" : elle donne la position géométrique de l'échantillon. On considère un quadrillage A,B,C,D.... sur 1,2,3,4..... La coordonnée est donnée sous le format "B-6" par exemple (respectez la syntaxe)
		# - colonne "analyse" : c'est la variable que l'on étudie et pour laquelle on calculera la moyenne des voisins, exemple :  note de symptôme mosaïque.

# UTILISATION : 
	# Pour charger la fonction dans R : 
	# source("Calculate_mean_neighbour.R")
	# Si on considère un tableau nommé "data", avec la colonne 3 donnant le cadrillage et la variable à étudier en colonne 5 alors on peut faire : 
	# new_data=calculate_neighbour( data , 3 , 5)
	# Le fichier new_data possède maintenant une nouvelle colonne nommée "voisinage" avec la valeur recherchée.

# Calcul à des ordres de voisinage supérieurs : 
	# Il est très aisé de calculer la moyenne des voisins d'ordres 2,3..etc
	# Pour ce faite changer la ligne du script :  cherche=c(-1,0,1) --> cherche=c(_2,-1,0,1,2)
#---------------------------------------------------------------------#


#---------------------------------------------------------------------#
# Fonction
calculate_neighbour=function(input,col_cadrillage,col_analyse){

	#J'impose la nature de mes colonnes d'intérêts
	input[,col_cadrillage]=as.character(input[,col_cadrillage])
	input[,col_analyse]=as.numeric(as.character(input[,col_analyse]))	

	#Fabrication de 2 petites fonctions qui récupèrent les coordonées x et y de la colonne de cadrillage
	fun_x=function(x){strsplit(x,"-")[[1]][1]}
	fun_y=function(x){strsplit(x,"-")[[1]][2] }

	# Récupération des informations concernant le cadrillage. (Je transforme le format A,B,C... en 1,2,3... :
	X=as.character(lapply(input[,col_cadrillage] , fun_x))
	X=match(X, unique(X))
	Y=as.numeric(lapply(input[,col_cadrillage] , fun_y))
	
	#Je créé une fonction qui va récupérer la valeur des voisins et calculer leur moyenne pour une ligne de input donnée
	myfun=function(line){
		around=c()
		cherche=c(-1,0,1)
		for (i in cherche){
			for (j in cherche){
				if(i!=0 || j!=0){
					a=X[line] + i
					b=Y[line] + j
					val=input[X==a & Y==b , col_analyse]
					around=c(around,val)
			}}}
		#Et on retourne la moyenne des voisins
		return(mean(around, na.rm=T))
		}
	
	#J'applique cette fonction à toutes les lignes de mon tableau
	for (line in c(1:nrow(input))){
		input$voisinage[line]=myfun(line)
		}
	return(input)
	}
#---------------------------------------------------------------------#




#---------------------------------------------------------------------#
#TEST

# Utilisation de la fonction sur un exemple : 
#a=paste(rep("sample",100) , seq(1,100) , sep="_")
#b=paste(sort(rep(c("A","B","C","D","E","F","G","H","I","J"),10)) , rep(seq(1,10) , 10) , sep="-")
#c=sample(0:5 , 100 , replace=T)
#data=data.frame(a,b,c)
#data=calculate_neighbour(data , 2 , 3)

# Vérification que le calcul est juste avec la valeur de la case C4 calculé a la main vs avec la fonction : 
#data[data[,2]=="B-5" , ]
#mean(data[data[,2] %in% c("A-4","A-5" , "A-6" , "B-4" , "B-6" , "C-4" , "C-5" , "C-6") , 3])
#---------------------------------------------------------------------#

