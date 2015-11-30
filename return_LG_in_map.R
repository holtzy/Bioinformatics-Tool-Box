


#-------------------------------------------------------------------------------------------
#	 Ce script R permet de retourner un groupe de liaison dans une carte génétique. 
#	 (Le dernier marqueur se retrouve en première position (à 0 cM). 
#	 Le nom du LG est en colonne 1, marqueur en 2, et position en 3
#--------------------------------------------------------------------------------------------

#A mettre dans un script R.


#Arguments
args <- commandArgs(trailingOnly = TRUE)
fic_map=args[1]
LG=args[2]

#Lecture du fichier d'entrée
data=read.table(fic_map , header=T , dec=".")

#Récupération longueur max du groupe de liaison a retourner
longueur=max(data[data[,1]==LG , 3])

#changement ordre marqueur
data[data[,1]==LG , 3] = longueur - data[data[,1]==LG , 3]

#je remets dans l'ordre
data=data[order(data[,1] , data[,3]) , ]

#J'enregistre la nouvelle carte obtenue
write.table(data , fic_map ,  quote = FALSE , sep="\t" , row.names=F)