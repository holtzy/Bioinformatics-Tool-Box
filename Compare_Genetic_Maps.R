#---------------------------------------------------------------------------------------------------------------#
#
#	 				Compare_Genetic_Maps.R
#
#	 					Script permettant de comparer n carte génétiques.
#	
#----------------------------------------------------------------------------------------------------------------#

# Format du fichier d'entrée : 
# Pour faire fonctionner ce script, taper : 


# --- Récupération des arguments : Arguments
args <- commandArgs(trailingOnly = TRUE)
nb_de_carte=length(args)
print(nb_de_carte)


# --- Je charge les n cartes dans une liste Et je fais une liste de toutes mes cartes
my_maps=list(read.table(args[1] , header=T , dec="." ))
for(i in c(2:nb_de_carte)){
	my_maps[[length(my_maps)+1]]=read.table(args[i] , header=T , dec="." )
}
# Donc si je veux des infos sur ma premiere carte je fais par example : nrow(my_maps[[1]])
# Et si je veux le nom de ma premiere carte : print(args[i])


# --- Merge the maps together
data=merge(my_maps[[1]] , my_maps[[2]], by.x=2 , by.y=2 , all=T)
colnames(data)=c("marker",paste("chromo",args[1],sep="_") , paste("pos",args[1],sep="_") , paste("chromo",args[2],sep="_") , paste("pos",args[2],sep="_"))
if(nb_de_carte>2){
	for(i in c(3:nb_de_carte)){
		data=merge(data , my_maps[[i]] , by.x=1 , by.y=2 , all=T)
		colnames(data)[c( ncol(data)-1 , ncol(data) )]= c( paste("chromo",args[i],sep="_") , paste("pos",args[i],sep="_") )
	}}
	

# --- Création d'une fonction qui a partir d'une partie de data fait le graph
my_function=function(data , nb_de_carte , chromo){
	
	#Initialiation de la carte
	par(mar=c(2,4,1,1))
	my_ylim=max(data[ , c(seq(3,ncol(data),2))] , na.rm=T) + 10
	plot(c(1:nb_de_carte) , data[1 , c(seq(3,ncol(data),2))] , type="l" , ylim=rev(c(-1,my_ylim)) , xlab="" , col=rgb(0.3,0.4,0.8,0.4) , axes=F , ylab="position en cM" , col.lab="grey" , xlim=c(0.75,nb_de_carte+0.25) )			
	axis(2 , las=2 , col="grey" , col.axis="grey"  )
	
	#J'ajoute un trait par carte
	for(i in c(1:nb_de_carte)){
		segments(i,0,i,max(data[,c((i-1)*2+3)],na.rm=T) , lwd=4)
		}
		
	#J'ajoute un trait par marqueur sur le chromosome
	for(i in c(1:nb_de_carte)){
		for(j in c(1:nrow(data))){
			segments(i-0.02 , data[j,c((i-1)*2+3)] , i+0.02 , data[j,c((i-1)*2+3)] )
		}}

	#Je relie les marqueurs communs
	for( i in c(2:nrow(data))){
		points(c(1:nb_de_carte) , data[i , c(seq(3,ncol(data),2)) ] , type="l" , col=rgb(0.3,0.4,0.8,0.4) )
		}

	#J'ajoute le nom des carte en dessous :
	for(i in c(1:nb_de_carte)){
		text(i,my_ylim+2,args[i] , col="orange" )
		}
		
	# J'ajoute le nom du chromosome + le nombre de marqueur
	nb_mark=
	text(ifelse(nb_de_carte==2,0.9,0.7) ,-4,chromo, col="orange" )

	}


# Réalisation du graph, chromosome par chromosome
pdf("my_Map_Comparison.pdf")
for(i in levels(data[ , 2] )){
	print(i)
	
	#Je dois récupérer un tableau avec tous les marqueurs ou une carte au moins a le bon chromosome ! pas facile...
	don=data[data[,2]==i & !is.na(data[,2]) , ]
	for(j in c(2:nb_de_carte)){
		temp=data[data[,c((j-1)*2+2)]==i & !is.na(data[,c((j-1)*2+2)]) , ]
		don=rbind(don,temp)
		}
	don=unique(don)
		
	
	#Je fais le graph en envoyant la fonction
	my_function(don , nb_de_carte , i)
	} 
	
	
dev.off()











