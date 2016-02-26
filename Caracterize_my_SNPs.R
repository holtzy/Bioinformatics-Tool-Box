#---------------------------------------------------------------------------------------------------------------#
#
#	 				CARACTERISE_MY_SNPs.R
#
#	 					Script permettant de caractériser un fichier de SNPs : graphiques, chiffres clés etc...
#	
#----------------------------------------------------------------------------------------------------------------#

#Format du fichier d'entrée : 


# Pour faire fonctionner ce script, taper : 



# Récupération des arguments : Arguments
args <- commandArgs(trailingOnly = TRUE)
fic_in=args[1]
fic_out=args[2]
fic_indiv=args[3]



#Lecture du fichier d'entrée
data=read.table(fic_in , dec="." , na.strings="-")

#Lecture du fichier donnant les individus
fic_indiv=read.table(fic_indiv )

#Ouverture du PDF de sortie
pdf(file=fic_out)

#Distribution des FIS.
hist(data[,5], xlab="FiS"  , ylab = "nbr of SNPs" , breaks=20 , col = "red" , main=paste("FIS moyen : ",round(mean(data[,5],na.rm=T),2)," | nbr of SNP : ",nrow(data),sep="")  )

#Relation entre He, FIS, et nbr d'individus génotypés :
plot(data[,1] , data[,4] , xlab="nbr d'individus homozygotes génotypés" , ylab="He")
plot(data[,4] , data[,5] , xlab="nbr d'individus homozygotes génotypés" , ylab="FIS")
plot(data[,1] , data[,5] , xlab="nbr d'individus homozygotes génotypés" , ylab="FIS")

#Nbr de snips restant selon le seuil de FIS choisis?
FIS=seq(-1,1,0.05) ; A=matrix(0,length(FIS),2)
line=0
for (i in FIS){line=line+1  ; A[line,1]=FIS[line] ; A[line,2]=nrow(data[data[,5]>i ,])}
plot(A[,2]~A[,1] , col="red" , xlab="FIS" , ylab="cumulated nbr of SNP" , main="nbr of SNP available with FIS > x" , lwd=3)

#Distribution du nombre d'individus par SNP
hist(data[,6], xlab="nombre d'individu génotypés pour le SNP"  , ylab = "nbr of SNPs" , breaks=20 , col = "grey" , main=paste("nbr individu moyen par SNP : ",round(mean(data[,6],na.rm=T),2)," | nbr of SNP : ",nrow(data),sep="")  )

#Distribution du nbr de snips par contig? --> toujours un par contig !
distr_nbr_snp_par_contig=aggregate(data[,1], by = list(contig = data[,7]), length)
hist(distr_nbr_snp_par_contig[,2] , col =22 ,   main ="distribution du nbr de snp par contig" , xlab="nbr de snp / contig" )

#Distribution du taux d'hétérozygotes
a=data[,2] / data[,6] *100
hist(a, xlab="Taux d'hétérozygotie (en %)"  , xlim=c(0,100) , ylab = "nbr of SNPs" , breaks=20 , col = "forestgreen" , main=paste("Hétérozygotie moyenne : ",round(mean(a,na.rm=T),2)," % | nbr of SNP : ",nrow(data),sep="")  )

#Distribution des He (représente l'équilibre allélique)
hist(data[,4] , xlab="He"  ,  ylab = "nbr of SNPs" , breaks=20 , col = "forestgreen" , main=paste("Hétérozygotie attendue He : ",round(mean(data[,4],na.rm=T),2)," | nbr of SNP : ",nrow(data),sep="")  )

#Nbr de SNP génotypé par individu --> Clairement certains individus qui foutent le bordel.
snp_par_indiv=data.frame(matrix(0,(ncol(data)-8)/2,2)) ; indiv=0 ; colnames(snp_par_indiv)=c("indiv" , "nbr de snp")
snp_par_indiv[,1]=fic_indiv[,1]
for (i in seq(9,ncol(data),2)){
	indiv=indiv+1
	snp_par_indiv[indiv,2]=as.numeric(length(data[!is.na(data[,i]),i]))
	}
barplot(snp_par_indiv[,2] ,names.arg =snp_par_indiv[,1] ,  col=c("Darkgreen")   , ylab="nbr de SNP génotypés pour l'individu" , xlab="individu de la pop" , main=paste("nbr de SNP genotypés par accession | moyenne = " , round(mean(as.numeric(snp_par_indiv[, 2])),2) , sep=""))
print(head(snp_par_indiv[order(snp_par_indiv[,2]) , ],10))

#Nbr de SNP disponibles avec au moins x individus
Missing=seq(0,(ncol(data)-8)/2,1) ; A=matrix(0,length(Missing),2) ; line=0
for (i in Missing){line=line+1  ; A[line,1]=Missing[line] ; A[line,2]=nrow(data[data[,6]>=i ,])}
plot(A[,2]~A[,1] , col="red" , xlab="Nbr of sample per SNP" , ylab="cumulated nbr of SNP" , main="nbr of SNP available at least x individual per SNP" , lwd=3)

#Nbr de locus couvert + Nbr de locus hétérozygote par individus
bilan=matrix(0, (ncol(data)-8)/2 , 6) ; bilan=as.data.frame(bilan) ;  colnames(bilan)=c("sample" , "tot_num_of_SNP" , "num_of_hetero" , "num_of_genot_SNP" , "%_hetero" , "%_genot")
bilan[,1]=fic_indiv[,1]
num=0
for(i in seq(10,ncol(data),2) ){
	num=num+1
	dd=data[,i]
	bilan[num,2]=length(dd)
	bilan[num,3]=length(dd[substr(dd,1,1)!=substr(dd,2,2) & !is.na(dd)])
	bilan[num,4]=length(dd[ !is.na(dd)])
	}
bilan[,5]=as.numeric(bilan[,3])/as.numeric(bilan[,4])*100
bilan[,6]=as.numeric(bilan[,4])/as.numeric(bilan[,2])*100
print(head(bilan[order(bilan[,5] , decreasing=T) , ] , 20))

#graphe nbr d'hétéro
barplot(bilan[,5] , names.arg =bilan[,1] , las=3 , main=paste("% of heterozygotie for all genotypes | mean = ",mean(bilan[,5])) , ylab="% of Heterozygotie" , col=as.numeric(bilan[,2])+1)


#Fermeture du pdf
dev.off()