#!/bin/bash
# arguments: 
#	$1 = chromosome of study
#	$2 = name of the map #1 
#	$3 = name of the map #2

#Format d'entree de la carte a respecter précisemment : 
#group   marker  position
#1A      wPt-5367        0.0
#1A      wPt-667153      0.6



#Je récupère MétaQTL
export CLASSPATH=${CLASSPATH}:/Users/holtz/Documents/PROGRAMMES/METAQTL/metaqtl.jar

#Je créé tout mon répertoir de travail
if [ -d TMP ] ; then rm -r TMP ; fi
if [ -d OUTPUT ] ; then rm -r OUTPUT ; fi
mkdir TMP
mkdir TMP/MAPS
mkdir TMP/xml
mkdir OUTPUT
cd TMP







# ----# SI J'AI UN ARGUMENT, ALORS J'AFFICHE JUSTE LA CARTE

if (($# == 2 )) ; then
	
	
	#message d'avis
	echo " Vous n'avez indiqué qu'une seule carte a afficher !!"	

	#Copie du fichier
	cp ../$2 MAPS/
	
	#Nom du fichier 
	fic1=$(echo $2 | sed 's/.*\///')

	#Fichier exp.txt
	echo -e "map""\t""mapping.function""\t""mapping.unit""\t""cross.type""\t""cross.name""\t""cross.size" > exp.txt
	echo -e $fic1"\t""kosambi""\t""cM""\t""SF2""\t""cross1""\t""150" >> exp.txt

	#Je transforme mes cartes en format xml
	java org.metaqtl.main.MetaDB -e exp.txt -m MAPS/ -o xml

	#Comparaison avant vs après insertion des baits
	java org.metaqtl.main.MMapView  -c $1  -m xml/$fic1.xml  -o ../OUTPUT/my_map  -p /Users/holtz/Documents/PROGRAMMES/METAQTL/comparison.par

	fi









# ----# SI J'AI 2 ARGUMENTS, ALORS JE FAIS UNE COMPARAISON


if (($# == 3 )) ; then


	#message d'avis
	echo "Vous avez donné 2 cartes --> réalisation d'une comparaison"

	#Copie des fichiers
	cp ../$2 MAPS
	cp $3 MAPS

	#Nom du fichier
	fic1=$(echo $2 | sed 's/.*\///')
	fic2=$(echo $3 | sed 's/.*\///')
	
	#Fichier exp.txt
	echo -e "map""\t""mapping.function""\t""mapping.unit""\t""cross.type""\t""cross.name""\t""cross.size" > exp.txt
	echo -e $fic1"\t""kosambi""\t""cM""\t""SF2""\t""cross1""\t""150" >> exp.txt
	echo -e $fic2"\t""kosambi""\t""cM""\t""SF2""\t""cross1""\t""150" >> exp.txt
	
	#Je transforme mes cartes en format xml
	java org.metaqtl.main.MetaDB -e exp.txt -m MAPS/ -o xml
	
	#Comparaison avant vs après insertion des baits
	java org.metaqtl.main.MMapView  -c $1  -r xml/$fic1.xml -m xml/$fic2.xml  -o ../OUTPUT/my_comparison -p /Users/holtz/Documents/PROGRAMMES/METAQTL/comparison.par
	
	cd ..
	#rm -r TMP
	fi
