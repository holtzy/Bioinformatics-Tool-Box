#------------------------------------------------------------------------------------------------------------#
#
# 		SCRIPT : MAKE A VENN DIAGRAMM
#					A partir de liste, fait un diagramm de Venn
#
#-------------------------------------------------------------------------------------------------------------#


# OBJECTIF DU PROGRAMME
# A partir de liste, fait un diagramm de Venn

#FICHIER D'ENTREE
	#des listes, peu importe leur format tant que les noms communs sont bien les mêmes !
	#Le nom de sortie demandé doit etre au format .png
	
	
#Arguments
args <- commandArgs(trailingOnly = TRUE)
fic_1=args[1]
fic_2=args[2]
fic_3=args[3]
output=args[4]

#Library de venndiagram
library(VennDiagram)

#Chargement des données
SNP_pop_1=read.table(fic_1)[,1]
SNP_pop_2=read.table(fic_2)[,1]
SNP_pop_3=read.table(fic_3)[,1]

#Venn-plot them !
venn.plot <- venn.diagram(
x = list(SNP_pop_1 , SNP_pop_2 , SNP_pop_3),
category.names = c(fic_1 , fic_2 , fic_3) ,
filename = output,
    output = TRUE ,
    height = 3000 ,
    width = 3000 ,
    resolution = 300,
    compression = "lzw",
    lwd = 5,
    lty = 'blank',
    fill = c('yellow', 'purple', 'green'),
    cex = 3.5,
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 3,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    rotation = 1
    )



