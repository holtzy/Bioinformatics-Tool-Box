# ===============
# CALCUL DES LODS SEUILS DES QTLS
# 	Ce script utilise rqtl pour calculer les LODs seuils par permutations (1000 permutations)
#	En entrée il faut lui donner une carte, un fichier de géno, et un fichier de phéno.
# ===============


# ===== RECUPERATION DE TROIS FICHIERS D'ENTREES

# --- récup fichier
args <- commandArgs(trailingOnly = TRUE)
fic_geno=args[1]
fic_pheno=args[2] 
fic_map=args[3]
col_to_study=as.numeric(args[4])

# --- Génotypage
geno <- read.table(fic_geno, sep = ";" , header = F, na.strings = "-")
geno=as.matrix(geno)
colnames(geno)=geno[1,]
geno=as.data.frame(geno[-1 , ])
print("--- Your genotyping matrix looks correct. Dimension of the matrix are :")
print(dim(geno))

# --- Carte 
# Format de la carte : 3 colonnes : LG, nom du marqueur, position dans le LG
map <- read.table(fic_map , header=T , dec = ".", na.strings = "-" , check.names=F)[ , c(1:3)]
colnames(map) <- c("LG", "marqueur", "Distance")
rownames(map) <- map$marqueur
map$LG <- as.factor(map$LG)
print("--- Your genetic map looks correct. Dimension of the map are :")
print(dim(map))

# --- Phénotypage
Y=read.table(file = fic_pheno, header = TRUE, sep = ";", dec = ".", na.strings = "NA")
my_carac=colnames(Y)[col_to_study]
Y=Y[ , c(1,col_to_study)]
colnames(Y)=c("geno","pheno")
print("--- Your Phenotyping matrix looks correct. Dimension of the matrix are :")
print(dim(Y))


# ===== FORMAT R-QTL
# fichier genfile
donnees <- merge(geno, Y , by.x=1 , by.y=1 )
indiv <- donnees[,1]
tmarq <- t(donnees[, -ncol(donnees)])[-1, ]
colnames(tmarq) <- donnees[,1]
tmarq <- data.frame(marqueur = rownames(tmarq), tmarq)
fich <- merge(map, tmarq, sort = FALSE)
tfich <- t(fich)
tfich <- data.frame(ID = c("ID", "", "", rownames(tfich)[4:length(rownames(tfich))]), tfich)
write.table(tfich, file = "tmpgen.csv", row.names = FALSE, sep = ";", col.names = FALSE,quote = FALSE)
# fichier phefile
phen <- donnees[, c(1, ncol(donnees) )]
colnames(phen)[1] <- "ID"
write.table(phen, file = "tmpphe.csv", row.names = FALSE, sep = ";", col.names = TRUE,  quote = FALSE)
	
# ===== CALCUL DU THRESHOLD VIA rQTL
library(qtl)
set.seed(123)
donnees <- read.cross("csvs", genfile = "tmpgen.csv", phefile = "tmpphe.csv" , sep=";") 
donnees <- convert2riself(donnees)
data=calc.genoprob(donnees)
res_permut <- scanone(cross = data, pheno.col = 2, method = "hk", n.perm = 1000)
a=summary( res_permut, alpha = 0.05 )
print(my_carac)
print(a)

