
#---------------------------------------------------------------------#
# -- REPRESENT A CORRELATION BETWEEN DISCRETE VARIABLES
#---------------------------------------------------------------------#


#---------------------------------------------------------------------#
# OBJECTIF
	#Because it is hard to study a correlation between 2 discrete values when dots overlap.
	
# FICHIER INPUT
	# A table with 2 discrete variables
	
# UTILISATION : 
	# Charge the function in R
	# source("Represent_discrete_variable_correlation.R")
	# If consider a table named "data", with the 2 discrete variables in 2 columns named "var1" and "var2", and we want the dots to be 4 times bigger when there are several data under it, let's make the plot like that :
	# represent_discrete_variable(data$var1, data$var2 , 3)
	
#---------------------------------------------------------------------#


#---------------------------------------------------------------------#
represent_discrete_variable=function(var1, var2 , coeff_bigger){
	#The "xyTable" function is going to count the occurence of each value couple
	AA=xyTable(var1,var2)
	#Let's plot it, representing the dots as big as the couple occurs
	plot(AA$x , AA$y , cex=2+AA$number*coeff_bigger  , pch=16 , col="chocolate1" , xlab= "rep1 sickness" , ylab="rep2 sickness" )
	text (AA$x , AA$y , AA$number )
	}
#---------------------------------------------------------------------#




