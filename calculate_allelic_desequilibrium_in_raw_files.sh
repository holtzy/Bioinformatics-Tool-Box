fic=$1
cat $fic | awk '{A=0 ; B=O ; for(i=1 ; i<=NF ; i++){ if($i=="A"){A=A+1} ; if($i=="B"){B=B+1} } ; if(A+B==0){He="-"}else{He=A/(A+B)} ; print $1,A,B,He }'

