awk -F '\t' '{ n11 = split($11, t11, ",")

n12 = split($12, t12, ",")

for(i = 0; ++i < n11 - 1; ){
  
  s12 = $2 + t12[i]
  
  print $1 "\t" s12 + t11[i] "\t" $2 + t12[i + 1] "\t" $4 "_intron_" i-1 "_0" "\t" $6 }

}' $1

