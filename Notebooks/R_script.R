library(SuperExactTest)
data("eqtls")
(length.gene.sets=sapply(cis.eqtls,length))
total=18196
(p=sapply(0:101,function(i) dpsets(i, length.gene.sets, n=total)))
print(p)