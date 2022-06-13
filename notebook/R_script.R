library(SuperExactTest)
data("eqtls")
(length.gene.sets=sapply(cis.eqtls,length))
total=18196
(p=sapply(0:101,function(i) dpsets(i, length.gene.sets, n=total)))
print(p) 

common.genes=intersect(cis.eqtls[[1]], cis.eqtls[[2]], cis.eqtls[[3]], cis.eqtls[[4]])
print(common.genes)

(num.observed.overlap=length(common.genes))
cset_test=cpsets(num.observed.overlap-1,length.gene.sets,total,lower.tail = FALSE)
print(cset_test)

fit=MSET(cis.eqtls, n=total, lower.tail=FALSE)
print(fit)

res=supertest(cis.eqtls, n=total)
summary(res)
plot(res, Layout="landscape", sort.by="size",degree=c(1,2,3,4),show.overlap.size = TRUE, margin=c(0.5,5,1,2),show.set.size = FALSE,show.elements=FALSE)