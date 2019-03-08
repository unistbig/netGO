simplifyLi = function(li1){
  l1 = list()
  for(i in 1:length(li1)){
    thisgene = names(li1)[i]
    thisgenes = li1[[i]]
    
    l1[[i]] = thisgenes[which(thisgene<thisgenes)]
  }
  names(l1) = names(li1)
  l1 = l1[-which(sapply(l1,length)==0)]
  return(l1)
}


l1 = simplifyLi(li1)
l2 = simplifyLi(li2)
l3 = simplifyLi(li3)
l4 = simplifyLi(li4)
l5 = simplifyLi(li5)
l6 = simplifyLi(li6)

save(l1,file='l1.RData')
save(l2,file='l2.RData')
save(l3,file='l3.RData')
save(l4,file='l4.RData')
save(l5,file='l5.RData')
save(l6,file='l6.RData')