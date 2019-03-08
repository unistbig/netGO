# Build StringNetworkSum
load("StrongString.RData")
StringNetworkSum = apply(string, 1, sum)
save(StringNetworkSum, file='StrongStringNetworkSum.RData')
