# Resample and return same Level Gene's Index

Resample = function(SampleGeneLV){
  res = c()
  for(i in 1:length(SampleGeneLV)){
    res = c(res, sample(which(GeneLevelIndex == names(SampleGeneLV)[i]), SampleGeneLV[i]))
  }
  return(unname(res))
}
