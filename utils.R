library('caret')

select.by.correlation <- function(input.expre,gene.list,cut.off=0.5){
  meanLargerThanAverage <-apply(input.expre,1,sum) > 
                            median(apply(input.expre,1,sum))
  varLargerThanAverage <- apply(input.expre,1,var) > 
                            median(apply(input.expre,1,var))
  input.expre <- input.expre[meanLargerThanAverage & varLargerThanAverage,]
  input.expre.meet.both.geneID <- 
                      gene.list[meanLargerThanAverage & varLargerThanAverage]
  expre.cor <- cor(t(input.expre))
  filter.index <- findCorrelation(expre.cor,cutoff = cut.off)
  retained.gene.ID <- input.expre.meet.both.geneID[-filter.index]
  retained.gene.expre <- input.expre[-filter.index,]
  return.list <- list()
  return.list[['expression']] <- retained.gene.expre
  return.list[['gene.ID']] <- retained.gene.ID
  return(return.list)
}


para_to_String <- function(fit.result,metrics){
  #fit.result: result of a caret fitted object
  #metrics: selection criteria
  #return: a string containing parameters
  selected.model.index <- which(fit.result[metrics] == max(fit.result[metrics]))
  num.para <- which(colnames(fit.result) == metrics) - 1
  return(paste0(colnames(fit.result[1:num.para]),
                fit.result[selected.model.index,1:num.para],
                collapse = "_"))
}