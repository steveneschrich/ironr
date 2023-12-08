
find_median.affybatch<-function(abatch) {

  # We do a Tukey's Biweight summarization without background, etc.
  eset<-affy::computeExprSet(abatch, pmcorrect.method="mas", summary.method="mas")
  return(find_median.matrix(exprs(eset), mean_center=FALSE, use_log=TRUE))

}

