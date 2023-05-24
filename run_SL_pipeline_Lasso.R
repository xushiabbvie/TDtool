### Lasso algorithm to identify SL pairs
### Data needed:
# gene_effect_scores: CERES scores from depmap portal and formatted to matrix with rows as genes and columns as cancer cell lines, genes with less than 5 essential and non-essetnial cell lines are fitlered out
# gene_mut_mat: damaging mutations from depmap portal and formatted to matrix with rows as genes and columns as cancer cell lines, genes with mut freq < 0.02 are filtered out

library(glmnet)
library(rlist)

genes_analyze = rownames(gene_effect_scores) ## genes to run for the pipeline
Lasso_model = mclapply(genes_analyze[1:5],mc.cores=10,function(x){
  x_model = cv.glmnet(t(gene_mut_mat),gene_effect_scores[x,],nfolds=10)
  return(as.matrix(coef(x_model,s="lambda.min"))[-1,1])
})

Lasso_coef = list.rbind(Lasso_model) ## Lasso coefficients used to get SL
rownames(Lasso_coef) = genes_analyze
