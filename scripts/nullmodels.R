###function for the computation of the null models
##network can be 'clics3', 'omegawiki' or 'freedict'
##lexicon can be 'NRC VAD' or 'WBK'
##nullmodel is either 'permutation' or 'random' for random neighbors
##niter is the number of iterations of the null model
##Nneigh is the maximum number of neighbors to consider as threshold 

null_models<-function(netw,lex,ratings,corr,nullmodel,niter,Nneigh){

###########################################permutation test
if (nullmodel=='permutation')
{
  network<-readRDS(paste0('datasets/',netw,'.Rda'))
  ratings$num_neigh<-NA
  for (k in seq(1,length(ratings$Word)))
  {
    word<-ratings$Word[k]
    neigh<-union(network$from_word[which(network$to_word==word)],network$to_word[which(network$from_word==word)])
    ratings$num_neigh[k]=length(neigh)
  }
  
  test<-data.frame(val_lang_perm=NA,aro_lang_perm=NA,dom_lang_perm=NA,
                   val_fam_perm=NA,aro_fam_perm=NA,dom_fam_perm=NA,
                   val_mean_perm=NA,aro_mean_perm=NA,dom_mean_perm=NA,pvalue=NA,numneigh=seq(1,Nneigh),stringsAsFactors = F)
  val_lang<-c()
  aro_lang<-c()
  dom_lang<-c()
  val_fam<-c()
  aro_fam<-c()
  dom_fam<-c()
  val_mean<-c()
  aro_mean<-c()
  dom_mean<-c()
  for (j in seq(1,Nneigh))
  {
    print(j)
  
    to_compute<-filter(ratings,num_neigh>=j)
    
    for (i in seq(1,niter))
    {
      permuted<-sample(nrow(to_compute))
    
      a<-cor.test(to_compute$Valence[permuted],to_compute$Valence_mean)
      val_mean<-c(val_mean,a$estimate[[1]])
      
      a<-cor.test(to_compute$Arousal[permuted],to_compute$Arousal_mean)
      aro_mean<-c(aro_mean,a$estimate[[1]])
      
      a<-cor.test(to_compute$Dominance[permuted],to_compute$Dominance_mean)
      dom_mean<-c(dom_mean,a$estimate[[1]])
      
      a<-cor.test(to_compute$Valence[permuted],to_compute$Valence_lang)
      val_lang<-c(val_lang,a$estimate[[1]])
      
      a<-cor.test(to_compute$Arousal[permuted],to_compute$Arousal_lang)
      aro_lang<-c(aro_lang,a$estimate[[1]])
      
      a<-cor.test(to_compute$Dominance[permuted],to_compute$Dominance_lang)
      dom_lang<-c(dom_lang,a$estimate[[1]])
      
      a<-cor.test(to_compute$Valence[permuted],to_compute$Valence_fam)
      val_fam<-c(val_fam,a$estimate[[1]])
      
      a<-cor.test(to_compute$Arousal[permuted],to_compute$Arousal_fam)
      aro_fam<-c(aro_fam,a$estimate[[1]])
      
      a<-cor.test(to_compute$Dominance[permuted],to_compute$Dominance_fam)
      dom_fam<-c(dom_fam,a$estimate[[1]])
    }
    test$val_lang_perm[test$numneigh==j]<-list(as.vector(quantile(val_lang[(length(val_lang)-niter):length(val_lang)],probs=c(0.25,0.5,0.75))))
    test$aro_lang_perm[test$numneigh==j]<-list(as.vector(quantile(aro_lang[(length(val_lang)-niter):length(val_lang)],probs=c(0.25,0.5,0.75))))
    test$dom_lang_perm[test$numneigh==j]<-list(as.vector(quantile(dom_lang[(length(val_lang)-niter):length(val_lang)],probs=c(0.25,0.5,0.75))))
    test$val_fam_perm[test$numneigh==j]<-list(as.vector(quantile(val_fam[(length(val_lang)-niter):length(val_lang)],probs=c(0.25,0.5,0.75))))
    test$aro_fam_perm[test$numneigh==j]<-list(as.vector(quantile(aro_fam[(length(val_lang)-niter):length(val_lang)],probs=c(0.25,0.5,0.75))))
    test$dom_fam_perm[test$numneigh==j]<-list(as.vector(quantile(dom_fam[(length(val_lang)-niter):length(val_lang)],probs=c(0.25,0.5,0.75))))
    test$val_mean_perm[test$numneigh==j]<-list(as.vector(quantile(val_mean[(length(val_lang)-niter):length(val_lang)],probs=c(0.25,0.5,0.75))))
    test$aro_mean_perm[test$numneigh==j]<-list(as.vector(quantile(aro_mean[(length(val_lang)-niter):length(val_lang)],probs=c(0.25,0.5,0.75))))
    test$dom_mean_perm[test$numneigh==j]<-list(as.vector(quantile(dom_mean[(length(val_lang)-niter):length(val_lang)],probs=c(0.25,0.5,0.75))))
    
    test$pvalue[test$numneigh==j]<-list(c(sum(val_fam>=corr$val_fam[corr$th_neigh==j])/iter_perm,
                                          sum(aro_fam>=corr$aro_fam[corr$th_neigh==j])/iter_perm,
                                          sum(dom_fam>=corr$dom_fam[corr$th_neigh==j])/iter_perm,
                                          sum(val_lang>=corr$val_lang[corr$th_neigh==j])/iter_perm,
                                          sum(aro_lang>=corr$aro_lang[corr$th_neigh==j])/iter_perm,
                                          sum(dom_lang>=corr$dom_lang[corr$th_neigh==j])/iter_perm,
                                          sum(val_mean>=corr$val_mean[corr$th_neigh==j])/iter_perm,
                                          sum(aro_mean>=corr$aro_mean[corr$th_neigh==j])/iter_perm,
                                          sum(dom_mean>=corr$dom_mean[corr$th_neigh==j])/iter_perm))
  }
}


#######################random neighbors
if (nullmodel=='random')
{
  source("scripts/extract.R")
  
  network<-readRDS(paste0('datasets/',netw,'.Rda'))
  ratings$num_neigh<-NA
  for (k in seq(1,length(ratings$Word)))
  {
    word<-ratings$Word[k]
    neigh<-union(network$from_word[which(network$to_word==word)],network$to_word[which(network$from_word==word)])
    ratings$num_neigh[k]=length(neigh)
  }
  
  test<-data.frame(val_lang_rnd=NA,aro_lang_rnd=NA,dom_lang_rnd=NA,val_fam_rnd=NA,aro_fam_rnd=NA,dom_fam_rnd=NA,val_mean_rnd=NA,
                   aro_mean_rnd=NA,dom_mean_rnd=NA,pvalue=NA,numneigh=seq(1,Nneigh),stringsAsFactors = F)
  
    mat_val_mean<-matrix(nrow=niter,ncol=length(ratings$Word))
    mat_aro_mean<-matrix(nrow=niter,ncol=length(ratings$Word))
    mat_dom_mean<-matrix(nrow=niter,ncol=length(ratings$Word))
    mat_val_lang<-matrix(nrow=niter,ncol=length(ratings$Word))
    mat_aro_lang<-matrix(nrow=niter,ncol=length(ratings$Word))
    mat_dom_lang<-matrix(nrow=niter,ncol=length(ratings$Word))
    mat_val_fam<-matrix(nrow=niter,ncol=length(ratings$Word))
    mat_aro_fam<-matrix(nrow=niter,ncol=length(ratings$Word))
    mat_dom_fam<-matrix(nrow=niter,ncol=length(ratings$Word))
    
    for (i in seq(1,length(ratings$Word))) 
    {
      print(i)
      word<-ratings$Word[i]
      neigh<-union(network$from_word[network$to_word==word],network$to_word[network$from_word==word])
      neigh.num<-which(network$from_word==word | network$to_word==word)
  
      lang_weight<-network$LanguageWeight[neigh.num]
      fam_weight<-network$FamilyWeight[neigh.num]
  
  
      N_neigh <- length(neigh)
  
      neighborComp <- function(df, indices)
      { 
        return(c(mean(df$Valence[indices][1:N_neigh]), mean(df$Arousal[indices][1:N_neigh]), mean(df$Dominance[indices][1:N_neigh]),
             weighted.mean(df$Valence[indices][1:N_neigh],lang_weight), weighted.mean(df$Arousal[indices][1:N_neigh],lang_weight),
             weighted.mean(df$Dominance[indices][1:N_neigh],lang_weight),weighted.mean(df$Valence[indices][1:N_neigh],fam_weight), 
             weighted.mean(df$Arousal[indices][1:N_neigh],fam_weight),weighted.mean(df$Dominance[indices][1:N_neigh],fam_weight)))
      }
      networkrest <- ratings[ratings$Word != word,]
      bootresults <- boot(data=networkrest, R=niter, sim ="ordinary", statistic = neighborComp, weights = networkrest$num_neigh)
  
      
      mat_val_mean[,i]<-bootresults$t[,1]
      mat_aro_mean[,i]<-bootresults$t[,2]
      mat_dom_mean[,i]<-bootresults$t[,3]
      mat_val_lang[,i]<-bootresults$t[,4]
      mat_aro_lang[,i]<-bootresults$t[,5]
      mat_dom_lang[,i]<-bootresults$t[,6]
      mat_val_fam[,i]<-bootresults$t[,7]
      mat_aro_fam[,i]<-bootresults$t[,8]
      mat_dom_fam[,i]<-bootresults$t[,9]
    }
    
    
    all_words<-ratings$Word
    for(j in seq(1,Nneigh))
    {
      print(j)
      
      to_compute<-filter(ratings,num_neigh>=j)
      words<-to_compute$Word
      num.words<-which(all_words %in% words)
      
      a<-apply(mat_val_mean[,num.words],1,cor.test,to_compute$Valence[num.words])
      val_mean<-c(extract(a,niter))
      a<-apply(mat_aro_mean[,num.words],1,cor.test,to_compute$Arousal[num.words])
      aro_mean<-c(extract(a,niter))
      a<-apply(mat_dom_mean[,num.words],1,cor.test,to_compute$Dominance[num.words])
      dom_mean<-c(extract(a,niter))
      a<-apply(mat_val_fam[,num.words],1,cor.test,to_compute$Valence[num.words])
      val_fam<-c(extract(a,niter))
      a<-apply(mat_aro_fam[,num.words],1,cor.test,to_compute$Arousal[num.words])
      aro_fam<-c(extract(a,niter))
      a<-apply(mat_dom_fam[,num.words],1,cor.test,to_compute$Dominance[num.words])
      dom_fam<-c(extract(a,niter))
      a<-apply(mat_val_lang[,num.words],1,cor.test,to_compute$Valence[num.words])
      val_lang<-c(extract(a,niter))
      a<-apply(mat_aro_lang[,num.words],1,cor.test,to_compute$Arousal[num.words])
      aro_lang<-c(extract(a,niter))
      a<-apply(mat_dom_lang[,num.words],1,cor.test,to_compute$Dominance[num.words])
      dom_lang<-c(extract(a,niter))
      test$val_lang_rnd[test$numneigh==j]<-list(as.vector(quantile(val_lang,probs=c(0.25,0.5,0.75))))
      test$aro_lang_rnd[test$numneigh==j]<-list(as.vector(quantile(aro_lang,probs=c(0.25,0.5,0.75))))
      test$dom_lang_rnd[test$numneigh==j]<-list(as.vector(quantile(dom_lang,probs=c(0.25,0.5,0.75))))
      test$val_fam_rnd[test$numneigh==j]<-list(as.vector(quantile(val_fam,probs=c(0.25,0.5,0.75))))
      test$aro_fam_rnd[test$numneigh==j]<-list(as.vector(quantile(aro_fam,probs=c(0.25,0.5,0.75))))
      test$dom_fam_rnd[test$numneigh==j]<-list(as.vector(quantile(dom_fam,probs=c(0.25,0.5,0.75))))
      test$val_mean_rnd[test$numneigh==j]<-list(as.vector(quantile(val_mean,probs=c(0.25,0.5,0.75))))
      test$aro_mean_rnd[test$numneigh==j]<-list(as.vector(quantile(aro_mean,probs=c(0.25,0.5,0.75))))
      test$dom_mean_rnd[test$numneigh==j]<-list(as.vector(quantile(dom_mean,probs=c(0.25,0.5,0.75))))
      test$pvalue[test$numneigh==j]<-list(c(sum(val_fam>=corr$val_fam[corr$th_neigh==j])/iter_random,
                                            sum(aro_fam>=corr$aro_fam[corr$th_neigh==j])/iter_random,
                                            sum(dom_fam>=corr$dom_fam[corr$th_neigh==j])/iter_random,
                                            sum(val_lang>=corr$val_lang[corr$th_neigh==j])/iter_random,
                                            sum(aro_lang>=corr$aro_lang[corr$th_neigh==j])/iter_random,
                                            sum(dom_lang>=corr$dom_lang[corr$th_neigh==j])/iter_random,
                                            sum(val_mean>=corr$val_mean[corr$th_neigh==j])/iter_random,
                                            sum(aro_mean>=corr$aro_mean[corr$th_neigh==j])/iter_random,
                                            sum(dom_mean>=corr$dom_mean[corr$th_neigh==j])/iter_random))
    }
}
  return(test)
}