###function for correlating computed ratings and real ratings on a 75/25 split cross-validation
###Inputs:
##network: the word network. It can be 'clics3', 'omegawiki' or 'freedict'
##lexicon1: the affective ratings lexicon. It can be either 'NRCVAD' or 'WBK'
##N: number of repetitions. It must be an integer
##perc_training_set: percentage of the size of the training set in decimals 

correlations<-function(network,lexicon1,perc_training_set,N){
  ##loading of the datasets
  lexicon<-readRDS(paste0('datasets/',network,'_',lexicon1,'.Rda'))
  netw<-readRDS(paste0('datasets/',network,'.Rda'))
    
   
  
  #reconstruction of the affective lexica
  netw1<-unique(select(lexicon,from_word,Valence.x,Arousal.x,Dominance.x))
  netw1 %>% rename(Valence = Valence.x, Arousal = Arousal.x, Dominance = Dominance.x, Word = from_word)->netw1
  netw2<-unique(select(lexicon,to_word,Valence.y,Arousal.y,Dominance.y))
  netw2 %>% rename(Valence = Valence.y, Arousal = Arousal.y, Dominance = Dominance.y, Word = to_word)->netw2
  ratings<-unique(rbind(netw1,netw2)) 
  rownames(ratings)<-seq(1,nrow(ratings))
  
  cor<-data.frame(stringsAsFactors = F)
    for (k in seq(1,N))
    {
      print(k)
      
     if (perc_training_set<1) ##if the words are split into training and test set
     {
       training_set<-sample(ratings$Word,round(length(ratings$Word)*perc_training_set))
       test_set<-setdiff(ratings$Word,training_set)
     }
     else
     {
       training_set<-NA
       test_set<-c()
     }
       

       source('scripts/expansion.R')
       rat_allnetw<-expansion(network,lexicon1,training_set)
       # if (k==1)
       # {
       #  saveRDS(rat_allnetw,paste0('/home/anatale/Documents/PhD/word networks/Word ratings and colexification/paper/GitHub repository/scripts nice/results/second/second_',k,'iter_omegawiki_warriner.Rda'))
       # 
       # }
       # rat_allnetw<-readRDS(paste0('results/first_',network,'_',lexicon1,'.Rda'))
       num_words_computed<-nrow(rat_allnetw[!is.na(rat_allnetw$Valence_lang),])
       
       if(length(test_set)>=1) ##if the words are split into training and test set
       {
         computed<-rat_allnetw[rat_allnetw$Word %in% test_set,]  #I select the words that belong to the test set
         # left_out<-length(test_set)-nrow(computed)
         left_out<-(length(ratings$Word)-round(length(ratings$Word)*perc_training_set))-nrow(computed)
         nneigh<-1
         computed<-left_join(select(computed,1:10),ratings,by='Word') #add the rela ratings of the words
       }
       else
       {
         computed<-rat_allnetw[!is.na(rat_allnetw$Valence),]  ##I select the words that have groud truth ratings
         left_out<-length(which(!union(netw$from_word,netw$to_word) %in% rat_allnetw$Word))
         nneigh<-10
       }
    
       computed$num_neigh<-NA
       for (i in seq(1,length(computed$Word)))
       {
         word<-computed$Word[i]
         neigh<-union(netw$from_word[which(netw$to_word==word)],netw$to_word[which(netw$from_word==word)])
         computed$num_neigh[i]=length(neigh)
       }
       
       
       
       iter_cor<-data.frame(stringsAsFactors = F)
       df<-data.frame(nwords=NA,th_neigh=NA,val_lang=NA,c.i1_val_lang=NA,c.i2_val_lang=NA,aro_lang=NA,c.i1_aro_lang=NA,c.i2_aro_lang=NA,
                      dom_lang=NA,c.i1_dom_lang=NA,c.i2_dom_lang=NA,val_fam=NA,c.i1_val_fam=NA,c.i2_val_fam=NA,aro_fam=NA,c.i1_aro_fam=NA,
                      c.i2_aro_fam=NA,dom_fam=NA,c.i1_dom_fam=NA,c.i2_dom_fam=NA,val_mean=NA,c.i1_val_mean=NA,c.i2_val_mean=NA,
                      aro_mean=NA,c.i1_aro_mean=NA,c.i2_aro_mean=NA,dom_mean=NA,c.i1_dom_mean=NA,c.i2_dom_mean=NA,
                      left_out=left_out,stringsAsFactors = F)
       for (i in seq(1,nneigh))
       {  
       to_compute<-filter(computed,num_neigh>=i)
       
       df$nwords<-length(to_compute$Word)
       df$th_neigh<-i
       a<-cor.test(to_compute$Valence,to_compute$Valence_lang)
       df$val_lang<-a$estimate[[1]]
       df$c.i1_val_lang<-a$conf.int[1]
       df$c.i2_val_lang<-a$conf.int[2]
       
       a<-cor.test(to_compute$Arousal,to_compute$Arousal_lang)
       df$aro_lang<-a$estimate[[1]]
       df$c.i1_aro_lang<-a$conf.int[1]
       df$c.i2_aro_lang<-a$conf.int[2]
       
       a<-cor.test(to_compute$Dominance,to_compute$Dominance_lang)
       df$dom_lang<-a$estimate[[1]]
       df$c.i1_dom_lang<-a$conf.int[1]
       df$c.i2_dom_lang<-a$conf.int[2]
       
       a<-cor.test(to_compute$Valence,to_compute$Valence_fam)
       df$val_fam<-a$estimate[[1]]
       df$c.i1_val_fam<-a$conf.int[1]
       df$c.i2_val_fam<-a$conf.int[2]
       
       a<-cor.test(to_compute$Arousal,to_compute$Arousal_fam)
       df$aro_fam<-a$estimate[[1]]
       df$c.i1_aro_fam<-a$conf.int[1]
       df$c.i2_aro_fam<-a$conf.int[2]
       
       a<-cor.test(to_compute$Dominance,to_compute$Dominance_fam)
       df$dom_fam<-a$estimate[[1]]
       df$c.i1_dom_fam<-a$conf.int[1]
       df$c.i2_dom_fam<-a$conf.int[2]
       
       a<-cor.test(to_compute$Valence,to_compute$Valence_mean)
       df$val_mean<-a$estimate[[1]]
       df$c.i1_val_mean<-a$conf.int[1]
       df$c.i2_val_mean<-a$conf.int[2]
       
       a<-cor.test(to_compute$Arousal,to_compute$Arousal_mean)
       df$aro_mean<-a$estimate[[1]]
       df$c.i1_aro_mean<-a$conf.int[1]
       df$c.i2_aro_mean<-a$conf.int[2]
       
       a<-cor.test(to_compute$Dominance,to_compute$Dominance_mean)
       df$dom_mean<-a$estimate[[1]]
       df$c.i1_dom_mean<-a$conf.int[1]
       df$c.i2_dom_mean<-a$conf.int[2]
       
      iter_cor<-rbind(iter_cor,df)
       }
       cor<-rbind(cor,iter_cor)
       cor$num_words_computed<-num_words_computed
    }
  return(cor)
}