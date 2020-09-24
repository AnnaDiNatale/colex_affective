expansion<-function(network,lexicon,training_set){
  ##Function that computes the ratings of the words in the networks from the original affective ratings of the NRC VAD and WBK lexica.
  ##It computes the ratings of the words as unweighted and weighted mean of the ratings of the neighbors. The weights are family and language
  ##weights.
  ##network can be "clics3", "omegawiki" or "freedict"
  ##lexicon can be "NRC VAD" or "WBK"
  ##training set is the set of words with real affective ratings. 
  
  library(dplyr)
  library(igraph)
  
  ##Loading the needed datasets
  lexicon<-readRDS(paste0('datasets/',network,'_',lexicon,'.Rda'))
  netw<-readRDS(paste0('datasets/',network,'.Rda'))
  
  
 
  
  
  #reconstructing the original affective ratings
  netw1<-unique(select(lexicon,from_word,Valence.x,Arousal.x,Dominance.x))
  netw1 %>% rename(Valence = Valence.x, Arousal = Arousal.x, Dominance = Dominance.x, Word = from_word)->netw1
  netw2<-unique(select(lexicon,to_word,Valence.y,Arousal.y,Dominance.y))
  netw2 %>% rename(Valence = Valence.y, Arousal = Arousal.y, Dominance = Dominance.y, Word = to_word)->netw2
  lexicon<-unique(rbind(netw1,netw2)) #reconstruction of the affective lexica
  rownames(lexicon)<-seq(1,nrow(lexicon))
  
  ##I reduce the groud truth lexicon keeping only the words in the training set
  if (!any(is.na(training_set)))
  {
    lexicon<-lexicon[which(lexicon$Word %in% training_set),]
  }
  
  ##Construction of the graph
  graph<-graph.data.frame(select(netw,to_word,from_word),directed=F)
  comp<-components(graph) ##connected components of the graph
  
  ##Checking that all the components have at least one word from the affective lexicon, otherwise delete the component
  to_delete<-c()
  for (i in seq(1,comp$no))
  {
    words_incomp<-groups(comp)[[i]]
    if (!any(lexicon$Word %in% words_incomp))
    {
      to_delete<-c(i,to_delete)
    }
  }

  if(length(to_delete>0))
  {
    for (i in to_delete)
    {
      words_incomp<-groups(comp)[[i]]
      netw<-netw[-which(netw$from_word %in% words_incomp | netw$to_word %in% words_incomp),]
    }
  }
  rownames(netw)<-seq(1,nrow(netw))
  
  words_network<-union(netw$from_word,netw$to_word) ##words in the filtered network
  lexicon_innetw<-intersect(lexicon$Word,words_network) ##words in the network that belong to the affective lexicon
  if (!any(is.na(training_set)))
  {
    words_tocompute<-setdiff(words_network,lexicon$Word)
  }
  else
  {
    words_tocompute<-words_network
  }

  
 ##Computing the ratings of all the words in the network 
  ratings<-data.frame(Word=words_network,
                      Valence_lang=0,Arousal_lang=0,Dominance_lang=0,
                      Valence_fam=0,Arousal_fam=0,Dominance_fam=0,
                      Valence_mean=0,Arousal_mean=0,Dominance_mean=0,stringsAsFactors = F)
  ratings<-left_join(ratings,lexicon,by='Word')
  
  
    heat<-1
    while(heat>0.001)
    {
      heat<-0
      for (i in seq(1,length(words_tocompute)))  
      {
        word<-words_tocompute[i]

        neigh<-union(netw$from_word[netw$to_word==word],netw$to_word[netw$from_word==word])
        neigh.num<-which(netw$from_word==word | netw$to_word==word)
        
        lang_weight<-netw$LanguageWeight[neigh.num]
        fam_weight<-netw$FamilyWeight[neigh.num]
        rat_val_lang<-c()
        rat_aro_lang<-c()
        rat_dom_lang<-c()

        rat_val_fam<-c()
        rat_aro_fam<-c()
        rat_dom_fam<-c()

        rat_val_mean<-c()
        rat_aro_mean<-c()
        rat_dom_mean<-c()

        for (j in neigh)
        {
          val<-unique(ratings$Valence[ratings$Word==j])
          if(!is.na(val)) ##if the word belongs to the original lexicon
          {
            aro<-unique(ratings$Arousal[ratings$Word==j])
            dom<-unique(ratings$Dominance[ratings$Word==j])
            rat_val_lang<-c(rat_val_lang,val)
            rat_aro_lang<-c(rat_aro_lang,aro)
            rat_dom_lang<-c(rat_dom_lang,dom)
            rat_val_fam<-c(rat_val_fam,val)
            rat_aro_fam<-c(rat_aro_fam,aro)
            rat_dom_fam<-c(rat_dom_fam,dom)
            rat_val_mean<-c(rat_val_mean,val)
            rat_aro_mean<-c(rat_aro_mean,aro)
            rat_dom_mean<-c(rat_dom_mean,dom)
          }
          else ##if the word does not belong to the original lexicon
          {
            val<-unique(ratings$Valence_lang[ratings$Word==j])
            aro<-unique(ratings$Arousal_lang[ratings$Word==j])
            dom<-unique(ratings$Dominance_lang[ratings$Word==j])
            rat_val_lang<-c(rat_val_lang,val)
            rat_aro_lang<-c(rat_aro_lang,aro)
            rat_dom_lang<-c(rat_dom_lang,dom)
            val<-unique(ratings$Valence_fam[ratings$Word==j])
            aro<-unique(ratings$Arousal_fam[ratings$Word==j])
            dom<-unique(ratings$Dominance_fam[ratings$Word==j])
            rat_val_fam<-c(rat_val_fam,val)
            rat_aro_fam<-c(rat_aro_fam,aro)
            rat_dom_fam<-c(rat_dom_fam,dom)
            val<-unique(ratings$Valence_mean[ratings$Word==j])
            aro<-unique(ratings$Arousal_mean[ratings$Word==j])
            dom<-unique(ratings$Dominance_mean[ratings$Word==j])
            rat_val_mean<-c(rat_val_mean,val)
            rat_aro_mean<-c(rat_aro_mean,aro)
            rat_dom_mean<-c(rat_dom_mean,dom)
            
          }
        }
        now_rat_val_lang<-weighted.mean(rat_val_lang,lang_weight)
        now_rat_aro_lang<-weighted.mean(rat_aro_lang,lang_weight)
        now_rat_dom_lang<-weighted.mean(rat_dom_lang,lang_weight)
        now_rat_val_fam<-weighted.mean(rat_val_fam,fam_weight)
        now_rat_aro_fam<-weighted.mean(rat_aro_fam,fam_weight)
        now_rat_dom_fam<-weighted.mean(rat_dom_fam,fam_weight)
        now_rat_val_mean<-mean(rat_val_mean)
        now_rat_aro_mean<-mean(rat_aro_mean)
        now_rat_dom_mean<-mean(rat_dom_mean)
        
        now_heat<-max(abs(ratings$Valence_lang[ratings$Word==word]-now_rat_val_lang),
                      abs(ratings$Arousal_lang[ratings$Word==word]-now_rat_aro_lang),
                      abs(ratings$Dominance_lang[ratings$Word==word]-now_rat_dom_lang))
        heat<-heat+now_heat
        if (now_heat>0)
        {
          ratings$Valence_lang[ratings$Word==word]<-now_rat_val_lang
          ratings$Arousal_lang[ratings$Word==word]<-now_rat_aro_lang
          ratings$Dominance_lang[ratings$Word==word]<-now_rat_dom_lang
          ratings$Valence_fam[ratings$Word==word]<-now_rat_val_fam
          ratings$Arousal_fam[ratings$Word==word]<-now_rat_aro_fam
          ratings$Dominance_fam[ratings$Word==word]<-now_rat_dom_fam
          ratings$Valence_mean[ratings$Word==word]<-now_rat_val_mean
          ratings$Arousal_mean[ratings$Word==word]<-now_rat_aro_mean
          ratings$Dominance_mean[ratings$Word==word]<-now_rat_dom_mean
        }
      }  
      print(heat)
    } 
  return(ratings)
  
}