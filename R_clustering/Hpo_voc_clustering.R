# https://cran.r-project.org/web/packages/textmineR/vignettes/

# https://cran.r-project.org/web/views/Cluster.html
# https://www.datanovia.com/en/blog/types-of-clustering-methods-overview-and-quick-start-r-code/
# ontology mathching + Bootstrapping + concept-recognition
# https://www.slideshare.net/samhati27/ontology-mapping-47703379


library(tidyverse)  # data manipulation
library(cluster, clustertend)    # clustering algorithms : agnes() + diana()
library(magrittr)   # Forward-Pipe Operator 
library(textmineR)  # ! Text Mining and Topic Modeling 
                    
library(NbClust)    # check cluster number
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms
library(wordcloud)

library(bgmm)       # Belief-Based Gaussian Mixture Modeling (soft-label and belief-based modeling, semi-supervised methods and fully supervised methods)

library(ggplot2)
library(ggdendro)
library(plotly)


##################     NLP and matrix construction    ##################

###  Extract all terms of defined group

# load translated hpo terms
# load("~/OneDrive/CarrierGene/2019.07.10-HPO/Hpo_ML/Hpo_reference/Jax_hpo_translated.RData")
load("C:/Users/Jingz/OneDrive/CarrierGene/2019.07.10-HPO (done)/hpo_import_python/Hpo_reference/Jax_hpo_translated_chpo.RData")

# add our own clinical information
res_all <-{}
cli <- read.csv("C:/Users/Jingz/OneDrive/CarrierGene/2019.07.10-HPO (done)/Cluster_R/analysis_output/test.csv") #clinical_informations_token.csv

cli$Panel_Type <- as.factor(cli$Panel_Type)


for (p in levels(cli$Panel_Type)) {
#res_panel <-{}
print(paste0('will analysis panel type : ', p))  


clisub <- cli[cli$Panel_Type==p,]   # select sample via panel type
row.names(clisub) <- 1:nrow(clisub)


# choose all hpo from selected class
if (clisub$Panel_Type=='精阅' | clisub$Panel_Type=='伊阅') {a=14; b=11}        # Abnormality of the genitourinary system (jingyue ok + yiyue 50% with wrong tranlation)
if (clisub$Panel_Type=='遗传性肿瘤') {a=25; b=6}                               # Neoplasm (ok)
if (clisub$Panel_Type=='子阅'| clisub$Panel_Type=='显阅')  {a=c(7,14); b=11}   # Abnormality of prenatal development or birth + Abnormality of the genitourinary system
if (clisub$Panel_Type=='全外'| clisub$Panel_Type=='嘉阅')  {a=1:25; b=11}      # {a=1:length(list.files('C:/Users/Jingz/OneDrive/CarrierGene/2019.07.10-HPO/Hpo_ML/Hpo_analysis_R')); b=1}

            
res_clisub<- {} # Abnormality of prenatal development or birth + Abnormality of the genitourinary system


# import reference hpo 

for (l in 1:length(a)){
  
setwd("C:/Users/Jingz/OneDrive/CarrierGene/2019.07.10-HPO (done)/hpo_import_python/Hpo_group_data/")
PATH<-toString(list.files()[a[l]])  # read HPO terms related to each panel type
setwd(PATH)
set.seed(300)
descendants <- read.csv(list.files()[1], sep="")
descendants <- na.omit(descendants)

if (l==1){Desc <- descendants[1,]}
Desc <-rbind(Desc, descendants)
}



df_subclass <- merge(Desc[,],Df_en, by.x='x',by.y='id',all.x = TRUE)  # merge with hpo database
rm(Desc)
df_subclass$full_text <- paste(df_subclass$name, df_subclass$synonyms, df_subclass$description, sep='.')
df_subclass <- df_subclass[c('x','name','full_text')] 
print(paste0('will analyze panel: ',p))

df_subclass <-Df_en
df_subclass$full_text <- paste(df_subclass$name, df_subclass$synonyms, df_subclass$description, sep='.')
df_subclass <-df_subclass [,c(1,2,5)]
colnames(df_subclass)<-c('x','name','full_text') 

# calculate distance and export into csv

for (i in 1:nrow(clisub)){
 
 df<- clisub[i,c(1,4)]
 sam<-df$Sample_ID
 
 colnames(df)<-c('x','full_text')
 df$name<-df$x
 tem<-merge(df_subclass,df, by=c('x','full_text','name'), all=TRUE)
 tem<-unique(tem)
 rownames(tem)<-tem$name 

 
# create a document term matrix (tokenizing)
 stopw<- c(stopwords::stopwords("en"),       # stopwords from tm
 stopwords::stopwords(source = "smart"))     # this is the default value

 stopw<- stopw[!grepl('^no$|^non$|^not$|^none$|^noone$',stopw)]

 dtm <- CreateDtm(doc_vec = tem$full_text,   # character vector of documents
                  doc_names = rownames(tem), # document names
                  stem_lemma_function = function(x) SnowballC::wordStem(x, "porter"),
                  ngram_window = c(1, 10),   # minimum and maximum n-gram length (!here setting lengths varied from one to ten tokens)
                  stopword_vec = stopw, 
                                             # stopword_vec = FALSE,
                  lower = TRUE,              # lowercase - this is the default value
                  remove_punctuation = TRUE, # convert all non-alpha numeric characters to spaces
                  remove_numbers = TRUE,     # onvert all numbers to spaces
                  verbose = FALSE,           # Turn off status bar for this demo
                  cpus = 2)                  # default is all available cpus on the system

 
 # construct the matrix of term counts to get the IDF vector
 tf_mat <- TermDocFreq(dtm)                  # term + term_freq + doc_freq + idf
 str(tf_mat)
 # head(tf_mat[ order(tf_mat$term_freq, decreasing = TRUE), ], 10)

 
 # look at the most frequent bigrams
 tf_bigrams <- tf_mat[ stringr::str_detect(tf_mat$term, "_") , ]
# knitr::kable(head(tf_bigrams[ order(tf_bigrams$term_freq, decreasing = TRUE) , ], 10),
#              caption = "Ten most frequent bi-grams") # look at the most frequent bigrams

 
# remove any tokens that were in 3 or fewer documents (option)
# dtm <- dtm[ , colSums(dtm > 0) > 3 ] # alternatively: dtm[ , tf_mat$term_freq > 3 ]
# tf_mat <- tf_mat[ tf_mat$term %in% colnames(dtm) , ]
# tf_bigrams <- tf_bigrams[ tf_bigrams$term %in% colnames(dtm) , ]
 
 
 # TF-IDF (several ways to calculate, need to check)
 tfidf <- t(dtm[ , tf_mat$term ]) * tf_mat$idf
 tfidf <- t(tfidf)

 
 # calculate cosine similarity / each pair and change it to a distance
 csim  <- tfidf / sqrt(rowSums(tfidf * tfidf))
 csim  <- csim %*% t(csim)
 tem <- as.data.frame(csim[1,])
 cdist <- 1 - tem
 #cdist <- as.dist(1 - csim)

 
 # Euclidean distance (need to change)
 d <- dist(cdist, method = 'euclidean')

 
 # find proximity hpo
 res<-as.matrix(d)
 res<-as.data.frame(res)
 res<-res[which(colnames(res)==sam),]
 res<-res[order(res[1,])]

 print(paste0(sam," is near ", colnames(res)[2]," with euclidean distance of ",res[1,2]))
 res<-as.data.frame(t(res[,2:b]))
 res['hpo_name']<-rownames(res)
 Id<-toString(sam)
 clinical <- toString(clisub[i,4])
 
 res_clisub[[clinical]] <- res
 rm(res)
 
}

tem <- data.frame(matrix(unlist(res_clisub), nrow=length(res_clisub), byrow=T))
tem <- cbind(names(res_clisub),tem)
res <- merge(clisub, tem, by.x = 'new_split', by.y = 'names(res_clisub)')

for (a in c(23:40)){
  for (b in 1:nrow(res)){
    content  = res[b,a] 
    content2 = toString(Df_cn[which(Df_cn$name==content),c(1,3)])
    content3 = paste(Df_cn[which(Df_cn$name==content),1],content,sep=', ')
    print(content2)
    if (content2!='character(0), character(0)') { res[b,a] = content2} else { res[b,a] = content3}
  }
}




for (a in c(5:22)){
  b = a + 18
  c = a + 36
  res[,c] = paste(res[,a], res[,b],sep=', ')}
res = res[,c(1:4,41:58)]
res <- transpose(res)


#rm(tem,dtm,tem,tf_bigrams,tf_mat,tfidf,df_subclass,res_clisub,a,b,cdist,d,i,Id,l,clinical,res_all, sam, csim, df)
write.csv(res,paste0('C:/Users/Jingz/OneDrive/CarrierGene/2019.07.10-HPO (done)/Cluster_R/analysis_output/','results_',p,'.csv'),row.names = FALSE)

}














# save into Rdata
# save(Res, file="clinical_results.RData")




# save into csv file
res_all<-res_all$全外


tem<-data.frame(1,2,3,4)
tem<-tem[-1,]
colnames(tem)<-c('distance','hpo','samid','class')
#RES<-data.frame(samid=1)


# for clutering in R

for (i in 2:length(res_all)){
  
  for (b in 1:length(res_all[[i]])){
  res<-as.data.frame(res_all[[i]][[b]])[1,]
  res$samid<-colnames(res)[1]
  colnames(res)<-c('distance','hpo','samid')
  res$class<-names(res_all)[i]
  Res<-rbind(Res,res)
  Res<-unique(Res)}
  
   if (i==1) {RES<-Res;rownames(RES)=1:nrow(RES)}
   if (i>1)  {RES<-merge(RES,Res,by='samid',all=TRUE)}
  
}

Res <- Res[c('samid','class','hpo','distance')]

RES <- merge(cli, Res[,1:10],by.x='Sample_ID', by.y='V1', all.y=TRUE)
write.csv(RES,'clinical_informations_all_results.csv',row.names = FALSE)






# transposed for trainning data (python)

data <- read_excel("C:/Users/Jingz/OneDrive/CarrierGene/2019.07.10-HPO (done)/NLP_NLP_python/trainning_data.xlsx")

library(dplyr)

for (i in 1:length(res_all)){
  
  res<-as.data.frame(res_all[[i]])
  res<-t(t(res)[1,])
  res[1,1]<-colnames(res)[1]
  colnames(res)[1]<-'Sample_id'
  for (b in 1:ncol(res)){if (b>1){res[,b]<-1}}
  
  if (i==1){Res<-res}
  if (i>1) {Res<-bind_rows(as.data.frame(Res), as.data.frame(res))}
}


colnames(Res)[2:ncol(Res)]<-paste0('Hp',1:(ncol(Res)-1))

cli <- read.csv("C:/Users/Jingz/OneDrive/CarrierGene/2019.07.10-HPO (done)/Hpo_cluster_R/clinical_data/clinical_informations.csv")
data<-merge(cli[,c(1,3,4)], Res, by.x='Sample_ID', by.y='Sample_id')
rownames(data)<-data$Sample_ID
write.table(data,'training_data.txt',sep='\t', na='0', row.names = FALSE, fileEncoding = 'utf-8')





# plot 
hc1<- df %>%
  dist(method = 'euclidean') %>% 
  hclust(method = 'average') 

dend1 <- as.dendrogram(hc1)
plot_dendro(dend1, height = 30000,width =3500,margin = 0.05) %>% 
  hide_legend() %>% 
  highlight(persistent = FALSE, dynamic = TRUE)







##################    Preparation before hc  ##################

### assessing clustering tendency

gradient.color <- list(low = "steelblue",  high = "white")


df %>%    # Remove column 5 (Species)
get_clust_tendency(n = 30, gradient = gradient.color)


### determining the optimal number of clusters

set.seed(123)
res.nbclust <- USArrests %>%
               scale() %>%
               NbClust(distance = "euclidean",
               min.nc = 2, max.nc = 10, 
               method = "complete", index ="all") 

fviz_nbclust(res.nbclust, ggtheme = theme_minimal())


### clustering validation statistics

set.seed(120)


### Model Based Clustering

library(mclust)
fit <- Mclust(df)
plot(fit) # plot results 
summary(fit) # display the best model










##################    Hierarchical clustering  ##################


####   import, remove NA and scale (for numeric data) 

df <- USArrests # Rows are observations and columns are variables
rownames(df)<-df$id
df<- df[,c(2:4)]
df <- na.omit(df)

#normalize to obtain homoscedastic data : the variance of an observable quantity does not depend on the mean
df<-scale(df)
head(df)




####   verify which method is better 

# compute correlation between the cophenetic distances and the original distance data generated by the dist()
# the closer the value of the correlation coefficient is to 1, the more accurately the clustering solution reflects data.
# Values above 0.75 are felt to be good

m <- c( "average", "single", "complete", "ward") # methods to assess
names(m) <- c( "average", "single", "complete", "ward")


for (i in names(m)){
  ac <- hclust(d, method = i)
  print(paste0(i,': ', cor(d, cophenetic(ac)))) # Correlation between cophenetic distance and the original distance
}


ac <- function(x) {agnes(df, method = x)$ac} # verify the cluster tree via agnes
map_dbl(m, ac) # ward has highest value





#### algorithm with optimal method

d<-dist(df, method = 'ward')

hc1<- df %>%
      dist(method = 'euclidean') %>%  # Dissimilarity matrix : "Manhattan", "correlation", "maximum"
      hclust(method = 'average')      # Linkage clustering. "ward.D", "ward.D2", "single", "complete", "average", "centroid" or "mcquitty", "median"
plot(hc1) 

hc2 <- agnes(df, metric = "euclidean", method = "ward", stand = TRUE) # compute via agnes (alternative)
hc2$ac # get the agglomerative coefficient, which measures the amount of clustering structure found. 
       # values closer to 1 suggest strong clustering structure


hc3 <- diana(df, metric = "euclidean") # compute divisive hierarchical clustering
hc3$dc # Divise coefficient; amount of clustering structure found








####   create a summary table of the top 5 words defining each cluster   

clustering <- cutree(hc1, 30)
p_words <- colSums(dtm) / sum(dtm)

cluster_words <- lapply(unique(clustering), function(x){
  rows <- dtm[ clustering == x , ]
  
  # for memory's sake, drop all words that don't appear in the cluster
  rows <- rows[ , colSums(rows) > 0 ]
  
  colSums(rows) / sum(rows) - p_words[ colnames(rows) ]
})



cluster_summary <- data.frame(cluster = unique(clustering),
                              size = as.numeric(table(clustering)),
                              top_words = sapply(cluster_words, function(d){
                                paste(
                                  names(d)[ order(d, decreasing = TRUE) ][ 1:5 ], 
                                  collapse = ", ")
                              }),
                              stringsAsFactors = FALSE)

cluster_summary
knitr::kable(cluster_summary, caption = "Cluster summary table")






###   plot a word cloud of one cluster as an example  

wordcloud::wordcloud(words = names(cluster_words[[ 5 ]]), 
                     freq = cluster_words[[ 5 ]], 
                     max.words = 50, 
                     random.order = FALSE, 
                     colors = c("red", "yellow", "blue"),
                     main = "Top words in cluster 100")







###   plot Dendrograms 


k=4
sub_grp <- cutree(hc1, k) # Cut tree into groups, for hclust
           # cutree(as.hclust(hc2), k) # cut tree for diana and agnes
sub_grp
table(sub_grp) # Number of members in each cluster
rownames(df)[sub_grp == 2]


# HCLUST plot
plot(hc1, cex = 0.6, hang = -1)
rect.hclust(hc1, k , border = 2:(k+1))


# AGNES plot using plot.hclust() + plot.dendrogram() + pltree()
plot(as.hclust(hc2), cex = 0.6, hang = -1)
plot(as.dendrogram(hc2), cex = 0.6, hang = -1)


# DIANA plot using pltree()
pltree(hc3,  main = "Dendrogram of diana")


# ALL plot
fviz_dend(hc1, 20, # Cut in four groups
          cex = 0.5, # label size
          k_colors = rainbow(20),
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE, # Add rectangle around groups
          palette = "jco" # Color palette
)


# plot cloud
fviz_cluster(list(data = df, cluster = sub_grp)) 



# interactive dendrogram via plotly
library(ggplot2)
library(ggdendro)
library(plotly)

dend1 <- as.dendrogram(hc1)
plot_dendro(dend1, height = 10000,width =3500,margin = 0.05) %>% 
  hide_legend() %>% 
  highlight(persistent = FALSE, dynamic = TRUE)



library(idendro)
idendro(hc1, tem)

library(collapsibleTree) # cute plot
library(networkD3)
radialNetwork(as.radialNetwork(hc1))

library(idendr0)





# ape plot
library(ape)

plot(as.phylo(hc1), cex = 0.6, label.offset = 0.5)
plot(as.phylo(hc1), type = "unrooted", cex = 0.6, no.margin = TRUE) # Unrooted
plot(as.phylo(hc1), type = "fan")

colors = rainbow(20)
clus4 = cutree(hc, 20)
plot(as.phylo(hc), type = "fan", tip.color = colors[clus4],
     label.offset = 1, cex = 0.7)





###  plot correlation (necessary for overlapped cluster: one HPO term appear in two different clusters)
























###    compare two dendrograms    

# Compute distance matrix
d <- dist(df, method = "euclidean")

# Compute 2 hierarchical clusterings
set.seed(100)
hc1 <- hclust(d, method = "complete")
set.seed(200)
hc2 <- hclust(d, method = "ward.D2")

# Create two dendrograms
dend1 <- as.dendrogram (hc1)
dend2 <- as.dendrogram (hc2)

tanglegram(dend1, dend2) #  unique nodes highlighted with dashed lines
                         #  Entanglement is a measure between 1 (full entanglement) and 0 (no entanglement)
                         #  A lower entanglement coefficient corresponds to a good alignment


# customized tanglegram
dend_list <- dendlist(dend1, dend2)
tanglegram(dend1, dend2,
           highlight_distinct_edges = FALSE,      # Turn-off dashed lines
           common_subtrees_color_lines = FALSE,   # Turn-off line colors
           common_subtrees_color_branches = TRUE, # Color common branches 
           main = paste("entanglement =", round(entanglement(dend_list), 2))
           )



cor.dendlist(dend_list, method = "cophenetic")    # Cophenetic correlation matrix
cor.dendlist(dend_list, method = "baker")         # Baker correlation matrix
cor_cophenetic(dend1, dend2)                      # Cophenetic correlation coefficient
cor_bakers_gamma(dend1, dend2)                    # Baker correlation coefficient

# http://www.sthda.com/english/wiki/print.php?id=237










###   Determining Optimal Clusters 

# Elbow Method
fviz_nbclust(df, FUN = hcut, method = "wss")

# Average Silhouette Method
fviz_nbclust(df, FUN = hcut, method = "silhouette")

# Gap Statistic Method
gap_stat <- clusGap(df, FUN = hcut, nstart = 25, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)




















##################    Semi - Hierarchical clustering  ##################
# http://bgmm.molgen.mpg.de/rapBGMM/
# available packages : hopach; hierinf; BHC; cola.


### method BHC (does not work with tfidf data)

library(BHC)

nDataItems <- nrow(dtm)
nFeatures <- ncol(dtm)
itemLabel <- tem$label

hc1 <- bhc(df, itemLabel, verbose=TRUE)
plot(hc1, axes=FALSE)
WriteOutClusterLabels(hc1, "labels.txt", verbose=TRUE)


newData <- data[] + rnorm(150, 0, 0.1);
percentiles <- FindOptimalBinning(newData, itemLabels, transposeData=TRUE, verbose=TRUE)
discreteData <- DiscretiseData(t(newData), percentiles=percentiles)
discreteData <- t(discreteData)
hc3 <- bhc(discreteData, itemLabels, verbose=TRUE)
plot(hc3, axes=FALSE)



















##################    K-Mean clustering (not suitable for HPO, but ok for HPO subgroup)  ##################



### data prepare (import + remove NA + scale variables)

# Load  and prepare the data
df<-matrix(sample(1:2),ncol = 2)
df<-as.data.frame(df)
colnames(df)<-c('x','label')
df<-df[-1,]

setwd("C:/Users/Jingz/OneDrive/CarrierGene/2019.07.10-HPO/Hpo_ML/Hpo_analysis_R")
a<-1:length(list.files())

for (l in a){
  setwd("C:/Users/Jingz/OneDrive/CarrierGene/2019.07.10-HPO/Hpo_ML/Hpo_analysis_R")
  Label<-toString(list.files()[l])
  PATH<-toString(list.files()[a[l]])  # read HPO terms related to each panel type
  setwd(PATH)
  set.seed(300)
  descendants <-read.csv(list.files()[1], sep="")
  descendants<-na.omit(descendants)
  descendants$label<-Label
  df<-rbind(df,descendants)
}

df<-df[unique(df$x),]

df_subclass<-merge(df[,],Df_en, by.x='x',by.y='id',all.x = TRUE)  # merge with hpo database
rm(df)
df_subclass$full_text<-paste(df_subclass$name,df_subclass$synonyms,df_subclass$description,sep='.')
df_subclass<-df_subclass[c('x','name','full_text')] 

write.csv(df_subclass,'Hpo_merged.csv')


tem<-df_subclass


# create a document term matrix (tokenizing)
stopw<- c(stopwords::stopwords("en"), # stopwords from tm
          stopwords::stopwords(source = "smart")) # this is the default value

stopw<- stopw[!grepl('^no$|^non$|^not$|^none$|^noone$',stopw)]

dtm <- CreateDtm(doc_vec =tem$full_text, # character vector of documents
                 doc_names = tem$name, # document names
                 stem_lemma_function = function(x) SnowballC::wordStem(x, "porter"),
                 ngram_window = c(1, 2), # minimum and maximum n-gram length
                 stopword_vec = stopw, 
                 lower = TRUE, # lowercase - this is the default value
                 remove_punctuation = FALSE, # convert all non-alpha numeric characters to spaces
                 remove_numbers = FALSE, # onvert all numbers to spaces
                 verbose = FALSE, # Turn off status bar for this demo
                 cpus = 2) # default is all available cpus on the system

# construct the matrix of term counts to get the IDF vector
tf_mat <- TermDocFreq(dtm) # term + term_freq + doc_freq + idf
#str(tf_mat)

#tf_bigrams <- tf_mat[ stringr::str_detect(tf_mat$term, "_") , ]
#knitr::kable(head(tf_bigrams[ order(tf_bigrams$term_freq, decreasing = TRUE) , ], 10),
#             caption = "Ten most frequent bi-grams") # look at the most frequent bigrams

# TF-IDF (several ways to calculate, need to check)
tfidf <- t(dtm[ , tf_mat$term ]) * tf_mat$idf
tfidf <- t(tfidf)

#  calculate cosine similarity / each pair and change it to a distance
csim <- tfidf / sqrt(rowSums(tfidf * tfidf))
csim <- csim %*% t(csim)
cdist <- as.dist(1 - csim)


# Euclidean distance
d<-dist(cdist, method = 'euclidean') 














df <- USArrests %>%
  na.omit() %>%          # Remove missing values (NA)
  scale()                # Scale variables

# View the firt 3 rows
head(df, n = 3)





### distance measures

# get_dist(): for computing a distance matrix between the rows of a data matrix. 
# Compared to the standard dist() function, it supports correlation-based distance measures including ???pearson???, ???kendall??? and ???spearman??? methods.
# fviz_dist(): for visualizing a distance matrix


res.dist <- get_dist(df, stand = TRUE, method = "pearson")
fviz_dist(res.dist, 
          gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) # heatmap - group correlation




### Partitioning clustering ( K-means clustering + K-medoids clustering/ Partitioning Around Medoids / PAM )


# method1 : K-Mean
# Determining the optimal number of clusters: use factoextra::fviz_nbclust()
fviz_nbclust(my_data, kmeans, method = "gap_stat")

# Compute and visualize k-means clustering
set.seed(123)
km.res <- kmeans(df, 5, nstart = 25)

# Visualize
fviz_cluster(km.res, data = df,
             ellipse.type = "convex",
             palette = "jco",
             ggtheme = theme_minimal())


# method 2 : PAM 
pam.res <- pam(df, 5)
# Visualize
fviz_cluster(pam.res)





