  #-------------------------------------------------------------

### Data exploration of Ptarmigan gut microbiome

#-------------------------------------------------------------

# Load libraries
library('topicmodels')
library('reshape2')
library('dplyr')
library('bbmle')
library('ggplot2')
library('tidytext')
library('data.table')

#-------------------------------------------------------------

#load genus level microbiom data
ptarmiganGenus <- read.csv("data/ptarmigan_taxonomy_genus.csv") # genus level resolution of microbiome
head(ptarmiganGenus)

#convert from long to wide format
ptarmiganGenus2 <- dcast(ptarmiganGenus, sample~genus, value.var="abundance")

#query Stephanie on sample id's
unique(ptarmiganGenus2$sample)

#IMPORTANT TO CLIP RECORDS TO EXCLUDE NON-MICROBE GENETIC INFORMATION
##format data
ptarmiganGenus3 <- ptarmiganGenus2[,-1] #drop sample id
rownames(ptarmiganGenus3) <- ptarmiganGenus2[,1]

#drop rows/sites that sum to zero = no fish recorded
ptarmiganGenus3$sum <- rowSums(ptarmiganGenus3)
ptarmiganGenus3 <- ptarmiganGenus3 %>%
  filter(sum>0) %>%
  select(-c(sum))

VEM2=LDA(ptarmiganGenus3,k=2,control=list(seed=240444)) #fit model with 2 communities
VEM3=LDA(ptarmiganGenus3,k=3,control=list(seed=240444)) #fit model with 3 communities
VEM4=LDA(ptarmiganGenus3,k=4,control=list(seed=240444)) #fit model with 4 communities
VEM5=LDA(ptarmiganGenus3,k=5,control=list(seed=240444)) #fit model with 5 communities
VEM6=LDA(ptarmiganGenus3,k=6,control=list(seed=240444)) #fit model with 6 communities
VEM7=LDA(ptarmiganGenus3,k=7,control=list(seed=240444)) #fit model with 6 communities


# model selection for the optimal # communities identified
AICtab(VEM3,VEM2,VEM4,VEM5,VEM6,VEM7)
# 5 communities results in the best model fit

#get parameter estimates
z=posterior(VEM5)
commun.plot= as.data.frame(z$topics) #community proportions in each sample bird
str(commun.plot)
commun.plot$samples <-ptarmiganGenus2$sample



commun.spp=z$terms #probability of a species belonging to a community
# look at the specific assignment of species to a community
ap_lda_td <- data.table(tidy(VEM5))
nrow(ap_lda_td[topic==1])


ap_top_terms <- ap_lda_td %>%
  group_by(topic) %>%
  top_n(10, beta) %>%
  ungroup() %>%
  arrange(topic, -beta)

ap_top_terms %>%
  mutate(term = reorder(term, beta)) %>%
  ggplot(aes(term, beta, fill = factor(topic))) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  facet_wrap(~ topic, scales = "free") +
  coord_flip() 


ap_lda_td[topic==4][order(beta, decreasing=TRUE)][1:10] #beta is the prob. of a specific sp. belonging to community 1
#modify the "topic==4" to see the different membership of microbiome to communities




#-------------------------------------------------------------

# look at the proportion of the communities in each unique sample

ap_gamma <- data.table(tidy(VEM5, matrix = "gamma"))
ap_gamma[document==16]
ap_gamma$birdID <- rep(ptarmiganGenus2$sample,5)

ap_gamma %>%
  ggplot(aes(topic, gamma, fill = factor(birdID))) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  facet_wrap(~ birdID, scales = "free") +
  coord_flip() 



