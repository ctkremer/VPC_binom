
###############################

## Use MCMCglmm to fit hierarchical binomial model with only random effects terms, and generate variance partition coefficients.

# Accompanies Royer et al. 2016, "Incomplete Loss of a Conserved Trait: Function, Latitudinal Cline, and Genetic Constraints", Evolution.

# Developed based on Browne et al. 2005, and Nakagawa & Schielzeth 2010

# Code written by C.T. Kremer

###############################

## Load tools:
library(MCMCglmm)
library(reshape2)
library(ggplot2)
library(dplyr)

expit<-function(x){exp(x)/(1+exp(x))}
logit<-function(x){log(x/(1-x))}

###############################

### Load data set:

###############################

### RIL data:

dat2<-read.csv("2010RILsFlowerAsRow_NoParents.csv",stringsAsFactors =F)

# organize data:
dat2$line<-as.factor(dat2$Line)
dat2$plant<-as.factor(dat2$Plant.nested.in.line)
dat2$short<-dat2$shorts

dat2<-dat2[,c('line','plant','short')]
dat2$id<-paste(dat2$line,dat2$plant)

dat2$stamens.lost <- 2 - dat2$short
dat2$stamens.kept <- dat2$short

uids<-unique(dat2$id)
flwr<-c()
for(i in 1:length(uids)){
	len<-length(dat2$id[dat2$id==uids[i]])
	flwr<-append(flwr,seq(1,len,1))
}
dat2$flower<-flwr


##################

## Fit the hierarchical model:

# specify priors
prior2<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002),G3=list(V=1,nu=0.002)))

#### Fit model 
#	- CAUTION: can take a day or more to run.
#	- run commented out code below if you actually intend to fit the MCMCglmm
#	- otherwise skip ahead to VPC plotting

#stamens.lost2<-dat2$stamens.lost
#stamens.kept2<-dat2$stamens.kept

#m1<- MCMCglmm(cbind(stamens.lost2,stamens.kept2)~1,
#	random=~line+plant:line+flower:plant:line, family="multinomial2",
#	data=dat2,verbose=T,
#	nitt=2*10^7+3*10^4,thin=5*10^3,burnin=3*10^4,prior=prior2,pr=F)	

plot(m1)
summary(m1)

#save(m1,file="MCMCglmm_RIL_chains_051616.Rdata")

####

 # Iterations = 30001:20025001
 # Thinning interval  = 5000
 # Sample size  = 4000 

 # DIC: 10349.26 

 # G-structure:  ~line

     # post.mean l-95% CI u-95% CI eff.samp
# line     2.159     1.74    2.622     1759

               # ~plant:line

           # post.mean l-95% CI u-95% CI eff.samp
# plant:line     1.358    1.115    1.619     1317

               # ~flower:plant:line

                  # post.mean  l-95% CI u-95% CI eff.samp
# flower:plant:line   0.03777 0.0001982   0.1485    431.6

 # R-structure:  ~units

      # post.mean l-95% CI u-95% CI eff.samp
# units   0.04356 0.000284   0.1618    927.1

 # Location effects: cbind(stamens.lost2, stamens.kept2) ~ 1 

            # post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)    -2.818   -2.993   -2.646     1408 <3e-04 ***
# ---
# Signif. codes:  
# 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

######

# Variance estimates:
#vars<-c(summary(m1)$Gcovariances[,1],summary(m1)$Rcovariances[1])
vars<-c(2.15945808,1.35814703,0.03777143,0.04355962)

# Intercept estimate:
#intr<-summary(m1)$solutions[,1]
intr<- c(0.196425,-2.818469,-5.449155)

# Intercepts = ~1.88, to match RIL mean... and use Italian and Swedish means for upper and lower bounds of of RIL plot.


### Set up variance partition calculation

get.vpc<-function(b0){
	
	P<-expit(b0)
	v.line<-vars[1]
	v.plant.line<-vars[2]
	v.flower.plant.line<-vars[3]
	v.e<-vars[4]

	scale<- (P^2)/(1+exp(b0))^2
	denom<- (v.line+v.plant.line+v.flower.plant.line+v.e)*scale+P*(1-P)
	
	p.v1<- (v.line*scale)/denom
	p.v2<- (v.plant.line*scale)/denom
	p.v3<- (v.flower.plant.line*scale)/denom
	p.v.e<- (v.e*scale)/denom
	res<-c(p.v1,p.v2,p.v3,p.v.e,P*(1-P)/denom,denom)
	names(res)<-c('line','plant','flower','error','distribution','total')
	return(res)
}


# Target mean short stamen number:
targets<-intr

# Try this out at each intercept position 
r2<-sapply(targets,get.vpc)

# Make a clean figure:
colnames(r2)<-c(paste("np = ",round(2-2*expit(intr[1]),2),sep=""),
				paste("np = ",round(2-2*expit(intr[2]),2),sep=""),
				paste("np = ",round(2-2*expit(intr[3]),2),sep=""))
r2<-as.data.frame(r2)
r2$component<-rownames(r2)
r3<-melt(r2)
names(r3)<-c("component","mean.short.stamen","value")
r3$comp<-as.numeric(r3$mean.short.stamen)

tmp<-r3[r3$component=='total',c('mean.short.stamen','value')]
names(tmp)[which(names(tmp)=="value")]<-"tot"

r4<-merge(r3,tmp)
r4$value2<-r4$value*r4$tot
r4<-r4[r4$component!="total",]
r4$component<-factor(r4$component,levels=rev(c('distribution','error','line','plant','flower')))
r4<-r4[order(r4$mean.short.stamen,r4$component),]

r4.ril<-r4
r4.ril$data.set<-"RIL"


### Plot results:

p2<-ggplot(r4.ril,aes(x=mean.short.stamen,y=value2,fill=component))+
	geom_bar(stat="identity")+
	scale_y_continuous("Total contribution")+
	scale_x_discrete("Mean short stamen number")+
	scale_fill_manual('Variance\ncomponent',
		values=c('#d7191c','#fdae61','#ffffbf','#abd9e9','#2c7bb6'))+
	guides(fill=guide_legend(reverse=T))+
	theme_bw()+
	ggtitle("RIL parent data")
p2

ggsave("VPC_RIL_PARENT_raw_051616.pdf",p2,width=5.5,height=5.5,units=c('in'))

####


# helper function for testing if there's variation in stamen number:
fl.mn<-function(x){
	(mean(x) %% 1) == 0
}

# % of flowers in RIL data have no variation
tmp3<-dat2 %>% group_by(line,plant,flower) %>% summarise(flg=fl.mn(stamens.kept))
func<-function(x) sum(x)/(length(x))
tmp3 %>% group_by() %>% summarise(func(flg))
# 0.9991273

# % of plants in RIL data have no variation
tmp3<-dat2 %>% group_by(line,plant) %>% summarise(flg=fl.mn(stamens.kept))
func<-function(x) sum(x)/(length(x))
tmp3 %>% group_by() %>% summarise(func(flg))
# 0.6731634

# % of lines in RIL data have no variation
tmp3<-dat2 %>% group_by(line) %>% summarise(flg=fl.mn(stamens.kept))
func<-function(x) sum(x)/(length(x))
tmp3 %>% group_by() %>% summarise(func(flg))
# 0.2640845

######################

## RIL parents analysis


dat2<-read.csv("RILparentsForColin_StamenNum.csv",stringsAsFactors =F)

# organize data:
dat2$population<-as.factor(dat2$Pop)
dat2$line<-as.factor(dat2$Line)
dat2$plant<-as.factor(dat2$Rep)
dat2$short<-dat2$ShortStamen

dat2<-dat2[,c('population','line','plant','short')]
dat2$id<-paste(dat2$population,dat2$line,dat2$plant)

dat2$stamens.lost <- 2 - dat2$short
dat2$stamens.kept <- dat2$short

uids<-unique(dat2$id)
flwr<-c()
for(i in 1:length(uids)){
	len<-length(dat2$id[dat2$id==uids[i]])
	flwr<-append(flwr,seq(1,len,1))
}
dat2$flower<-flwr


##################

## Fit the hierarchical model:

# specify priors
prior2<-list(R=list(V=1,nu=0.002),
			G=list(G1=list(V=1,nu=0.002),
					G2=list(V=1,nu=0.002),
					G3=list(V=1,nu=0.002)))

# fit model
#	- CAUTION! Takes several days on a fast computer to run ths model.
#	- run commented out code below if you actually intend to fit the MCMCglmm
#	- otherwise skip ahead to VPC plotting

#stamens.lost2<-dat2$stamens.lost
#stamens.kept2<-dat2$stamens.kept

#m1<- MCMCglmm(cbind(stamens.lost2,stamens.kept2)~1+population,
#	random=~line:population+plant:line:population+flower:plant:line:population, 
#	family="multinomial2",
#	data=dat2,verbose=T,
#	nitt=6*10^7+3*10^4,thin=6*10^3,burnin=3*10^4,
#	prior=prior2,pr=F)
#)

plot(m1)
summary(m1)

# Save output:
#save(m1,file="MCMCglmm_RIL_PARENT_chains_051216.Rdata")


#################

 # Iterations = 30001:60024001
 # Thinning interval  = 6000
 # Sample size  = 10000 

 # DIC: 1223.158 

 # G-structure:  ~line:population

                # post.mean  l-95% CI u-95% CI eff.samp
# line:population    0.1935 0.0001976   0.6385     9370

               # ~plant:line:population

                      # post.mean l-95% CI u-95% CI eff.samp
# plant:line:population     2.522     1.56    3.638     5921

               # ~flower:plant:line:population

                             # post.mean  l-95% CI u-95% CI
# flower:plant:line:population   0.04094 0.0001483   0.1764
                             # eff.samp
# flower:plant:line:population     7587

 # R-structure:  ~units

      # post.mean  l-95% CI u-95% CI eff.samp
# units   0.04152 0.0001739   0.1796     9426

 # Location effects: cbind(stamens.lost2, stamens.kept2) ~ 1 + population 

             # post.mean l-95% CI u-95% CI eff.samp  pMCMC
# (Intercept)     0.1964  -0.1959   0.6265    10000  0.317
# populationSW   -5.4492  -6.6298  -4.3733     4395 <1e-04
                
# (Intercept)     
# populationSW ***
# ---
# Signif. codes:  
# 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Variance estimates:
#vars<-c(summary(m1)$Gcovariances[,1],summary(m1)$Rcovariances[1])
vars<-c(0.19348549,2.52234258,0.04093539,0.04152305)

# Intercept estimate:
#intr<-summary(m1)$solutions[,1]
intr<- c(0.196425,-2.818469,-5.449155)

# Intercepts = ~1.88, to match RIL mean... and use Italian and Swedish means for upper and lower bounds of of RIL plot.



####

# Now use this for variance partitioning.

get.vpc<-function(b0){
	
	P<-expit(b0)
	v.line<-vars[1]
	v.plant.line<-vars[2]
	v.flower.plant.line<-vars[3]
	v.e<-vars[4]

	scale<- (P^2)/(1+exp(b0))^2
	denom<- (v.line+v.plant.line+v.flower.plant.line+v.e)*scale+P*(1-P)
	
	p.v1<- (v.line*scale)/denom
	p.v2<- (v.plant.line*scale)/denom
	p.v3<- (v.flower.plant.line*scale)/denom
	p.v.e<- (v.e*scale)/denom
	res<-c(p.v1,p.v2,p.v3,p.v.e,P*(1-P)/denom,denom)
	names(res)<-c('line','plant','flower','error','distribution','total')
	return(res)
}

# Target mean short stamen number:
targets<-intr

# Try this out at each intercept position 
r2<-sapply(targets,get.vpc)

# Make a clean figure:
colnames(r2)<-c(paste("np = ",round(2-2*expit(intr[1]),2),sep=""),
				paste("np = ",round(2-2*expit(intr[2]),2),sep=""),
				paste("np = ",round(2-2*expit(intr[3]),2),sep=""))
r2<-as.data.frame(r2)
r2$component<-rownames(r2)
r3<-melt(r2)
names(r3)<-c("component","mean.short.stamen","value")
r3$comp<-as.numeric(r3$mean.short.stamen)

tmp<-r3[r3$component=='total',c('mean.short.stamen','value')]
names(tmp)[which(names(tmp)=="value")]<-"tot"

r4<-merge(r3,tmp)
r4$value2<-r4$value*r4$tot
r4<-r4[r4$component!="total",]
r4$component<-factor(r4$component,levels=rev(c('distribution','error','line','plant','flower')))
r4<-r4[order(r4$mean.short.stamen,r4$component),]

# Save results:
r4.ril.p<-r4
r4.ril.p$data.set<-"RIL.Parent"


## Plot the result:
p2<-ggplot(r4,aes(x=mean.short.stamen,y=value2,fill=component))+
	geom_bar(stat="identity")+
	scale_y_continuous("Total contribution")+
	scale_x_discrete("Mean short stamen number")+
	scale_fill_manual('Variance\ncomponent',
		values=c('#d7191c','#fdae61','#ffffbf','#abd9e9','#2c7bb6'))+
	guides(fill=guide_legend(reverse=T))+
	theme_bw()+
	ggtitle("RIL parent data")
p2

ggsave("VPC_RIL_PARENT_raw_051216.pdf",p2,width=5.5,height=5.5,units=c('in'))


################

# helper function for testing if there's variation in stamen number:
fl.mn<-function(x){
	(mean(x) %% 1) == 0
}

# % of plants in RIL Parent data have no variation
tmp3<-dat2 %>% group_by(population,line,plant,flower) %>% summarise(flg=fl.mn(stamens.kept))
func<-function(x) sum(x)/(length(x))
tmp3 %>% group_by(population) %>% summarise(func(flg))

# % of plants in RIL Parent data have no variation
tmp3<-dat2 %>% group_by(population,line,plant) %>% summarise(flg=fl.mn(stamens.kept))
func<-function(x) sum(x)/(length(x))
tmp3 %>% group_by(population) %>% summarise(func(flg))

# % of lines in RIL Parent data have no variation
tmp3<-dat2 %>% group_by(population,line) %>% summarise(flg=fl.mn(stamens.kept))
func<-function(x) sum(x)/(length(x))
tmp3 %>% group_by(population) %>% summarise(func(flg))



################################


### Merged graphic across RILs and RIL parents:

# combine data sets
r4.j<-rbind(r4.ril.p,r4.ril)
r4.j$data.set<-factor(r4.j$data.set,levels=c("RIL.Parent","RIL"))

p3<-ggplot(r4.j,aes(x=mean.short.stamen,y=value2,fill=component))+
	geom_bar(stat="identity")+
	facet_wrap(~data.set)+
	scale_y_continuous("Total contribution")+
	scale_x_discrete("\nMean short stamen number")+
	scale_fill_manual('Variance\ncomponent',
		values=c('#d7191c','#fdae61','#ffffbf','#abd9e9','#2c7bb6'))+
	guides(fill=guide_legend(reverse=T))+
	theme_bw()+
	ggtitle("VPC Analysis")
p3

ggsave("VPC_full_raw_051616.pdf",p3,width=9,height=5.5,units=c('in'))


