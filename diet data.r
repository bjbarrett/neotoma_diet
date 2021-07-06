
library(rethinking)
library(rstan)
library(reshape)
data <- read.csv("/home/bjbarrettadmin/Dropbox/Quail Ridge Woodrat Data/Diet Study/arielle_feeding.csv", header=TRUE)
data <- read.csv("~/Dropbox/Quail Ridge Woodrat Data/Diet Study/arielle_feeding.csv", header=TRUE, strip.white=TRUE)
d <- data
d <- reshape(data,varying=list(names(data)[6:12]),direction='long' )
names(d)[names(d) == 'time'] <- 'species'
d$species <- as.factor(d$species)
levels(d$species) <-  c("Poison" ,"Toyon","Bay","Buckeye","Cedar","chamise","Oak")
names(d)[names(d) == 'Poison'] <- 'grams_eaten'
d$log_distance <- log(d$Distance)
d$species_index <- as.integer(d$species)


geo <- read.csv("~/Dropbox/Quail Ridge Woodrat Data/WoodratHousesUTM.csv", header=TRUE , strip.white=TRUE)
house <- geo[1:23,]
chemline <- geo[24:40,]
house$dist2d <- 0
house$dist3d <- 0

for (i in 1:23){
house$dist2d[i] <- min(sqrt( (house$LatUTM[i] - chemline$LatUTM )^2 + (house$LongUTM[i] - chemline$LongUTM)^2 ))
house$dist3d[i] <- min(sqrt( (house$LatUTM[i] - chemline$LatUTM )^2 + (house$LongUTM[i] - chemline$LongUTM)^2  + (house$Elevation[i] - chemline$Elevation)^2 ))

}
house$Location <- house$Name
house <- subset(house, select=c(Location,dist2d,dist3d))

d<- merge(d,house, by="Location")
d0 <- merge(data,house, by="Location")
d$dist3d.c <- d$dist3d- mean(d$dist3d)
d$dist2d.c <- d$dist2d- mean(d$dist2d)
d$cell <- 0
for (i in 1:nrow(d)){
 d$cell[i] <- i
}
d$logdist2d <- log(d$dist2d)
d$logdist3d <- log(d$dist3d)
d$totalgramssum <- 0
 for (i in 1:max(d$id)){
    d$totalgramssum[d$id==i] <- sum(d$grams_eaten[d$id==i])
    
}

d$consumed <- ifelse(d$grams_eaten>0 , 1 , 0 )

#d <- subset(d,trial_indiv==1)

#m5521 wrong in arielle data i changed it
plot(totalgramssum~logdist2d, data=subset(d, Adult==0) , col="blue")
points(totalgramssum~logdist2d, data=subset(d, Adult==1) , col="red")

 par(mfrow=c(3, 3),oma = c(3, 3, 0, 0) , cex=0.7 , mar=c(1.1,1.1,1,1))
plant_key <- c("Poison Oak","Toyon" , "Bay" , "Buckeye" , "Cedar" , "chamise" , "Oak")
for (i in 1:7){
dens(d$grams_eaten[d$species_index==i], main=plant_key[i])
}
dens(d$grams_eaten, main="all species")

par(mfrow=c(4, 4),oma = c(3, 3, 0, 0) , cex=0.7 , mar=c(1.1,1.1,1,1))
plant_key <- c("Poison Oak","Toyon" , "Bay" , "Buckeye" , "Cedar" , "chamise" , "Oak")
for (i in 1:7){
    dens(d$grams_eaten[d$species_index==i & d$Adult==0], main=plant_key[i] , col="blue")
    lines(dens(d$grams_eaten[d$species_index==i & d$Adult==1], main=plant_key[i] , col="red"))
}
dens(d$grams_eaten, main="all species")

par(mfrow=c(3, 3),oma = c(3, 3, 0, 0) , cex=0.7 , mar=c(1.1,1.1,1,1))
plant_key <- c("Poison Oak","Toyon" , "Bay" , "Buckeye" , "Cedar" , "chamise" , "Oak")
for (i in 1:7){
    dens(d$grams_eaten[d$species_index==i & d$Adult==1], main=plant_key[i] , col="red")
}
dens(d$grams_eaten, main="all species")

###MODELS
m0 <- map2stan(
    alist(
        Distance ~ dgamma2(mu,scale),
        log(mu) ~ a + b_adult*adult,
        a ~ dnorm(0,2),
        b_adult ~ dnorm(0,2),
        scale ~ dexp(2)
    ),
    data=list(
        #Male=as.integer(d$Male),
        Distance=d0$dist2d,
        adult=as.integer(d0$Adult)
    ),
       cores=2 , warmup=4000 , iter=8000 , WAIC=TRUE, constraints=list(scale="lower=0")
)

m0l <- map2stan(
alist(
Distance ~ dgamma2(mu,scale),
log(mu) ~ a + b_adult*adult,
a ~ dnorm(0,2),
b_adult ~ dnorm(0,2),
scale ~ dexp(2)
),
data=list(
#Male=as.integer(d$Male),
Distance=d0$dist2d,
adult=as.integer(d0$Adult)
),
cores=2 , warmup=4000 , iter=8000 , WAIC=TRUE, constraints=list(scale="lower=0")

m1l <- map2stan(
alist(
Distance ~ dgamma2(mu,scale),
log(mu) ~ a + a_ratid + b_adult*adult,
a_ratid[rat_index] ~ dnorm(0,sigma_ratid),
a ~ dnorm(0,2),
b_adult ~ dnorm(0,2),
sigma_ratid ~ dcauchy(0,2),
scale ~ dexp(2)
),
data=list(
#Male=as.integer(d$Male),
Distance=d0$dist2d,
adult=as.integer(d0$Adult),
rat_index=as.integer(d0$Rat_ID)

),
cores=2 , chains =2 , warmup=4000 , iter=8000 , WAIC=TRUE, constraints=list(scale="lower=0")
)

png("DistaneDensDietPref.png",res=300,height=7,width=7, units = 'in')

post <- extract.samples(m0)
plot(density(d0$dist2d[d0$Adult==0], adjust=1.8), xlim=c(0,90)  , col="white" , ylim=c(-0.004,0.08) , main="" , xlab="Distance (m) from habitat edge")
#lines(density(d0$dist2d[d0$Adult==1]) )
#lines(density(d0$dist2d[d0$Adult==0] , adjust=1.8)  , lw=5 , lty=3, col="red" )
#lines(density(d0$dist2d[d0$Adult==1] , adjust=1.8)  , lw=5 , lty=3 , col="blue" )

for ( i in 1:100 ) {
    curve(dgamma(x , shape=post$a[i] , scale=post$scale[i])  , add=TRUE ,  col=col.alpha("red",0.1))
    curve(dgamma(x , shape=(post$a[i] + post$b_adult[i]) , scale=post$scale[i] ) , add=TRUE ,  col=col.alpha("blue",0.1) , lty=1)
}

curve(dgamma(x , shape=mean(post$a) , scale=median(post$scale))  , add=TRUE , col="red" , lw=5)
curve(dgamma(x , shape=mean((post$a + post$b_adult)) , scale=median(post$scale) ) , add=TRUE , col="blue" , lty=2 , lw=5)

points( d0$dist2d[d0$Adult==0] , (d0$Adult[d0$Adult==0]-0.0011) , pch="J" , col="red")
points( d0$dist2d[d0$Adult==1] , (d0$Adult[d0$Adult==1] -1.0034) , pch="A" , col="blue")
text(45,.075, expression(log(mu)%~% 2.92 + 0.35%*%Adult))
#text(45,.077, expression(scale==5.95))
text(45,.08, expression(dgamma(mu,scale==5.95)))

#
legend("topright", c("Adult", "Juvenile") , lty=c(2,1) , col=c("blue","red"), horiz=FALSE, bty='n' , lw=3)
 dev.off()

##triptych supp
par(mfrow=c(1, 3),oma = c(2.5, 2.5, 1, 1) , cex=0.7 )
#a) juve
plot(density(d0$dist2d[d0$Adult==0], adjust=1.8), xlim=c(0,90) , col="white" , ylim=c(-0.003,0.08) , main="a) Juveniles" , xlab="Distance (m) from habitat edge")

curve(dgamma(x , shape=median(post$a) , scale=median(post$scale))  , add=TRUE , col="black" , lw=5)


for ( i in 1:100 ) {
    curve(dgamma(x , shape=post$a[i] , scale=post$scale[i])  , add=TRUE ,  col=col.alpha("black",0.1))
}

points( d0$dist2d[d0$Adult==0] , (d0$Adult[d0$Adult==0]-0.0011) , pch="J" , col="black")

#b) post medians
plot(density(d0$dist2d[d0$Adult==0], adjust=1.8), xlim=c(0,90) , col="white" , ylim=c(-0.003,0.08) , main="b) Adults & Juveniles" , xlab="Distance (m) from habitat edge")

curve(dgamma(x , shape=median(post$a) , scale=median(post$scale))  , add=TRUE , col="gray" , lw=5)
curve(dgamma(x , shape=median((post$a + post$b_adult)) , scale=median(post$scale) ) , add=TRUE , col="black" , lty=3 , lw=5)


points( d0$dist2d[d0$Adult==0] , (d0$Adult[d0$Adult==0]-0.0011) , pch="J" , col="black")
points( d0$dist2d[d0$Adult==1] , (d0$Adult[d0$Adult==1] -1.0034) , pch="A" , col="black")
legend("topright", c("Adult", "Juvenile") , lty=c(3,1) , col=c("black","gray"), horiz=FALSE, bty='n' , lw=3)

#c_ adults
plot(density(d0$dist2d[d0$Adult==1], adjust=1.8), xlim=c(0,90) , col="white" , ylim=c(-0.003,0.08) , main="c) Adults" , xlab="Distance (m) from habitat edge")
curve(dgamma(x , shape=median((post$a + post$b_adult)) , scale=median(post$scale) ) , add=TRUE , col="black" , lty=1 , lw=5)
points( d0$dist2d[d0$Adult==1] , (d0$Adult[d0$Adult==1] -1.0011) , pch="A" , col="black")

for ( i in 1:100 ) {
    curve(dgamma(x , shape=(post$a[i] + post$b_adult[i]) , scale=post$scale[i] ) , add=TRUE ,  col=col.alpha("black",0.1) , lty=1)
}

mtext("Distance (m) from edge", outer =TRUE, cex = 1, side=1)
mtext("Density", outer =TRUE, cex = 1, side=2 , line=1)



plot(dist2d ~ Adult , data=d0)


m1 <- map2stan(
alist(
Distance ~ dgamma2(mu,scale),
log(mu) ~ a + b_male*male,
a ~ dnorm(0,2),
b_male ~ dnorm(0,2),
scale ~ dexp(2)
),
data=list(
#Male=as.integer(d$Male),
Distance=d$logdist2d,
male=as.integer(d$Male)
),
cores=2 , chains=2 , warmup=5000, iter=10000, constraints=list(scale="lower=0") , WAIC=TRUE
)


m1zag <- map2stan(
alist(
grams_eaten ~ dzagamma2( p, mu , scale ),
logit(p) ~ a + a_species + a_ratid + a_date + (b_Male_species + b_Male)*Male + a_order ,
log(mu) ~  a + a_species + a_ratid + a_date + (b_Male_species + b_Male)*Male + a_order ,
c(a_species , b_Male_species)[species_index] ~ dmvnorm2(0, sigma_species , Rho_species),
a_ratid[rat_index] ~ dnorm(0,sigma_ratid),
a_date[date_index] ~ dnorm(0,sigma_date),
a_order[order_index] ~ dnorm(0,sigma_order),
a ~ dnorm(0,1),
b_Male ~ dnorm(0,1),
sigma_species ~ dcauchy(0,2),
sigma_ratid ~ dcauchy(0,2),
sigma_date ~ dcauchy(0,2),
Rho_species ~ dlkjcorr(2),
sigma_order ~ dcauchy(0,2),
scale ~ dcauchy(0,2)
),
data=list(
grams_eaten=as.numeric(d$grams_eaten),
species_index=as.integer(d$species),
rat_index=as.integer(d$Rat_ID),
Male=as.integer(d$Male),
date_index=d$Date_Index,
order_index=as.integer(d$trial_indiv)
),
cores=2 , chains=2 , warmup=5000, iter=10000, constraints=list(scale="lower=0") , WAIC=TRUE
)

m1zag <- map2stan(
alist(
grams_eaten ~ dzagamma2( p, mu , scale ),
logit(p) ~ a + a_species + a_ratid + a_date + (b_Male_species + b_Male)*Male + a_order ,
log(mu) ~  a + a_species + a_ratid + a_date + (b_Male_species + b_Male)*Male + a_order ,
c(a_species , b_Male_species)[species_index] ~ dmvnorm2(0, sigma_species , Rho_species),
a_ratid[rat_index] ~ dnorm(0,sigma_ratid),
a_date[date_index] ~ dnorm(0,sigma_date),
a_order[order_index] ~ dnorm(0,sigma_order),
a ~ dnorm(0,1),
b_Male ~ dnorm(0,1),
sigma_species ~ dcauchy(0,2),
sigma_ratid ~ dcauchy(0,2),
sigma_date ~ dcauchy(0,2),
Rho_species ~ dlkjcorr(2),
sigma_order ~ dcauchy(0,2),
scale ~ dcauchy(0,2)
),
data=list(
grams_eaten=as.numeric(d$grams_eaten),
species_index=as.integer(d$species),
rat_index=as.integer(d$Rat_ID),
Male=as.integer(d$Male),
date_index=d$Date_Index,
order_index=as.integer(d$trial_indiv)
),
cores=2 , chains=2 , warmup=5000, iter=10000, constraints=list(scale="lower=0") , WAIC=TRUE
)

m2zag <- map2stan(
list(
grams_eaten ~ dzagamma2( p, mu , scale),
logit(p) ~ a + a_species + a_ratid + a_date + a_order +  b_distance*Distance + (b_species_Adult + b_Adult)*Adult ,
log(mu) ~ a + a_species + a_ratid + a_date + a_order +  b_distance*Distance + (b_species_Adult + b_Adult)*Adult ,
c(a_species,b_species_Adult)[species_index] ~ dmvnorm2(0, sigma_species , Rho_species),
a ~ dnorm(0,1),
a_ratid[rat_index] ~ dnorm(0,sigma_ratid),
a_date[date_index] ~ dnorm(0,sigma_date),
a_order[order_index] ~ dnorm(0,sigma_order),
b_distance ~ dnorm(0,1),
b_Adult ~ dnorm(0,1),
sigma_species ~ dcauchy(0,2),
sigma_ratid ~ dcauchy(0,2),
sigma_date ~ dcauchy(0,2),
sigma_order ~ dcauchy(0,2),
Rho_species ~ dlkjcorr(2) ,
scale ~ dcauchy(0,2)
) ,
data=list(
grams_eaten=as.numeric(d$grams_eaten),
species_index=as.integer(d$species),
Distance=d$logdist2d,
Adult=as.integer(d$Adult) ,
rat_index=as.integer(d$Rat_ID),
date_index=d$Date_Index,
order_index=as.integer(d$trial_indiv)
),
cores=2 , chains=2 , warmup=5000, iter=10000, constraints=list(scale="lower=0") , WAIC=TRUE
)



m3zag <- map2stan(
list(
grams_eaten ~ dzagamma2(p, mu , scale ) ,
logit(p) ~ a + a_species + a_ratid + a_date + a_order +  (b_distance_species + b_distance)*Distance ,
log(mu) ~ a + a_species + a_ratid + a_date + a_order +  (b_distance_species + b_distance)*Distance ,
c(a_species,b_distance_species)[species_index] ~ dmvnorm2(0, sigma_species , Rho_species),
a ~ dnorm(0,1),
a_date[date_index] ~ dnorm(0,sigma_date),
a_ratid[rat_index] ~ dnorm(0,sigma_ratid),
a_order[order_index] ~ dnorm(0,sigma_order),
b_distance ~ dnorm(0,1),
scale ~ dcauchy(0,2),
sigma_species ~ dcauchy(0,2),
sigma_ratid ~ dcauchy(0,2),
sigma_date ~ dcauchy(0,2),
sigma_order ~ dcauchy(0,2),
Rho_species ~ dlkjcorr(2)
) ,
data=list(
grams_eaten=as.numeric(d$grams_eaten),
species_index=as.integer(d$species),
Distance=d$dist2d,
rat_index=as.integer(d$Rat_ID),
date_index=d$Date_Index,
order_index=as.integer(d$trial_indiv)
),

cores=2 , chains=2 , warmup=5000, iter=10000, constraints=list(scale="lower=0") , WAIC=TRUE
)



m4zag <- map2stan(
list(
grams_eaten ~ dzagamma2( p, mu , scale ) ,
logit(p) ~ a + a_species + a_ratid + a_date + a_order + (b_species_Adult + b_Adult)*Adult  ,
log(mu) ~ a + a_species + a_ratid + a_date + a_order + (b_species_Adult + b_Adult)*Adult ,
c(a_species,b_species_Adult)[species_index] ~ dmvnorm2(0, sigma_species , Rho_species),
a ~ dnorm(0,1),
a_ratid[rat_index] ~ dnorm(0,sigma_ratid),
a_date[date_index] ~ dnorm(0,sigma_date),
a_order[order_index] ~ dnorm(0,sigma_order),
b_Adult ~ dnorm(0,1),
scale ~ dcauchy(0,2),
sigma_species ~ dcauchy(0,2),
sigma_ratid ~ dcauchy(0,2),
sigma_order ~ dcauchy(0,2),
sigma_date ~ dcauchy(0,2),
Rho_species ~ dlkjcorr(2)
) ,
data=list(
grams_eaten=as.numeric(d$grams_eaten),
species_index=as.integer(d$species),
Adult=as.integer(d$Adult) ,
rat_index=as.integer(d$Rat_ID),
date_index=d$Date_Index,
order_index=as.integer(d$trial_indiv)
),
cores=2 , chains=2 , warmup=5000, iter=10000, constraints=list(scale="lower=0") , WAIC=TRUE
)




m5zag <- map2stan(
list(
grams_eaten ~ dzagamma2(p, mu , scale ) ,
logit(p) ~ a + a_species + a_ratid + a_date + a_order + (b_distance_species + b_distance)*Distance + (b_species_Adult + b_Adult)*Adult  ,
log(mu) ~ a + a_species + a_ratid + a_date + a_order +  (b_distance_species + b_distance)*Distance + (b_species_Adult + b_Adult)*Adult ,
c(a_species,b_species_Adult,b_distance_species)[species_index] ~ dmvnorm2(0, sigma_species , Rho_species),
a ~ dnorm(0,1),
a_ratid[rat_index] ~ dnorm(0,sigma_ratid),
a_date[date_index] ~ dnorm(0,sigma_date),
a_order[order_index] ~ dnorm(0,sigma_order),
b_distance ~ dnorm(0,1),
b_Adult ~ dnorm(0,1),
scale ~ dcauchy(0,2),
sigma_species ~ dcauchy(0,2),
sigma_ratid ~ dcauchy(0,2),
sigma_date ~ dcauchy(0,2),
sigma_order ~ dcauchy(0,2),
Rho_species ~ dlkjcorr(2)
) ,
data=list(
grams_eaten=as.numeric(d$grams_eaten),
species_index=as.integer(d$species),
Distance=d$logdist2d,
Adult=as.integer(d$Adult) ,
rat_index=as.integer(d$Rat_ID),
date_index=as.integer(d$Date_Index),
order_index=as.integer(d$trial_indiv)
),

cores=2 , chains=2 , warmup=5000, iter=10000, constraints=list(scale="lower=0") , WAIC=TRUE
)



m6zag <-  map2stan(
list(
    grams_eaten ~ dzagamma2( p , mu , scale ) ,
    logit(p) ~ a + a_species + a_ratid + a_date + a_order +  (b_distance_species + b_distance)*Distance + (b_species_Adult + b_Adult)*Adult + b_adultXdist*Distance*Adult ,
    log(mu) ~ a + a_species + a_ratid + a_date + a_order +  (b_distance_species + b_distance)*Distance + (b_species_Adult + b_Adult)*Adult + b_adultXdist*Distance*Adult ,
    c(a_species,b_species_Adult,b_distance_species)[species_index] ~ dmvnorm2(0, sigma_species , Rho_species),
    a ~ dnorm(0,1),
    a_ratid[rat_index] ~ dnorm(0,sigma_ratid),
    a_date[date_index] ~ dnorm(0,sigma_date),
    a_order[order_index] ~ dnorm(0,sigma_order),
    b_distance ~ dnorm(0,10),
    b_Adult ~ dnorm(0,10),
    b_adultXdist ~ dnorm(0,10),
    scale ~ dcauchy(0,2),
    sigma_species ~ dcauchy(0,2),
    sigma_ratid ~ dcauchy(0,2),
    sigma_order ~ dcauchy(0,2),
    sigma_date ~ dcauchy(0,2),
    Rho_species ~ dlkjcorr(2)
) ,
data=list(
    grams_eaten=as.numeric(d$grams_eaten),
    species_index=as.integer(d$species),
    Distance=d$logdist2d,
    Adult=as.integer(d$Adult) ,
    rat_index=as.integer(d$Rat_ID),
    date_index=d$Date_Index,
    order_index=as.integer(d$trial_indiv)
),

cores=2 , chains=2 , warmup=5000, iter=10000, constraints=list(scale="lower=0") , WAIC=TRUE , adapt.delta=0.9
)

m6zagfixedprior <-  map2stan(
list(
    grams_eaten ~ dzagamma2( p , mu , scale ) ,
    logit(p) ~ a + a_species + a_ratid + a_date + a_order +  (b_distance_species + b_distance)*Distance + (b_species_Adult + b_Adult)*Adult + b_adultXdist*Distance*Adult ,
    log(mu) ~ a + a_species + a_ratid + a_date + a_order +  (b_distance_species + b_distance)*Distance + (b_species_Adult + b_Adult)*Adult + b_adultXdist*Distance*Adult ,
    a_species[species_index]~ dnorm(0,1),
    b_species_Adult[species_index]~ dnorm(0,1),
    b_distance_species[species_index] ~  dnorm(0,1),
    a ~ dnorm(0,1),
    a_ratid[rat_index] ~ dnorm(0,sigma_ratid),
    a_date[date_index] ~ dnorm(0,sigma_date),
    a_order[order_index] ~ dnorm(0,sigma_order),
    b_distance ~ dnorm(0,10),
    b_Adult ~ dnorm(0,10),
    b_adultXdist ~ dnorm(0,10),
    scale ~ dcauchy(0,2),
   #sigma_species ~ dcauchy(0,2),
    sigma_ratid ~ dcauchy(0,2),
    sigma_order ~ dcauchy(0,2),
    sigma_date ~ dcauchy(0,2)
    #Rho_species ~ dlkjcorr(2)
) ,
data=list(
    grams_eaten=as.numeric(d$grams_eaten),
    species_index=as.integer(d$species),
    Distance=d$logdist2d,
    Adult=as.integer(d$Adult) ,
    rat_index=as.integer(d$Rat_ID),
    date_index=d$Date_Index,
    order_index=as.integer(d$trial_indiv)
),

cores=2 , chains=2 , warmup=5000, iter=10000, constraints=list(scale="lower=0") , WAIC=TRUE
)



m7zag <-  map2stan(
list(
grams_eaten ~ dzagamma2( p , mu , scale ) ,
logit(p) ~ a + a_species + a_ratid + a_date + a_order + (b_distance_species + b_distance)*Distance + (b_species_Adult + b_Adult)*Adult + b_adultXdist*Distance*Adult + (b_Male_species + b_Male)*Male ,
log(mu) ~ a + a_species + a_ratid + a_date + a_order + (b_distance_species + b_distance)*Distance + (b_species_Adult + b_Adult)*Adult + b_adultXdist*Distance*Adult + (b_Male_species + b_Male)*Male ,
c(a_species,b_species_Adult,b_distance_species,b_Male_species)[species_index] ~ dmvnorm2(0, sigma_species , Rho_species),
a ~ dnorm(0,1),
a_ratid[rat_index] ~ dnorm(0,sigma_ratid),
a_date[date_index] ~ dnorm(0,sigma_date),
a_order[order_index] ~ dnorm(0,sigma_order),
b_distance ~ dnorm(0,10),
b_Adult ~ dnorm(0,10),
b_Male ~ dnorm(0,10),
b_adultXdist ~ dnorm(0,10),
scale ~ dcauchy(0,2),
sigma_species ~ dcauchy(0,2),
sigma_ratid ~ dcauchy(0,2),
sigma_date ~ dcauchy(0,2),
sigma_order ~ dcauchy(0,2),
Rho_species ~ dlkjcorr(2)
) ,
data=list(
grams_eaten=as.numeric(d$grams_eaten),
species_index=as.integer(d$species),
Distance=d$logdist2d,
Adult=as.integer(d$Adult) ,
Male=as.integer(d$Male),
rat_index=as.integer(d$Rat_ID),
date_index=d$Date_Index,
order_index=as.integer(d$trial_indiv)
),

cores=2 , chains=2 , warmup=5000 , iter=10000 , constraints=list(scale="lower=0") , WAIC=TRUE
)


m8zag <-  map2stan(
list(
    grams_eaten ~ dzagamma2( p , mu , scale ) ,
    logit(p) ~ a + a_species + a_ratid + a_date + a_order +  (b_distance_species + b_distance)*Distance + (b_species_Adult + b_Adult)*Adult + (b_adultXdist + b_adultXdist_species)*Distance*Adult ,
    log(mu) ~ a + a_species + a_ratid + a_date + a_order +  (b_distance_species + b_distance)*Distance + (b_species_Adult + b_Adult)*Adult + (b_adultXdist + b_adultXdist_species)*Distance*Adult ,
    c(a_species,b_species_Adult,b_distance_species,b_adultXdist_species)[species_index] ~ dmvnorm2(0, sigma_species , Rho_species),
    a ~ dnorm(0,1),
    a_ratid[rat_index] ~ dnorm(0,sigma_ratid),
    a_date[date_index] ~ dnorm(0,sigma_date),
    a_order[order_index] ~ dnorm(0,sigma_order),
    b_distance ~ dnorm(0,10),
    b_Adult ~ dnorm(0,10),
    b_adultXdist ~ dnorm(0,10),
    scale ~ dcauchy(0,2),
    sigma_species ~ dcauchy(0,2),
    sigma_ratid ~ dcauchy(0,2),
    sigma_order ~ dcauchy(0,2),
    sigma_date ~ dcauchy(0,2),
    Rho_species ~ dlkjcorr(2)
) ,
data=list(
    grams_eaten=as.numeric(d$grams_eaten),
    species_index=as.integer(d$species),
    Distance=d$logdist2d,
    Adult=as.integer(d$Adult) ,
    rat_index=as.integer(d$Rat_ID),
    date_index=d$Date_Index,
    order_index=as.integer(d$trial_indiv)
),

cores=2 , chains=2 , warmup=5000, iter=10000, constraints=list(scale="lower=0") , WAIC=TRUE
)

###calc for paper
post <- extract.samples(m6zag)
mean(exp(post$a + post$b_Adult))
###############GRAPHS BELOW


WAICvalues <- compare(m1zag,m2zag,m3zag,m4zag,m5zag,m6zag)
write.csv(WAICvalues@output, "~/Dropbox/Quail Ridge Woodrat Data/woodratdietWAICweights.csv")
woodratdietcoef <- coeftab(m1zag,m2zag,m3zag,m4zag,m5zag,m6zag)
write.csv(woodratdietcoef@coefs, "~/Dropbox/Quail Ridge Woodrat Data/woodratdietcoefs.csv")
m6output <- precis(m6zag, depth=2 , digits=2 , prob=.90)
write.csv(m6output@output, "~/Dropbox/Quail Ridge Woodrat Data/m6output.csv")
m5output <- precis(m5zag, depth=2 , digits=2 , prob=.90)
write.csv(m5output@output, "~/Dropbox/Quail Ridge Woodrat Data/m5output.csv")
m4output <- precis(m4zag, depth=2 , digits=2 , prob=.90)
write.csv(m4output@output, "~/Dropbox/Quail Ridge Woodrat Data/m4output.csv")
m3output <- precis(m3zag, depth=2 , digits=2 , prob=.90)
write.csv(m3output@output, "~/Dropbox/Quail Ridge Woodrat Data/m3output.csv")
m2output <- precis(m2zag, depth=2 , digits=2 , prob=.90)
write.csv(m2output@output, "~/Dropbox/Quail Ridge Woodrat Data/m2output.csv")
m1output <- precis(m1zag, depth=2 , digits=2 , prob=.90)
write.csv(m1output@output, "~/Dropbox/Quail Ridge Woodrat Data/m1output.csv")

png("DietPrefMu.png",res=300,height=7,width=7, units = 'in')

###FOR PAPER
a_species_zeros <- matrix(0,1000,7)
a_rat_zeros <- matrix(0,1000,22)
a_date_zeros <- matrix(0,1000,6)
a_order_zeros <- matrix(0,1000,3)
b_Male_species_zeros <- matrix(0,1000,7)
b_species_Adult_zeros<- matrix(0,1000,7)
b_distance_species_zeros <- matrix(0,1000,7)

plantlist <- c("a. poison oak","b. toyon","c. bay laurel","d. buckeye","e. incense cedar","f. chamise","g. scrub oak" , "h. all plants")
#plantlist <- c("a. Toxicodendron diversilobum","b. Heteromeles arbutifolia","c. Umbellularia californica","d. Aesculus californica","e. Calocedrus decurrens","f. Aesculus californica","g. Quercus berberidifolia" , "h. All Plants")

dist_seq <- seq(from=min(d$logdist2d),to=max(d$logdist2d),length=30)
old.par <- par(mfrow=c(3, 3),oma = c(3, 3, 0, 0) , cex=0.7 , mar=c(1.1,1.1,1,1) )

for (i in 1:6){
    d.pred <- list(
    date_index=rep(1,30),
    rat_index=rep(1,30),
    order_index=rep(1,30),
    species_index=rep(i,30),
    Adult=rep(0,30),
    Distance = dist_seq
    )
    
    link2 <- link(m6zag , data=d.pred,replace=list(a_ratid=a_rat_zeros, a_species=a_species_zeros,a_date=a_date_zeros, b_species_Adult=b_species_Adult_zeros, b_distance_species=b_distance_species_zeros ), WAIC=TRUE)

    pred.p.med <- apply(link2$mu , 2 , median)
    pred.p.PI <- apply(link2$mu , 2 , PI , prob=0.90)
    plot(grams_eaten~dist2d,data=subset(d,species_index==i & Adult==0 & grams_eaten > 0), xlab='' , ylab='', ann=FALSE, ylim=c(0,0.8), xlim=c(1,60) , col="red", pch=1, main=plantlist[i])
    lines(pred.p.med ~ exp(dist_seq) , lw=2, col="red" , lty=1)
    shade(pred.p.PI , exp(dist_seq), col=col.alpha("red", alpha=0.2))
    title(plantlist[i], line = -1 )
    
    
    ##adults
    dist_seq <- seq(from=min(d$logdist2d),to=max(d$logdist2d),length=30)
    
    d.pred <- list(
    date_index=rep(1,30),
    rat_index=rep(1,30),
    species_index=rep(i,30),
    order_index=rep(1,30),
    Adult=rep(1,30),
    Distance = dist_seq
    )
    
    link2 <- link(m6zag , data=d.pred,replace=list(a_ratid=a_rat_zeros, a_species=a_species_zeros,a_date=a_date_zeros, b_species_Adult=b_species_Adult_zeros, b_distance_species=b_distance_species_zeros ), WAIC=TRUE)
    points(grams_eaten~dist2d, data=subset(d,species_index==i & Adult==1 & grams_eaten > 0), pch=17, col="blue")
    pred.p.med <- apply(link2$mu , 2 , mean)
    pred.p.PI <- apply(link2$mu , 2 , PI , prob=0.90)
    lines(pred.p.med ~ exp(dist_seq) , lw=2, col="blue" , lty=2)
    shade(pred.p.PI , exp(dist_seq), col=col.alpha("blue", alpha=0.2))
    #legend("topleft", c("Adult", "Juvenile") , lty=c(1,2) , pch=c(17,1), col=c("black","black"), horiz=FALSE, bty='n')
}

#frame()
i <- 7
d.pred <- list(
date_index=rep(1,30),
rat_index=rep(1,30),
order_index=rep(1,30),
species_index=rep(i,30),
Adult=rep(0,30),
Distance = dist_seq
)

link2 <- link(m6zag , data=d.pred,replace=list(a_ratid=a_rat_zeros, a_species=a_species_zeros,a_date=a_date_zeros, b_species_Adult=b_species_Adult_zeros, b_distance_species=b_distance_species_zeros ), WAIC=TRUE)

pred.p.med <- apply(link2$mu , 2 , mean)
pred.p.PI <- apply(link2$mu , 2 , PI , prob=0.90)
plot(grams_eaten~dist2d,data=subset(d,species_index==i & Adult==0 & grams_eaten > 0), xlab='n' , ylab='n' , ann=FALSE, ylim=c(0,2), xlim=c(1,60) , col="red", pch=1, main=plantlist[i])
lines(pred.p.med ~ exp(dist_seq) , lw=2, col="red" , lty=1)
shade(pred.p.PI , exp(dist_seq), col=col.alpha("red", alpha=0.2))
title(plantlist[i], line = -1)
##adults
dist_seq <- seq(from=min(d$logdist2d),to=max(d$logdist2d),length=30)

d.pred <- list(
date_index=rep(1,30),
rat_index=rep(1,30),
species_index=rep(i,30),
order_index=rep(1,30),
Adult=rep(1,30),
Distance = dist_seq
)

link2 <- link(m6zag , data=d.pred,replace=list(a_ratid=a_rat_zeros, a_species=a_species_zeros,a_date=a_date_zeros, b_species_Adult=b_species_Adult_zeros, b_distance_species=b_distance_species_zeros ), WAIC=TRUE)
points(grams_eaten~dist2d, data=subset(d,species_index==i & Adult==1 & grams_eaten > 0), pch=17, col="blue")
pred.p.med <- apply(link2$mu , 2 , mean)
pred.p.PI <- apply(link2$mu , 2 , PI , prob=0.90)
lines(pred.p.med ~ exp(dist_seq) , lw=2, col="blue" , lty=2)
shade(pred.p.PI , exp(dist_seq), col=col.alpha("blue", alpha=0.2))

#all plants

a_species_zeros <- matrix(0,1000,7)
a_rat_zeros <- matrix(0,1000,22)
a_date_zeros <- matrix(0,1000,6)
a_order_zeros <- matrix(0,1000,3)
b_Male_species_zeros <- matrix(0,1000,7)
b_species_Adult_zeros<- matrix(0,1000,7)
b_distance_species_zeros <- matrix(0,1000,7)

d.pred <- list(
date_index=rep(1,30),
rat_index=rep(1,30),
order_index=rep(1,30),
species_index=rep(1,30),
Adult=rep(0,30),
Distance = dist_seq
)

link2 <- link(m6zag , n=1000 , data=d.pred, replace=list(a_ratid=a_rat_zeros, a_species=a_species_zeros,a_date=a_date_zeros, b_species_Adult=b_species_Adult_zeros, b_distance_species=b_distance_species_zeros , a_order=a_order_zeros ), WAIC=TRUE)

#blue juvenile
pred.p.med <- apply(link2$mu , 2 , mean)
pred.p.PI <- apply(link2$mu  , 2 , PI , prob=0.90)
plot(grams_eaten~dist2d,data=subset(d, Adult==0 & grams_eaten>0), xlab='n' , ylab=
'n' , ann=FALSE, ylim=c(0,2), xlim=c(1,60) , col="red", pch=1 , main="h. All Plants")
lines(pred.p.med ~ exp(dist_seq) , lw=2, col="red" , lty=1)
shade(pred.p.PI , exp(dist_seq), col=col.alpha("red", alpha=0.2))
title(plantlist[8], line = -1)



d.pred <- list(
date_index=rep(1,30),
rat_index=rep(1,30),
species_index=rep(i,30),
order_index=rep(1,30),
Adult=rep(1,30),
Distance = dist_seq
)

link2 <- link(m6zag , n=1000 , data=d.pred, replace=list(a_ratid=a_rat_zeros, a_species=a_species_zeros,a_date=a_date_zeros, b_species_Adult=b_species_Adult_zeros, b_distance_species=b_distance_species_zeros ), WAIC=TRUE)
points(grams_eaten~dist2d, data=subset(d, Adult==1 & grams_eaten>0), pch=17, col="blue")
pred.p.med <- apply(link2$mu , 2 , median)
pred.p.PI <- apply(link2$mu , 2 , PI , prob=0.90)
lines(pred.p.med ~ exp(dist_seq) , lw=2, col="blue" , lty=2 )
shade(pred.p.PI , exp(dist_seq), col=col.alpha("blue", alpha=0.2))

#legend
plot(2,2 , xlab='n' , ylab='n' , ann=FALSE, ylim=c(0,1) , xlim=c(0,1) , col="white", pch=0 ,axes=FALSE)
legend("center", c("Adult", "Juvenile") , lty=c(2,1) , pch=c(17,1), col=c("blue","red"), horiz=FALSE, bty='n', cex=2, lw=2)
#mtext(expression("Figure 2. Model averaged predictions for" ~mu), outer =TRUE, cex = 1.2 , side=3)
mtext("Grams eaten", outer =TRUE, cex = 1.2, side=2, line=1.5)
mtext("Distance (m) from habitat edge", outer =TRUE, cex = 1.2, side=1 , line=1.5)
par(old.par)

dev.off()


#######p

png("DietPrefP.png",res=300,height=7,width=7, units = 'in')

plantlist <- c("a. poison oak","b. toyon","c. bay laurel","d. buckeye","e. incense cedar","f. chamise","g. scrub oak" , "h. all plants")
dist_seq <- seq(from=min(d$logdist2d),to=max(d$logdist2d),length=30)
old.par <- par(mfrow=c(3, 3),oma = c(3, 3, 0, 0) , cex=0.7 , mar=c(1.1,1.1,1,1))

for (i in 1:6){
    d.pred <- list(
    date_index=rep(1,30),
    rat_index=rep(1,30),
    order_index=rep(1,30),
    species_index=rep(i,30),
    Adult=rep(0,30),
    Distance = dist_seq
    )
    
    link2 <- link(m6zag , data=d.pred)
    
    pred.p.med <- apply(link2$p , 2 , mean)
    pred.p.PI <- apply(link2$p , 2 , PI , prob=0.90)
    plot(consumed~dist2d,data=subset(d, Adult==0), xlab='' , ylab='', ann=FALSE, ylim=c(0,.8) , xlim=c(1,60) , col="white", pch=1, main=plantlist[i])
    lines(pred.p.med ~ exp(dist_seq), col="red" , lw=2 , lty=1)
    shade(pred.p.PI , exp(dist_seq) , col=col.alpha("red", alpha=0.2))
    title(plantlist[i], line = -1)
    
    
    ##adults
    dist_seq <- seq(from=min(d$logdist2d),to=max(d$logdist2d),length=30)
    
    d.pred <- list(
    date_index=rep(1,30),
    rat_index=rep(1,30),
    species_index=rep(i,30),
    order_index=rep(1,30),
    Adult=rep(1,30),
    Distance = dist_seq
    )
    
    link2 <- link(m6zag , data=d.pred)
    # points(consumed~logdist2d,data=subset(d, Adult==1), pch=17, col="black")
    pred.p.med <- apply(link2$p , 2 , mean)
    pred.p.PI <- apply(link2$p , 2 , PI , prob=0.90)
    lines(pred.p.med ~ exp(dist_seq) , lw=2, col="blue" , lty=2)
    shade(pred.p.PI , exp(dist_seq), col=col.alpha("blue", alpha=0.2))
    #legend("topleft", c("Adult", "Juvenile") , lty=c(1,2) , pch=c(17,1), col=c("black","black"), horiz=FALSE, bty='n')
}

#frame()
i <- 7
d.pred <- list(
date_index=rep(1,30),
rat_index=rep(1,30),
order_index=rep(1,30),
species_index=rep(i,30),
Adult=rep(0,30),
Distance = dist_seq
)

link2 <- link(m6zag , data=d.pred)

pred.p.med <- apply(link2$p , 2 , mean)
pred.p.PI <- apply(link2$p , 2 , PI , prob=0.90)
plot(consumed~dist2d,data=subset(d,species_index==i & Adult==0 ), xlab='n' , ylab='n' , ann=FALSE, ylim=c(0,0.8), xlim=c(1,60) , col="white", pch=1, main=plantlist[i])
lines(pred.p.med ~ exp(dist_seq) , lw=2, col="red" , lty=1)
shade(pred.p.PI , exp(dist_seq), col=col.alpha("red", alpha=0.2))
title(plantlist[i], line = -1)
##adults

d.pred <- list(
date_index=rep(1,30),
rat_index=rep(1,30),
species_index=rep(i,30),
order_index=rep(1,30),
Adult=rep(1,30),
Distance = dist_seq
)

link2 <- link(m6zag , data=d.pred)
#points(consumed~logdist2d, data=subset(d,species_index==i & Adult==1), pch=17, col="black")
pred.p.med <- apply(link2$p , 2 , mean)
pred.p.PI <- apply(link2$p , 2 , PI , prob=0.90)
lines(pred.p.med ~ exp(dist_seq) , lw=2, col="blue" , lty=2)
shade(pred.p.PI , exp(dist_seq), col=col.alpha("blue", alpha=0.2))
#legend("topleft", c("Adult", "Juvenile") , lty=c(1,2) , pch=c(17,1), col=c("black","black"), horiz=FALSE, bty='n')


##
#all plants

a_species_zeros <- matrix(0,1000,7)
a_rat_zeros <- matrix(0,1000,22)
a_date_zeros <- matrix(0,1000,6)
a_order_zeros <- matrix(0,1000,3)
b_Male_species_zeros <- matrix(0,1000,7)
b_species_Adult_zeros<- matrix(0,1000,7)
b_distance_species_zeros <- matrix(0,1000,7)

d.pred <- list(
date_index=rep(1,30),
rat_index=rep(1,30),
order_index=rep(1,30),
species_index=rep(1,30),
Adult=rep(0,30),
Distance = dist_seq
)

link2 <- link(m6zag , n=1000 , data=d.pred, replace=list(a_ratid=a_rat_zeros, a_species=a_species_zeros,a_date=a_date_zeros, b_species_Adult=b_species_Adult_zeros, b_distance_species=b_distance_species_zeros , a_order=a_order_zeros ), WAIC=TRUE)

#blue juvenile
pred.p.med <- apply(link2$p , 2 , mean)
pred.p.PI <- apply(link2$p  , 2 , PI , prob=0.90)
plot(grams_eaten~dist2d,data=subset(d, Adult==0 & grams_eaten>0), xlab='n' , ylab=
'n' , ann=FALSE, ylim=c(0,0.8), xlim=c(1,60) , col="white", pch=1 , main="h. All Plants")
lines(pred.p.med ~ exp(dist_seq) , lw=2, col="red" , lty=1)
shade(pred.p.PI , exp(dist_seq), col=col.alpha("red", alpha=0.2))
title(plantlist[8], line = -1)



d.pred <- list(
date_index=rep(1,30),
rat_index=rep(1,30),
species_index=rep(1,30),
order_index=rep(1,30),
Adult=rep(1,30),
Distance = dist_seq
)

link2 <- link(m6zag , n=1000 , data=d.pred, replace=list(a_ratid=a_rat_zeros, a_species=a_species_zeros,a_date=a_date_zeros, b_species_Adult=b_species_Adult_zeros, b_distance_species=b_distance_species_zeros ), WAIC=TRUE)
#points(grams_eaten~logdist2d, data=subset(d, Adult==1 & grams_eaten>0), pch=17, col="black")
pred.p.med <- apply(link2$p , 2 , mean)
pred.p.PI <- apply(link2$p , 2 , PI , prob=0.90)
lines(pred.p.med ~ exp(dist_seq) , lw=2, col="blue" , lty=2 )
shade(pred.p.PI , exp(dist_seq), col=col.alpha("blue", alpha=0.2))

#legend
plot(2,2 , xlab='n' , ylab='n' , ann=FALSE, ylim=c(0,1) , xlim=c(0,1) , col="white", pch=0 ,axes=FALSE)
legend("center", c("Adult", "Juvenile") , lty=c(2,1) , lw=2,col=c("blue","red"), horiz=FALSE, bty='n', cex=2)

#mtext(expression("Figure 1. Model averaged predictions for" ~p), outer =TRUE, cex = 1.2 , side=3)
###
mtext("Probabilty of sampling plant", outer =TRUE, cex = 1.2, side=2, line=1.5)
mtext("Distance (m) from habitat edge", outer =TRUE, cex = 1.2, side=1, line=1.5)
par(old.par)
dev.off()


col_index <- c("red" , "orange" , "gold" , "green" ,"blue" , "slateblue" , "violet")
plot(grams_eaten~Rat_ID , data=d, col="white" , pch=19 )
for (i in 1:7){
    points(grams_eaten~Rat_ID , data=subset(d , species_index==i), col=col_index[i] , pch=19 )
}

##########calculate numbers
dist_seq <- seq(from=min(d$logdist2d),to=max(d$logdist2d),length=30)

a_species_zeros <- matrix(0,1000,7)
a_rat_zeros <- matrix(0,1000,22)
a_date_zeros <- matrix(0,1000,6)
a_order_zeros <- matrix(0,1000,3)
b_Male_species_zeros <- matrix(0,1000,7)
b_species_Adult_zeros<- matrix(0,1000,7)
b_distance_species_zeros <- matrix(0,1000,7)

d.pred <- list(
date_index=c(1,1,1),
rat_index=c(1,1,1),
order_index=c(1,2,3),
species_index=c(5,6,7),
Adult=rep(mean(d$Adult),3),
#Adult=c(0,0,1),
Distance = rep(mean(d$logdist2d),3),
Male = 0
)

#link2 <- link(m6zag , n=1000 , data=d.pred, replace=list(a_ratid=a_rat_zeros, a_species=a_species_zeros,a_date=a_date_zeros, b_species_Adult=b_species_Adult_zeros, b_distance_species=b_distance_species_zeros , a_order=a_order_zeros ), WAIC=TRUE)
#link2 <- link(m6zag , n=1000 , data=d.pred, replace=list(a_ratid=a_rat_zeros,a_date=a_date_zeros, a_order=a_order_zeros ), WAIC=TRUE)
link2 <- link(m6zag , n=1000 , data=d.pred, replace=list(a_ratid=a_rat_zeros, a_species=a_species_zeros,a_date=a_date_zeros, b_species_Adult=b_species_Adult_zeros, b_distance_species=b_distance_species_zeros  ), WAIC=TRUE)

###
###

pred.p.med <- apply(link2$p , 2 , mean)
pred.p.PI <- apply(link2$p , 2 , PI , prob=0.90)
pred.mu.med <- apply(link2$mu , 2 , mean)
pred.mu.PI <- apply(link2$mu , 2 , PI , prob=0.90)

pred.med <- apply(link2$p*link2$mu , 2 , mean)
pred.PI <- apply(link2$p*link2$mu , 2 , PI , prob=0.90)

pred.mu.med 
pred.mu.PI 



##########I SHOULD FORGET ABOUT THE BELOW BUT I AM A CODE HOARDER###################
########PLOTS FOR ZA gamma using ensemble
###JOINT POSTERIORZ
dist_seq <- seq(from=min(d$logdist2d),to=max(d$logdist2d),length=30)

a_species_zeros <- matrix(0,1000,7)
a_rat_zeros <- matrix(0,1000,22)
a_date_zeros <- matrix(0,1000,6)
a_order_zeros <- matrix(0,1000,3)
b_Male_species_zeros <- matrix(0,1000,7)
b_species_Adult_zeros<- matrix(0,1000,7)
b_distance_species_zeros <- matrix(0,1000,7)

d.pred <- list(
date_index=rep(1,30),
rat_index=rep(1,30),
order_index=rep(1,30),
species_index=rep(1,30),
Adult=rep(0,30),
Distance = dist_seq,
Male = rep(0,30)
)

link2 <- link(m6zag , n=1000 , data=d.pred, replace=list(a_ratid=a_rat_zeros, a_species=a_species_zeros,a_date=a_date_zeros, b_species_Adult=b_species_Adult_zeros, b_distance_species=b_distance_species_zeros , a_order=a_order_zeros ), WAIC=TRUE)

#blue juvenile
pred.joint <- link2$p*link2$mu
pred.p.med <- apply(pred.joint  , 2 , median)
pred.p.PI <- apply(pred.joint  , 2 , PI , prob=0.90)
plot(grams_eaten~logdist2d,data=subset(d, Adult==0), xlab="2D Distance (log m)" , ylab=
"Grams Eaten" , ylim=c(0,2) , col="blue", pch=19 , main="Joint Predictions for ZAG Mixture Model")
lines(pred.p.med ~ dist_seq , lw=2, col="blue" , lty=2)
shade(pred.p.PI , dist_seq, col=col.alpha("blue", alpha=0.05))


d.pred <- list(
date_index=rep(1,30),
rat_index=rep(1,30),
species_index=rep(1,30),
order_index=rep(1,30),
Adult=1,
Distance = rep(0,30),
Male = rep(0,30)
)

link2 <- link(m6zag , n=1000 , data=d.pred, replace=list(a_ratid=a_rat_zeros, a_species=a_species_zeros,a_date=a_date_zeros, b_species_Adult=b_species_Adult_zeros, b_distance_species=b_distance_species_zeros ), WAIC=TRUE)
points(grams_eaten~logdist2d, data=subset(d, Adult==1), pch=19, col="red")
pred.joint <- link2$p*link2$mu
pred.p.med <- apply(pred.joint  , 2 , median)
pred.p.PI <- apply(pred.joint  , 2 , PI , prob=0.90)
lines(pred.p.med ~ dist_seq , lw=2, col="red" )
shade(pred.p.PI , dist_seq, col=col.alpha("red", alpha=0.05))
legend("topleft", inset=.01, title="Sex", c("Juvenile" , "Adult"), pch=19, col=c("blue" , "red") , lty=c(2,1), horiz=FALSE )

###JUST MU
##START HERE
dist_seq <- seq(from=min(d$logdist2d),to=max(d$logdist2d),length=30)

a_species_zeros <- matrix(0,1000,7)
a_rat_zeros <- matrix(0,1000,22)
a_date_zeros <- matrix(0,1000,6)
a_order_zeros <- matrix(0,1000,3)
b_Male_species_zeros <- matrix(0,1000,7)
b_species_Adult_zeros<- matrix(0,1000,7)
b_distance_species_zeros <- matrix(0,1000,7)

d.pred <- list(
date_index=rep(1,30),
rat_index=rep(1,30),
order_index=rep(1,30),
species_index=rep(1,30),
Adult=rep(0,30),
Distance = dist_seq,
Male = rep(0,30)
)

link2 <- link(m6zag , n=1000 , data=d.pred, replace=list(a_ratid=a_rat_zeros, a_species=a_species_zeros,a_date=a_date_zeros, b_species_Adult=b_species_Adult_zeros, b_distance_species=b_distance_species_zeros , a_order=a_order_zeros ), WAIC=TRUE)

#blue juvenile
pred.p.med <- apply(link2$mu , 2 , median)
pred.p.PI <- apply(link2$mu  , 2 , PI , prob=0.90)
plot(grams_eaten~logdist2d,data=subset(d, Adult==0 & grams_eaten>0), xlab="2D Distance (log m)" , ylab=
"Grams Eaten" , ylim=c(0,2) , col="blue", pch=19 , main="Predictions for ZAG Mixture Model")
lines(pred.p.med ~ dist_seq , lw=2, col="blue" , lty=2)
shade(pred.p.PI , dist_seq, col=col.alpha("blue", alpha=0.05))


d.pred <- list(
date_index=rep(1,30),
rat_index=rep(1,30),
species_index=rep(1,30),
order_index=rep(1,30),
Adult=1,
Distance = rep(0,30),
Male = rep(0,30)
)

link2 <- link(m6zag , n=1000 , data=d.pred, replace=list(a_ratid=a_rat_zeros, a_species=a_species_zeros,a_date=a_date_zeros, b_species_Adult=b_species_Adult_zeros, b_distance_species=b_distance_species_zeros ), WAIC=TRUE)
points(grams_eaten~logdist2d, data=subset(d, Adult==1 & grams_eaten>0), pch=19, col="red")
pred.joint <- link2$p*link2$mu
pred.p.med <- apply(link2$mu , 2 , median)
pred.p.PI <- apply(link2$mu , 2 , PI , prob=0.90)
lines(pred.p.med ~ dist_seq , lw=2, col="red" )
shade(pred.p.PI , dist_seq, col=col.alpha("red", alpha=0.05))
legend("topleft", inset=.01, title="Sex", c("Juvenile" , "Adult"), pch=19, col=c("blue" , "red") , lty=c(2,1), horiz=FALSE )


###JUST P

##START HERE
dist_seq <- seq(from=min(d$logdist2d),to=max(d$logdist2d),length=30)

a_species_zeros <- matrix(0,1000,7)
a_rat_zeros <- matrix(0,1000,22)
a_date_zeros <- matrix(0,1000,6)
a_order_zeros <- matrix(0,1000,3)
b_Male_species_zeros <- matrix(0,1000,7)
b_species_Adult_zeros<- matrix(0,1000,7)
b_distance_species_zeros <- matrix(0,1000,7)

d.pred <- list(
date_index=rep(1,30),
rat_index=rep(1,30),
order_index=rep(1,30),
species_index=rep(1,30),
Adult=rep(0,30),
Distance = dist_seq,
Male = rep(0,30)
)

link2 <- link(m6zag , n=1000 , data=d.pred, replace=list(a_ratid=a_rat_zeros, a_species=a_species_zeros,a_date=a_date_zeros, b_species_Adult=b_species_Adult_zeros, b_distance_species=b_distance_species_zeros , a_order=a_order_zeros ), WAIC=TRUE)

#blue juvenile
pred.p.med <- apply(link2$p , 2 , median)
pred.p.PI <- apply(link2$p  , 2 , PI , prob=0.90)
plot(consumed~logdist2d,data=subset(d, Adult==0), xlab="log distance from edge (m)" , ylab="Prob Eating Plant" , ylim=c(0,0.8) , col=col.alpha("white", alpha=0.05), pch=19 , main="Predictions p (1/0) for ZAG Mixture Model")
lines(pred.p.med ~ dist_seq, col="blue" , lw=2 , lty=2)
shade(pred.p.PI , dist_seq , col=col.alpha("blue", alpha=0.05))

#red adult
d.pred <- list(
date_index=rep(1,30),
rat_index=rep(1,30),
species_index=rep(1,30),
order_index=rep(1,30),
Adult=1,
Distance = rep(0,30),
Male = rep(0,30)
)

link2 <- link(m6zag , n=1000 , data=d.pred, replace=list(a_ratid=a_rat_zeros, a_species=a_species_zeros,a_date=a_date_zeros, b_species_Adult=b_species_Adult_zeros, b_distance_species=b_distance_species_zeros ), WAIC=TRUE)
pred.p.med <- apply(link2$p , 2 , median)
pred.p.PI <- apply(link2$p  , 2 , PI , prob=0.90)
#points(consumed~logdist2d,data=subset(d, Adult==1), col=col.alpha("red", alpha=0.05), pch=19 )
lines(pred.p.med ~ dist_seq, col="red" , lw=2 , lty=1)
shade(pred.p.PI , dist_seq , col=col.alpha("red", alpha=0.05))
legend("topleft", inset=.01, title="Sex", c("Juvenile" , "Adult"), pch=19, col=c("blue" , "red") , lty=c(2,1), horiz=FALSE )





#######ENSEMBLE###########
dist_seq <- seq(from=min(d$logdist2d),to=max(d$logdist2d),length=30)

d.pred <- list(
date_index=rep(1,30),
rat_index=rep(1,30),
species_index=rep(1,30),
Adult=rep(0,30),
Distance = dist_seq,
Male=rep(0,30)
)

rzagamma2 <- function (n, p,mu,scale) {
    z <- rbinom(n, size = 1, prob = p)
    y <- (1 - z) * rgamma(n, shape=mu/scale, scale = scale)
    return(y)
}

compare(m6zag,m2zag)
coeftab(m6zag,m2zag)
ens <- ensemble(m6zag,m5zag,m4zag,m3zag,m1zag,m2zag, n=1000 , data=d.pred ,replace=list(a_ratid=a_rat_zeros, a_species=a_species_zeros, a_date=a_date_zeros, b_species_Adult=b_species_Adult_zeros, b_distance_species=b_distance_species_zeros ), WAIC=TRUE)



################Plots for ZA Gamma
###posterior densities for
post <- extract.samples(m6zag)
plot(density(exp(post$a + post$b_Adult)), main="Posterior Distributions for Plant Mass Eaten- Age" , xlab="Mass plants eaten (g)")
lines(density(exp(post$a)), lty=3 , lw=2)
legend("topright", inset=.01, title="Age", c("Juvenile" , "Adult"), lty=c(3,1), horiz=FALSE)
#abline(v=exp(median(post$a)))
#abline(v=exp(median(post$a + post$b_Adult)))


post <- extract.samples(m6zag)
plot(density(exp(post$a + post$b_Male)), main="Posterior Distributions for Plant Mass Eaten- Sex" , xlab="Mass plants eaten (g)")
lines(density(exp(post$a)), lty=3 , lw=3)
legend("topright", inset=.01, title="Age", c("Female" , "Male"), lty=c(3,1), horiz=FALSE)
abline(v=exp(median(post$a)))
abline(v=exp(median(post$a + post$b_Male)))

#plot for all foods over log distance
post <- extract.samples(m6zag)
col_index <- c("red", "orange" , "gold" , "green" , "blue" , "violet" , "black")
mu <- function(distance, species_index=0)
    exp( post$a + post$a_species[,species_index] + (post$b_distance + post$b_distance_species[,species_index])*distance )


plot(grams_eaten~logdist2d, data=subset(d, species_index==1), pch=19, col="white", ylim=c(0,1.5) ,  xlim=c(1,max(d$logdist2d)) , main="Grams Eaten of Each Food Type vs. Distance of Nest Site to Habitat Boundary" , cex.main=1 , ylab="Grams Eaten" , xlab="Log Distance (m) from Habitat Boundary")


for (i in 1:7){
mu.i <- sapply(dist_seq , function(log_dist) median(mu(log_dist,i)))
mu.ci.i <- sapply(dist_seq , function(log_dist) HPDI(mu(log_dist,i)))
lines(dist_seq , mu.i , lty=1 , col=col_index[i] , lw=2  )
#lines(log_dist_seq , mu.ci.i[1,] , lty=2 , col=col_index[i])
#lines(log_dist_seq , mu.ci.i[2,] , lty=2 , col=col_index[i])
points(grams_eaten~logdist2d, data=subset(d, species_index==i), pch=1, col=col_index[i])
}
legend("topright", inset=.01, title="Plant Species", c("Poison Oak","Toyon" , "Bay" , "Buckeye" , "Cedar" , "chamise" , "Oak"), pch=19, col=col_index, horiz=FALSE)


##plot age differences for each plant using var slops of the effect of being adult
post <- extract.samples(m7zag)

mu <- function(Adult, species_index=0)
    exp( post$a + post$a_species[,species_index] + (post$b_Adult + post$b_species_Adult[,species_index])*Adult )

Adult_seq <- c(0,1)

##create staggered species index for plot visualization
d$species_index_stag <- ifelse(d$Adult==0,d$species_index + 0.1, d$species_index - 0.1)

plot(grams_eaten~species_index_stag, data=subset(d, Adult==0) , pch=1, col=col.alpha( "red" , alpha = 0.2 ), ylim=c(0,2.1) , main="Age Differences for Each Plant Type Eaten" , cex.main=1 , ylab="Grams Eaten" , xlab="Species" , xaxt="n" , xlim =c(.75,7.25))

points(grams_eaten ~ species_index_stag, data=subset(d, Adult==1) , pch=1, col=col.alpha( "blue" , alpha = 0.2 ))

for (i in 1:7){
mu.i <- sapply(Adult_seq ,function(Adult) mean(mu(Adult,i)))
mu.ci.i <- sapply(Adult_seq , function(Adult) HPDI(mu(Adult,i)))
points( i+0.1 , mu.i[1] , col="red" , pch=17 , cex=1.1  )
points( i-0.1 , mu.i[2] , col="blue" , pch=19, cex=1.1  )
segments(i+0.1, mu.ci.i[1,1],i+0.1, mu.ci.i[2,1],  , lty=2 , col="red")
segments(i-0.1, mu.ci.i[1,2],i-0.1, mu.ci.i[2,2],  , lty=2 , col="blue")
}
axis( 1 , at=1:7 , labels=c("Poison Oak","Toyon" , "Bay" , "Buckeye" , "Cedar" , "chamise" , "Oak*") , font=6 , , las=1 ,cex.axis=.7 )

legend("topleft", inset=.01, title="Age", c("Juvenile","Adult"), pch=c(17,19), col=c("red","blue"), horiz=FALSE)

plantlist <- c("Poison Oak","Toyon","Bay","Buckeye","Cedar","chamise","Oak")

old.par <- par(mfrow=c(3, 3))

for (i in 1:6){
plot(density(exp(post$a + post$a_species[,i] + (post$b_Adult + post$b_species_Adult[,i])*1)), main=plantlist[i] , xlab="Mass eaten (g)", col=col_index[i])
lines(density(exp(post$a + post$a_species[,i])), lty=3 , lw=2, col=col_index[i])
legend("topright", inset=.01, title="Age", c("Juvenile" , "Adult"), lty=c(3,1),col=col_index[i], horiz=FALSE)
}
plot(1~1, pch=0 , col="white", xaxt="n",yaxt="n")

plot(density(exp(post$a + post$a_species[,7] + (post$b_Adult + post$b_species_Adult[,7])*1)), main=plantlist[7], xlab="Mass eaten (g)", xlim=c(0,3), col=col_index[7])
lines(density(exp(post$a + post$a_species[,7])), lty=3 , lw=2, col=col_index[7])
legend("topright", inset=.01, title="Age", c("Juvenile" , "Adult"), lty=c(3,1),col=col_index[7], horiz=FALSE)

plot(1~1, pch=0 , col="white", xaxt="n",yaxt="n")
#legend("left", inset=.01, title="Plant Species", c("Poison Oak","Toyon" , "Bay" , "Buckeye" , "Cedar" , "chamise" , "Oak"), pch=19, col=col_index, horiz=FALSE)
#legend("right", inset=.01, title="Age", c("Juvenile","Adult"), lty=c(3,1), col=1, horiz=FALSE)
par(old.par)

i <- 7
plot(density(d$grams_eaten[d$Adult==1 & d$species_index==i]), main="Posteriors for Oak Mass Eaten-Age", xlim=c(0,2.5) , xlab="Mass eaten (g)")
lines(density(d$grams_eaten[d$Adult==0 & d$species_index==i]), main="Posteriors for Oak Mass Eaten-Age" , lty=3, xlim=c(0,2.5))
lines(density(exp(post$a + post$a_species[,i] + (post$b_Adult + post$b_species_Adult[,i])*1)), lty=1 , col="green")
lines(density(exp(post$a + post$a_species[,i])), lty=3 , lw=2, col="green")
legend("topright", inset=.01, title="Age", c("Juvenile" , "Adult"), lty=c(3,1), horiz=FALSE)


#######plot sex differences for each plant using var slops of the effect of sex
post <- extract.samples(m6zag)
mu <- function(Male, species_index=0)
    exp( post$a + post$a_species[,species_index] + (post$b_Male + post$b_Male_species[,species_index])*Male )
Male_seq <- c(0,1)
##create staggered species index for plot visualization
d$species_index_stag <- ifelse(d$Male==0,d$species_index + 0.1, d$species_index - 0.1)

plot(grams_eaten~species_index_stag, data=subset(d, Male==0) , pch=1, col=col.alpha( "magenta" , alpha = 0.2 ), ylim=c(0,2.1) , main="Sex Differences for Each Plant Type Eaten" , cex.main=1 , ylab="Grams Eaten" , xlab="Species" , xaxt="n" , xlim =c(.75,7.25))
points(grams_eaten ~ species_index_stag, data=subset(d, Male==1) , pch=1, col=col.alpha( "slateblue" , alpha = 0.2 ))

for (i in 1:7){
mu.i <- sapply(Male_seq ,function(Male) median(mu(Male,i)))
mu.ci.i <- sapply(Male_seq , function(Male) HPDI(mu(Male,i)))
points( i+0.1 , mu.i[1] , col="magenta" , pch=17 , cex=1.1  )
points( i-0.1 , mu.i[2] , col="slateblue" , pch=19, cex=1.1  )
segments(i+0.1, mu.ci.i[1,1],i+0.1, mu.ci.i[2,1],  , lty=2 , col="magenta")
segments(i-0.1, mu.ci.i[1,2],i-0.1, mu.ci.i[2,2],  , lty=2 , col="slateblue")
}
axis( 1 , at=1:7 , labels=c("Poison Oak","Toyon" , "Bay" , "Buckeye" , "Cedar" , "chamise" , "Oak*") , font=6 , , las=1 ,cex.axis=.7 )

legend("topleft", inset=.01, title="Age", c("Female","Male"), pch=c(17,19), col=c("magenta","slateblue"), horiz=FALSE)



##plot interaction

col_index <- c("red", "orange" , "gold" , "green" , "blue" , "violet" , "black")
mu <- function(distance, species_index=0)
    exp( post$a + post$a_species[,species_index] + (post$b_distance + post$b_distance_species[,species_index])*distance + (post$b_species_Adult[,species_index] + post$b_Adult)*Adult + post$b_adultXdist*distance*Adult )

log_dist_seq <- seq(from=0, to=5, by=0.1)

old.par <- par(mfrow=c(3, 3))

for (i in 1:7){
plot(grams_eaten~logdist2d, data=subset(d, species_index==i & Adult==0), pch=1, col=col_index[i], ylim=c(0,1) ,  xlim=c(0,max(d$logdist2d))  , cex.main=1 , ylab='n' , xlab='n' )

Adult <- 0
mu.i <- sapply(log_dist_seq , function(log_dist) median(mu(log_dist,i)))
mu.ci.i <- sapply(log_dist_seq , function(log_dist) HPDI(mu(log_dist,i)))
lines(log_dist_seq , mu.i , lty=5 , col=col_index[i] , lw=2  )
lines(log_dist_seq , mu.ci.i[1,] , lty=3 , col=col_index[i])
lines(log_dist_seq , mu.ci.i[2,] , lty=3 , col=col_index[i])
points(grams_eaten~logdist2d, data=subset(d, species_index==i & Adult==1), pch=19, col=col_index[i])

Adult <- 1
mu.i <- sapply(log_dist_seq , function(log_dist) median(mu(log_dist,i)))
mu.ci.i <- sapply(log_dist_seq , function(log_dist) HPDI(mu(log_dist,i)))
lines(log_dist_seq , mu.i , lty=1 , col=col_index[i] , lw=2  )
lines(log_dist_seq , mu.ci.i[1,] , lty=1 , col=col_index[i])
lines(log_dist_seq , mu.ci.i[2,] , lty=1 , col=col_index[i])
points(grams_eaten~logdist2d, data=subset(d, species_index==i), pch=1, col=col_index[i])
}
plot(1~1, pch=0 , col="white", xaxt=0,yaxt=0)
legend("left", inset=.01, title="Plant Species", c("Poison Oak","Toyon" , "Bay" , "Buckeye" , "Cedar" , "Chamise" , "Oak"), pch=19, col=col_index, horiz=FALSE)
legend("right", inset=.01, title="Age", c("Adult", "Juvenile") , lty=c(1,3) , pch=c(19,1), col=1, horiz=FALSE)


par(old.par)


####plot interaction of age species and distance using link

##Just MU
plantlist <- c("Poison Oak","Toyon","Bay","Buckeye","Cedar","Chamise","Oak")
dist_seq <- seq(from=min(d$logdist2d),to=max(d$logdist2d),length=30)
old.par <- par(mfrow=c(3, 3),oma = c(2, 2, 2, 0) , cex=0.5)

for (i in 1:6){
    d.pred <- list(
    date_index=rep(1,30),
    rat_index=rep(1,30),
    order_index=rep(1,30),
    species_index=rep(i,30),
    Adult=rep(0,30),
    Distance = dist_seq
    )


    link2 <- link(m6zag , data=d.pred)

    pred.p.med <- apply(link2$mu , 2 , median)
    pred.p.PI <- apply(link2$mu , 2 , PI , prob=0.90)
    plot(grams_eaten~dist2d,data=subset(d,species_index==i & Adult==0 & grams_eaten > 0), xlab='' , ylab='', ann=FALSE, ylim=c(0,0.8) , col="black", pch=1, main=plantlist[i])
    lines(pred.p.med ~ exp(dist_seq) , lw=2, col="black" , lty=2)
    shade(pred.p.PI , exp(dist_seq), col=col.alpha("black", alpha=0.05))
    title(plantlist[i], line = -1)


    ##adults
    dist_seq <- seq(from=min(d$logdist2d),to=max(d$logdist2d),length=30)

    d.pred <- list(
    date_index=rep(1,30),
    rat_index=rep(1,30),
    species_index=rep(i,30),
    order_index=rep(1,30),
    Adult=rep(1,30),
    Distance = dist_seq
    )

    link2 <- link(m6zag , data=d.pred)
    points(grams_eaten~dist2d, data=subset(d,species_index==i & Adult==1 & grams_eaten > 0), pch=19, col="red")
    pred.p.med <- apply(link2$mu , 2 , median)
    pred.p.PI <- apply(link2$mu , 2 , PI , prob=0.95)
    lines(pred.p.med ~ exp(dist_seq) , lw=2, col="black" , lty=1)
    shade(pred.p.PI , exp(dist_seq), col=col.alpha("black", alpha=0.05))
    legend("topleft", c("Adult", "Juvenile") , lty=c(1,2) , pch=c(19,1), col=c("black"), horiz=FALSE, bty='n')
}

frame()
i <- 7
    d.pred <- list(
    date_index=rep(1,30),
    rat_index=rep(1,30),
    order_index=rep(1,30),
    species_index=rep(i,30),
    Adult=rep(0,30),
    Distance = dist_seq
    )

    link2 <- link(m6zag , data=d.pred)

    pred.p.med <- apply(link2$mu , 2 , median)
    pred.p.PI <- apply(link2$mu , 2 , PI , prob=0.90)
    plot(grams_eaten~logdist2d,data=subset(d,species_index==i & Adult==0 & grams_eaten > 0), xlab='n' , ylab='n' , ann=FALSE, ylim=c(0,2) , col="blue", pch=1, main=plantlist[i])
    lines(pred.p.med ~ dist_seq , lw=2, col="black" , lty=2)
    shade(pred.p.PI , dist_seq, col=col.alpha("black", alpha=0.05))
    title(plantlist[i], line = -1)
    ##adults
    dist_seq <- seq(from=min(d$logdist2d),to=max(d$logdist2d),length=30)

    d.pred <- list(
    date_index=rep(1,30),
    rat_index=rep(1,30),
    species_index=rep(i,30),
    order_index=rep(1,30),
    Adult=rep(1,30),
    Distance = dist_seq
    )

    link2 <- link(m6zag , data=d.pred)
    points(grams_eaten~logdist2d, data=subset(d,species_index==i & Adult==1 & grams_eaten > 0), pch=19, col="red")
    pred.p.med <- apply(link2$mu , 2 , median)
    pred.p.PI <- apply(link2$mu , 2 , PI , prob=0.95)
    lines(pred.p.med ~ dist_seq , lw=2, col="black" , lty=1)
    shade(pred.p.PI , dist_seq, col=col.alpha("black", alpha=0.05))
    legend("topleft", c("Adult", "Juvenile") , lty=c(1,2) , pch=c(19,1), col=c("black"), horiz=FALSE, bty='n')
mtext("Effect of Distance and Age on Diet Preference- MU", outer =TRUE, cex = 1.5 , side=3)
mtext("Grams Eaten", outer =TRUE, cex = 1.2, side=2)
mtext("log distance (m) from edge", outer =TRUE, cex = 1.2, side=1)
par(old.par)

###JOint
plantlist <- c("Poison Oak","Toyon","Bay","Buckeye","Cedar","Chamise","Oak")
dist_seq <- seq(from=min(d$logdist2d),to=max(d$logdist2d),length=30)
old.par <- par(mfrow=c(3, 3),oma = c(2, 2, 2, 0) , cex=0.5)

for (i in 1:6){
    d.pred <- list(
    date_index=rep(1,30),
    rat_index=rep(1,30),
    order_index=rep(1,30),
    species_index=rep(i,30),
    Adult=rep(0,30),
    Distance = dist_seq
    )
    
    link2 <- link(m6zag , data=d.pred)
    pred.joint <- link2$p*link2$mu
    pred.p.med <- apply(pred.joint , 2 , median)
    pred.p.PI <- apply(pred.joint , 2 , PI , prob=0.90)
    plot(grams_eaten~logdist2d,data=subset(d,species_index==i & Adult==0 ), xlab='' , ylab='', ann=FALSE, ylim=c(0,0.8) , col="blue", pch=1, main=plantlist[i])
    lines(pred.p.med ~ dist_seq , lw=2, col="blue" , lty=2)
    shade(pred.p.PI , dist_seq, col=col.alpha("blue", alpha=0.05))
    title(plantlist[i], line = -1)
    
    
    ##adults
    dist_seq <- seq(from=min(d$logdist2d),to=max(d$logdist2d),length=30)
    
    d.pred <- list(
    date_index=rep(1,30),
    rat_index=rep(1,30),
    species_index=rep(i,30),
    order_index=rep(1,30),
    Adult=rep(1,30),
    Distance = dist_seq
    )
    
    link2 <- link(m6zag , data=d.pred)
    points(grams_eaten~logdist2d, data=subset(d,species_index==i & Adult==1 ), pch=19, col="red")
    pred.joint <- link2$p*link2$mu
    pred.p.med <- apply(pred.joint, 2 , median)
    pred.p.PI <- apply(pred.joint , 2 , PI , prob=0.90)
    lines(pred.p.med ~ dist_seq , lw=2, col="red" , lty=1)
    shade(pred.p.PI , dist_seq, col=col.alpha("red", alpha=0.05))
    legend("topleft", c("Adult", "Juvenile") , lty=c(1,2) , pch=c(19,1), col=c("red","blue"), horiz=FALSE, bty='n')
}

frame()
i <- 7
d.pred <- list(
date_index=rep(1,30),
rat_index=rep(1,30),
order_index=rep(1,30),
species_index=rep(i,30),
Adult=rep(0,30),
Distance = dist_seq
)

link2 <- link(m6zag , data=d.pred)
pred.joint <- link2$p*link2$mu

pred.p.med <- apply(pred.joint , 2 , median)
pred.p.PI <- apply(pred.joint , 2 , PI , prob=0.90)
plot(grams_eaten~logdist2d,data=subset(d,species_index==i & Adult==0 ), xlab='n' , ylab='n' , ann=FALSE, ylim=c(0,2) , col="blue", pch=1, main=plantlist[i])
lines(pred.p.med ~ dist_seq , lw=2, col="blue" , lty=2)
shade(pred.p.PI , dist_seq, col=col.alpha("blue", alpha=0.05))
title(plantlist[i], line = -1)
##adults
dist_seq <- seq(from=min(d$logdist2d),to=max(d$logdist2d),length=30)

d.pred <- list(
date_index=rep(1,30),
rat_index=rep(1,30),
species_index=rep(i,30),
order_index=rep(1,30),
Adult=rep(1,30),
Distance = dist_seq
)

link2 <- link(m6zag , data=d.pred)
points(grams_eaten~logdist2d, data=subset(d,species_index==i & Adult==1 ), pch=19, col="red")
pred.joint <- link2$p*link2$mu
pred.p.med <- apply(pred.joint , 2 , median)
pred.p.PI <- apply(pred.joint , 2 , PI , prob=0.90)
lines(pred.p.med ~ dist_seq , lw=2, col="red" , lty=1)
shade(pred.p.PI , dist_seq, col=col.alpha("red", alpha=0.05))
legend("topleft", c("Adult", "Juvenile") , lty=c(1,2) , pch=c(19,1), col=c("red","blue"), horiz=FALSE, bty='n')
mtext("Effect of Distance and Age on Diet Preference", outer =TRUE, cex = 1.5 , side=3)
mtext("Grams Eaten", outer =TRUE, cex = 1.2, side=2)
mtext("log distance (m) from edge", outer =TRUE, cex = 1.2, side=1)
par(old.par)


#prob only

plot(consumed~logdist2d,data=subset(d, Adult==0), xlab="log distance from edge (m)" , ylab="Prob Sampling Plant" , ylim=c(0,0.8) , col=col.alpha("white", alpha=0.05), pch=19 , main="Predictions p (1/0) for ZAG Mixture Model")
lines(pred.p.med ~ dist_seq, col="blue" , lw=2 , lty=2)
shade(pred.p.PI , dist_seq , col=col.alpha("blue", alpha=0.05))

plantlist <- c("Poison Oak","Toyon","Bay","Buckeye","Cedar","Chamise","Oak")
dist_seq <- seq(from=min(d$logdist2d),to=max(d$logdist2d),length=30)
old.par <- par(mfrow=c(3, 3),oma = c(2, 2, 2, 0) , cex=0.5)

for (i in 1:6){
    d.pred <- list(
    date_index=rep(1,30),
    rat_index=rep(1,30),
    order_index=rep(1,30),
    species_index=rep(i,30),
    Adult=rep(0,30),
    Distance = dist_seq
    )
    
    link2 <- link(m6zag , data=d.pred)
    
    pred.p.med <- apply(link2$p , 2 , mean)
    pred.p.PI <- apply(link2$p , 2 , PI , prob=0.90)
    plot(consumed~logdist2d,data=subset(d, Adult==0), xlab='' , ylab='', ann=FALSE, ylim=c(0,.8) , col="white", pch=1, main=plantlist[i])
    lines(pred.p.med ~ dist_seq, col="blue" , lw=2 , lty=2)
    shade(pred.p.PI , dist_seq , col=col.alpha("blue", alpha=0.05))
    title(plantlist[i], line = -1.5)
    
    
    ##adults
    dist_seq <- seq(from=min(d$logdist2d),to=max(d$logdist2d),length=30)
    
    d.pred <- list(
    date_index=rep(1,30),
    rat_index=rep(1,30),
    species_index=rep(i,30),
    order_index=rep(1,30),
    Adult=rep(1,30),
    Distance = dist_seq
    )
    
    link2 <- link(m6zag , data=d.pred)
    # points(consumed~logdist2d,data=subset(d, Adult==1), pch=19, col="red")
    pred.p.med <- apply(link2$p , 2 , mean)
    pred.p.PI <- apply(link2$p , 2 , PI , prob=0.90)
    lines(pred.p.med ~ dist_seq , lw=2, col="red" , lty=1)
    shade(pred.p.PI , dist_seq, col=col.alpha("red", alpha=0.05))
    legend("topleft", c("Adult", "Juvenile") , lty=c(1,2) , pch=c(19,1), col=c("red","blue"), horiz=FALSE, bty='n')
}

frame()
i <- 7
d.pred <- list(
date_index=rep(1,30),
rat_index=rep(1,30),
order_index=rep(1,30),
species_index=rep(i,30),
Adult=rep(0,30),
Distance = dist_seq
)

link2 <- link(m6zag , data=d.pred)

pred.p.med <- apply(link2$p , 2 , mean)
pred.p.PI <- apply(link2$p , 2 , PI , prob=0.90)
plot(consumed~logdist2d,data=subset(d,species_index==i & Adult==0 ), xlab='n' , ylab='n' , ann=FALSE, ylim=c(0,.8) , col="white", pch=1, main=plantlist[i])
lines(pred.p.med ~ dist_seq , lw=2, col="blue" , lty=2)
shade(pred.p.PI , dist_seq, col=col.alpha("blue", alpha=0.05))
title(plantlist[i], line = -1)
##adults

d.pred <- list(
date_index=rep(1,30),
rat_index=rep(1,30),
species_index=rep(i,30),
order_index=rep(1,30),
Adult=rep(1,30),
Distance = dist_seq
)

link2 <- link(m6zag , data=d.pred)
#points(consumed~logdist2d, data=subset(d,species_index==i & Adult==1), pch=19, col="red")
pred.p.med <- apply(link2$p , 2 , mean)
pred.p.PI <- apply(link2$p , 2 , PI , prob=0.90)
lines(pred.p.med ~ dist_seq , lw=2, col="red" , lty=1)
shade(pred.p.PI , dist_seq, col=col.alpha("red", alpha=0.05))
legend("topleft", c("Adult", "Juvenile") , lty=c(1,2) , pch=c(19,1), col=c("red","blue"), horiz=FALSE, bty='n')
mtext("Effect of Distance and Age on Diet Preference- P", outer =TRUE, cex = 1.5 , side=3)
mtext("Probabilty of sampling plant", outer =TRUE, cex = 1.2, side=2)
mtext("log distance (m) from edge", outer =TRUE, cex = 1.2, side=1)
par(old.par)

###JOint
plantlist <- c("Poison Oak","Toyon","Bay","Buckeye","Cedar","Chamise","Oak")
dist_seq <- seq(from=min(d$logdist2d),to=max(d$logdist2d),length=30)
old.par <- par(mfrow=c(3, 3),oma = c(2, 2, 2, 0) , cex=0.5)

for (i in 1:6){
    d.pred <- list(
    date_index=rep(1,30),
    rat_index=rep(1,30),
    order_index=rep(1,30),
    species_index=rep(i,30),
    Adult=rep(0,30),
    Distance = dist_seq
    )
    
    link2 <- link(m6zag , data=d.pred)
    pred.joint <- link2$p*link2$mu
    pred.p.med <- apply(pred.joint , 2 , median)
    pred.p.PI <- apply(pred.joint , 2 , PI , prob=0.90)
    plot(grams_eaten~logdist2d,data=subset(d,species_index==i & Adult==0 ), xlab='' , ylab='', ann=FALSE, ylim=c(0,0.8) , col="blue", pch=1, main=plantlist[i])
    lines(pred.p.med ~ dist_seq , lw=2, col="blue" , lty=2)
    shade(pred.p.PI , dist_seq, col=col.alpha("blue", alpha=0.05))
    title(plantlist[i], line = -1)
    
    
    ##adults
    dist_seq <- seq(from=min(d$logdist2d),to=max(d$logdist2d),length=30)
    
    d.pred <- list(
    date_index=rep(1,30),
    rat_index=rep(1,30),
    species_index=rep(i,30),
    order_index=rep(1,30),
    Adult=rep(1,30),
    Distance = dist_seq
    )
    
    link2 <- link(m6zag , data=d.pred)
    points(grams_eaten~logdist2d, data=subset(d,species_index==i & Adult==1 ), pch=19, col="red")
    pred.joint <- link2$p*link2$mu
    pred.p.med <- apply(pred.joint, 2 , median)
    pred.p.PI <- apply(pred.joint , 2 , PI , prob=0.90)
    lines(pred.p.med ~ dist_seq , lw=2, col="red" , lty=1)
    shade(pred.p.PI , dist_seq, col=col.alpha("red", alpha=0.05))
    legend("topleft", c("Adult", "Juvenile") , lty=c(1,2) , pch=c(19,1), col=c("red","blue"), horiz=FALSE, bty='n')
}

frame()
i <- 7
d.pred <- list(
date_index=rep(1,30),
rat_index=rep(1,30),
order_index=rep(1,30),
species_index=rep(i,30),
Adult=rep(0,30),
Distance = dist_seq
)

link2 <- link(m6zag , data=d.pred)
pred.joint <- link2$p*link2$mu

pred.p.med <- apply(pred.joint , 2 , median)
pred.p.PI <- apply(pred.joint , 2 , PI , prob=0.90)
plot(grams_eaten~logdist2d,data=subset(d,species_index==i & Adult==0 ), xlab='n' , ylab='n' , ann=FALSE, ylim=c(0,2) , col="blue", pch=1, main=plantlist[i])
lines(pred.p.med ~ dist_seq , lw=2, col="blue" , lty=2)
shade(pred.p.PI , dist_seq, col=col.alpha("blue", alpha=0.05))
title(plantlist[i], line = -1)
##adults
dist_seq <- seq(from=min(d$logdist2d),to=max(d$logdist2d),length=30)

d.pred <- list(
date_index=rep(1,30),
rat_index=rep(1,30),
species_index=rep(i,30),
order_index=rep(1,30),
Adult=rep(1,30),
Distance = dist_seq
)

link2 <- link(m6zag , data=d.pred)
points(grams_eaten~logdist2d, data=subset(d,species_index==i & Adult==1 ), pch=19, col="red")
pred.joint <- link2$p*link2$mu
pred.p.med <- apply(pred.joint , 2 , median)
pred.p.PI <- apply(pred.joint , 2 , PI , prob=0.90)
lines(pred.p.med ~ dist_seq , lw=2, col="red" , lty=1)
shade(pred.p.PI , dist_seq, col=col.alpha("red", alpha=0.05))
legend("topleft", c("Adult", "Juvenile") , lty=c(1,2) , pch=c(19,1), col=c("red","blue"), horiz=FALSE, bty='n')
mtext("Effect of Distance and Age on Diet Preference", outer =TRUE, cex = 1.5 , side=3)
mtext("Grams Eaten", outer =TRUE, cex = 1.2, side=2)
mtext("log distance (m) from edge", outer =TRUE, cex = 1.2, side=1)
par(old.par)

##plot sex differences
post <- extract.samples(m6g)

#with distance
col_index <- c("red", "orange" , "gold" , "green" , "blue" , "violet" , "black")
mu <- function(distance, species_index=0)
    exp( post$a + post$a_species[,species_index] + (post$b_distance + post$b_distance_species[,species_index])*distance + (post$b_Male_species[,species_index] + post$b_Male)*Male )

log_dist_seq <- seq(from=0, to=5, by=0.1)

old.par <- par(mfrow=c(3, 3))

for (i in 1:7){
plot(grams_eaten~distance, data=subset(d, species_index==i & Male==0), pch=1, col=col_index[i], ylim=c(0,1) ,  xlim=c(0,max(d$distance))  , cex.main=1 , ylab="Grams Eaten" , xlab="Log Distance (m) from Habitat Boundary")

Male <- 0
mu.i <- sapply(log_dist_seq , function(log_dist) median(mu(log_dist,i)))
mu.ci.i <- sapply(log_dist_seq , function(log_dist) HPDI(mu(log_dist,i)))
lines(log_dist_seq , mu.i , lty=5 , col=col_index[i] , lw=2  )
lines(log_dist_seq , mu.ci.i[1,] , lty=3 , col=col_index[i])
lines(log_dist_seq , mu.ci.i[2,] , lty=3 , col=col_index[i])
points(grams_eaten~distance, data=subset(d, species_index==i & Male==1), pch=19, col=col_index[i])

Male <- 1
mu.i <- sapply(log_dist_seq , function(log_dist) median(mu(log_dist,i)))
mu.ci.i <- sapply(log_dist_seq , function(log_dist) HPDI(mu(log_dist,i)))
lines(log_dist_seq , mu.i , lty=1 , col=col_index[i] , lw=2  )
lines(log_dist_seq , mu.ci.i[1,] , lty=3 , col=col_index[i])
lines(log_dist_seq , mu.ci.i[2,] , lty=3 , col=col_index[i])
points(grams_eaten~distance, data=subset(d, species_index==i), pch=1, col=col_index[i])
}
plot(1~1, pch=0 , col="white")
legend("center", inset=.01, title="Plant Species", c("Poison Oak","Toyon" , "Bay" , "Buckeye" , "Cedar" , "chamise" , "Oak"), pch=19, col=col_index, horiz=FALSE)

par(old.par)


colnums <- function(x){
    print( ( cbind( as.list(names(x)) , as.list(1:ncol(x)) ) ) )
            }


######all graphs on 1 window
###MU
plantlist <- c("a. Poison Oak","b. Toyon","c. Bay","d. Buckeye","e. Cedar","f. Chamise","g. Oak" , "h. All Plants")
dist_seq <- seq(from=min(d$logdist2d),to=max(d$logdist2d),length=30)
old.par <- par(mfrow=c(3, 3),oma = c(2, 2, 2, 0) , cex=0.5)

for (i in 1:6){
    d.pred <- list(
    date_index=rep(1,30),
    rat_index=rep(1,30),
    order_index=rep(1,30),
    species_index=rep(i,30),
    Adult=rep(0,30),
    Distance = dist_seq
    )
    
    link2 <- link(m6zag , data=d.pred)
    
    pred.p.med <- apply(link2$mu , 2 , median)
    pred.p.PI <- apply(link2$mu , 2 , PI , prob=0.90)
    plot(grams_eaten~logdist2d,data=subset(d,species_index==i & Adult==0 & grams_eaten > 0), xlab='' , ylab='', ann=FALSE, ylim=c(0,0.8), xlim=c(1,4.2) , col="blue", pch=1, main=plantlist[i])
    lines(pred.p.med ~ dist_seq , lw=2, col="blue" , lty=2)
    shade(pred.p.PI , dist_seq, col=col.alpha("blue", alpha=0.05))
    title(plantlist[i], line = -1)
    
    
    ##adults
    dist_seq <- seq(from=min(d$logdist2d),to=max(d$logdist2d),length=30)
    
    d.pred <- list(
    date_index=rep(1,30),
    rat_index=rep(1,30),
    species_index=rep(i,30),
    order_index=rep(1,30),
    Adult=rep(1,30),
    Distance = dist_seq
    )
    
    link2 <- link(m6zag , data=d.pred)
    points(grams_eaten~logdist2d, data=subset(d,species_index==i & Adult==1 & grams_eaten > 0), pch=19, col="red")
    pred.p.med <- apply(link2$mu , 2 , mean)
    pred.p.PI <- apply(link2$mu , 2 , PI , prob=0.90)
    lines(pred.p.med ~ dist_seq , lw=2, col="red" , lty=1)
    shade(pred.p.PI , dist_seq, col=col.alpha("red", alpha=0.05))
    #legend("topleft", c("Adult", "Juvenile") , lty=c(1,2) , pch=c(19,1), col=c("red","blue"), horiz=FALSE, bty='n')
}

#frame()
i <- 7
d.pred <- list(
date_index=rep(1,30),
rat_index=rep(1,30),
order_index=rep(1,30),
species_index=rep(i,30),
Adult=rep(0,30),
Distance = dist_seq
)

link2 <- link(m6zag , data=d.pred)

pred.p.med <- apply(link2$mu , 2 , mean)
pred.p.PI <- apply(link2$mu , 2 , PI , prob=0.90)
plot(grams_eaten~logdist2d,data=subset(d,species_index==i & Adult==0 & grams_eaten > 0), xlab='n' , ylab='n' , ann=FALSE, ylim=c(0,2), xlim=c(1,4.2) , col="blue", pch=1, main=plantlist[i])
lines(pred.p.med ~ dist_seq , lw=2, col="blue" , lty=2)
shade(pred.p.PI , dist_seq, col=col.alpha("blue", alpha=0.05))
title(plantlist[i], line = -1)
##adults
dist_seq <- seq(from=min(d$logdist2d),to=max(d$logdist2d),length=30)

d.pred <- list(
date_index=rep(1,30),
rat_index=rep(1,30),
species_index=rep(i,30),
order_index=rep(1,30),
Adult=rep(1,30),
Distance = dist_seq
)

link2 <- link(m6zag , data=d.pred)
points(grams_eaten~logdist2d, data=subset(d,species_index==i & Adult==1 & grams_eaten > 0), pch=19, col="red")
pred.p.med <- apply(link2$mu , 2 , mean)
pred.p.PI <- apply(link2$mu , 2 , PI , prob=0.90)
lines(pred.p.med ~ dist_seq , lw=2, col="red" , lty=1)
shade(pred.p.PI , dist_seq, col=col.alpha("red", alpha=0.05))

#all plants

a_species_zeros <- matrix(0,1000,7)
a_rat_zeros <- matrix(0,1000,22)
a_date_zeros <- matrix(0,1000,6)
a_order_zeros <- matrix(0,1000,3)
b_Male_species_zeros <- matrix(0,1000,7)
b_species_Adult_zeros<- matrix(0,1000,7)
b_distance_species_zeros <- matrix(0,1000,7)

d.pred <- list(
date_index=rep(1,30),
rat_index=rep(1,30),
order_index=rep(1,30),
species_index=rep(1,30),
Adult=rep(0,30),
Distance = dist_seq
)

link2 <- link(m6zag , n=1000 , data=d.pred, replace=list(a_ratid=a_rat_zeros, a_species=a_species_zeros,a_date=a_date_zeros, b_species_Adult=b_species_Adult_zeros, b_distance_species=b_distance_species_zeros , a_order=a_order_zeros ), WAIC=TRUE)

#blue juvenile
pred.p.med <- apply(link2$mu , 2 , mean)
pred.p.PI <- apply(link2$mu  , 2 , PI , prob=0.90)
plot(grams_eaten~logdist2d,data=subset(d, Adult==0 & grams_eaten>0), xlab='n' , ylab=
'n' , ann=FALSE, ylim=c(0,2), xlim=c(1,4.2) , col="blue", pch=1 , main="h. All Plants")
lines(pred.p.med ~ dist_seq , lw=2, col="blue" , lty=2)
shade(pred.p.PI , dist_seq, col=col.alpha("blue", alpha=0.05))
title(plantlist[8], line = -1)



d.pred <- list(
date_index=rep(1,30),
rat_index=rep(1,30),
species_index=rep(i,30),
order_index=rep(1,30),
Adult=rep(1,30),
Distance = dist_seq
)

link2 <- link(m6zag , n=1000 , data=d.pred, replace=list(a_ratid=a_rat_zeros, a_species=a_species_zeros,a_date=a_date_zeros, b_species_Adult=b_species_Adult_zeros, b_distance_species=b_distance_species_zeros ), WAIC=TRUE)
points(grams_eaten~logdist2d, data=subset(d, Adult==1 & grams_eaten>0), pch=19, col="red")
pred.p.med <- apply(link2$mu , 2 , median)
pred.p.PI <- apply(link2$mu , 2 , PI , prob=0.90)
lines(pred.p.med ~ dist_seq , lw=2, col="red" )
shade(pred.p.PI , dist_seq, col=col.alpha("red", alpha=0.05))

#legend
plot(2,2 , xlab='n' , ylab='n' , ann=FALSE, ylim=c(0,1) , xlim=c(0,1) , col="white", pch=0 ,axes=FALSE)
legend("center", c("Adult", "Juvenile") , lty=c(1,2) , pch=c(19,1), col=c("red","blue"), horiz=FALSE, bty='n', cex=2.5)
mtext(expression("Figure 2. Model averaged predictions for" ~mu), outer =TRUE, cex = 1.2 , side=3)
mtext("grams eaten", outer =TRUE, cex = 1.2, side=2)
mtext("log distance (m) from edge", outer =TRUE, cex = 1.2, side=1)
par(old.par)

###P

plantlist <- c("a. Poison Oak","b. Toyon","c. Bay","d. Buckeye","e. Cedar","f. Chamise","g. Oak" , "h. All Plants")
dist_seq <- seq(from=min(d$logdist2d),to=max(d$logdist2d),length=30)
old.par <- par(mfrow=c(3, 3),oma = c(2, 2, 2, 0) , cex=0.5)

for (i in 1:6){
    d.pred <- list(
    date_index=rep(1,30),
    rat_index=rep(1,30),
    order_index=rep(1,30),
    species_index=rep(i,30),
    Adult=rep(0,30),
    Distance = dist_seq
    )
    
    link2 <- link(m6zag , data=d.pred)
    
    pred.p.med <- apply(link2$p , 2 , mean)
    pred.p.PI <- apply(link2$p , 2 , PI , prob=0.90)
    plot(consumed~logdist2d,data=subset(d, Adult==0), xlab='' , ylab='', ann=FALSE, ylim=c(0,.8) , xlim=c(1,4.2) , col="white", pch=1, main=plantlist[i])
    lines(pred.p.med ~ dist_seq, col="blue" , lw=2 , lty=2)
    shade(pred.p.PI , dist_seq , col=col.alpha("blue", alpha=0.05))
    title(plantlist[i], line = -1.5)
    
    
    ##adults
    dist_seq <- seq(from=min(d$logdist2d),to=max(d$logdist2d),length=30)
    
    d.pred <- list(
    date_index=rep(1,30),
    rat_index=rep(1,30),
    species_index=rep(i,30),
    order_index=rep(1,30),
    Adult=rep(1,30),
    Distance = dist_seq
    )
    
    link2 <- link(m6zag , data=d.pred)
    # points(consumed~logdist2d,data=subset(d, Adult==1), pch=19, col="red")
    pred.p.med <- apply(link2$p , 2 , mean)
    pred.p.PI <- apply(link2$p , 2 , PI , prob=0.90)
    lines(pred.p.med ~ dist_seq , lw=2, col="red" , lty=1)
    shade(pred.p.PI , dist_seq, col=col.alpha("red", alpha=0.05))
    #legend("topleft", c("Adult", "Juvenile") , lty=c(1,2) , pch=c(19,1), col=c("red","blue"), horiz=FALSE, bty='n')
}

#frame()
i <- 7
d.pred <- list(
date_index=rep(1,30),
rat_index=rep(1,30),
order_index=rep(1,30),
species_index=rep(i,30),
Adult=rep(0,30),
Distance = dist_seq
)

link2 <- link(m6zag , data=d.pred)

pred.p.med <- apply(link2$p , 2 , mean)
pred.p.PI <- apply(link2$p , 2 , PI , prob=0.90)
plot(consumed~logdist2d,data=subset(d,species_index==i & Adult==0 ), xlab='n' , ylab='n' , ann=FALSE, ylim=c(0,.8), xlim=c(1,4.2) , col="white", pch=1, main=plantlist[i])
lines(pred.p.med ~ dist_seq , lw=2, col="blue" , lty=2)
shade(pred.p.PI , dist_seq, col=col.alpha("blue", alpha=0.05))
title(plantlist[i], line = -1)
##adults

d.pred <- list(
date_index=rep(1,30),
rat_index=rep(1,30),
species_index=rep(i,30),
order_index=rep(1,30),
Adult=rep(1,30),
Distance = dist_seq
)

link2 <- link(m6zag , data=d.pred)
#points(consumed~logdist2d, data=subset(d,species_index==i & Adult==1), pch=19, col="red")
pred.p.med <- apply(link2$p , 2 , mean)
pred.p.PI <- apply(link2$p , 2 , PI , prob=0.90)
lines(pred.p.med ~ dist_seq , lw=2, col="red" , lty=1)
shade(pred.p.PI , dist_seq, col=col.alpha("red", alpha=0.05))
#legend("topleft", c("Adult", "Juvenile") , lty=c(1,2) , pch=c(19,1), col=c("red","blue"), horiz=FALSE, bty='n')


##
#all plants

a_species_zeros <- matrix(0,1000,7)
a_rat_zeros <- matrix(0,1000,22)
a_date_zeros <- matrix(0,1000,6)
a_order_zeros <- matrix(0,1000,3)
b_Male_species_zeros <- matrix(0,1000,7)
b_species_Adult_zeros<- matrix(0,1000,7)
b_distance_species_zeros <- matrix(0,1000,7)

d.pred <- list(
date_index=rep(1,30),
rat_index=rep(1,30),
order_index=rep(1,30),
species_index=rep(1,30),
Adult=rep(0,30),
Distance = dist_seq
)

link2 <- link(m6zag , n=1000 , data=d.pred, replace=list(a_ratid=a_rat_zeros, a_species=a_species_zeros,a_date=a_date_zeros, b_species_Adult=b_species_Adult_zeros, b_distance_species=b_distance_species_zeros , a_order=a_order_zeros ), WAIC=TRUE)

#blue juvenile
pred.p.med <- apply(link2$p , 2 , mean)
pred.p.PI <- apply(link2$p  , 2 , PI , prob=0.90)
plot(grams_eaten~logdist2d,data=subset(d, Adult==0 & grams_eaten>0), xlab='n' , ylab=
'n' , ann=FALSE, ylim=c(0,0.8), xlim=c(1,4.2) , col="white", pch=1 , main="h. All Plants")
lines(pred.p.med ~ dist_seq , lw=2, col="blue" , lty=2)
shade(pred.p.PI , dist_seq, col=col.alpha("blue", alpha=0.05))
title(plantlist[8], line = -1)



d.pred <- list(
date_index=rep(1,30),
rat_index=rep(1,30),
species_index=rep(1,30),
order_index=rep(1,30),
Adult=rep(1,30),
Distance = dist_seq
)

link2 <- link(m6zag , n=1000 , data=d.pred, replace=list(a_ratid=a_rat_zeros, a_species=a_species_zeros,a_date=a_date_zeros, b_species_Adult=b_species_Adult_zeros, b_distance_species=b_distance_species_zeros ), WAIC=TRUE)
#points(grams_eaten~logdist2d, data=subset(d, Adult==1 & grams_eaten>0), pch=19, col="red")
pred.p.med <- apply(link2$p , 2 , mean)
pred.p.PI <- apply(link2$p , 2 , PI , prob=0.90)
lines(pred.p.med ~ dist_seq , lw=2, col="red" )
shade(pred.p.PI , dist_seq, col=col.alpha("red", alpha=0.05))

#legend
plot(2,2 , xlab='n' , ylab='n' , ann=FALSE, ylim=c(0,1) , xlim=c(0,1) , col="white", pch=0 ,axes=FALSE)
legend("center", c("Adult", "Juvenile") , lty=c(1,2) , col=c("red","blue"), horiz=FALSE, bty='n', cex=2.5)
mtext(expression("Figure 1. Model averaged predictions for" ~p), outer =TRUE, cex = 1.2 , side=3)
###
mtext("probabilty of sampling plant", outer =TRUE, cex = 1.2, side=2)
mtext("log distance (m) from edge", outer =TRUE, cex = 1.2, side=1)
par(old.par)




