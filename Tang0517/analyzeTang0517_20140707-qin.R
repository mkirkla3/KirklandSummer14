rm( list=ls())
list.files()

file = "AgingData_0517.csv"
tb = read.csv( file );
head(tb)
summary(tb)

t.test(tb$wt, tb$wt.CR, alt='less')#this group CR effect is not significant, p 0.07
t.test(tb$ybr053cD, tb$Ybr053c.CR, alt='gr') #p 0.029
t.test(tb$izh4D, tb$izh4.CR) #p 0.003
t.test(tb$Pvac8ERG2, tb$Pvac8ERG2.CR)


##### log likelihood function, simple gompertz mortality model
   #s = exp( (I/G) *(1 - exp(G* my.data)) )  ;
   #m = I * exp( G * my.data );   
 llh.gompertz.single.run <- function( IG, lifespan, trace=0 ) {
   #print(IG)
   I = IG[1]; G = IG[2]; 
   my.data = lifespan[!is.na(lifespan)];
   log_s = (I/G) *(1 - exp(G* my.data))
   if( I< 0 ) { I = 1E-10 }
   log_m = log(I) +  G * my.data ; 
   my.lh = sum(log_s)  + sum(log_m);
   if (trace) {
     print (IG ); #trace the convergence
   }
   ret = - my.lh # because optim seems to maximize 
 }

#####

#ret1a = optim ( c(0.0034, 0.1), fn=llh.gompertz.single.run, lifespan=tb[,1], lower=c(1E-10, 1E-5), method="L-BFGS-B" );
#ret1b = optim ( c(0.05, 0.2), fn=llh.gompertz.single.run, lifespan=tb[,1], lower=c(1E-10, 1E-5), method="L-BFGS-B" );
#ret1c = optim ( c(0.039, 0.1), fn=llh.gompertz.single.run, lifespan=tb[,1]  );
#ret1d = optim ( c(0.05, 0.15), fn=llh.gompertz.single.run, lifespan=tb[,1]  );
#ret1c$par / ret1d$par
#ret1e = optim ( c(0.1, 0.1), fn=llh.gompertz.single.run, lifespan=tb[,1]  );
#ret1e$par / ret1d$par

#set up a template for fitting parameters
label = c("R","G", "mean", "median", "stddev", "CV", 'p_wilcox', 'p_ks')
allFitPar = data.frame( matrix(NA, nrow=length(label), ncol=length(tb[1,])) )
rownames(allFitPar) = label
names(allFitPar) = names(tb)

  
#loop over 
for (col in names(tb)) {
  cur_lifespan = tb[,col]
  retTmp = optim ( c(0.039, 0.1), fn=llh.gompertz.single.run, lifespan=cur_lifespan  );
  allFitPar[1:2, col] = retTmp$par
  allFitPar["mean", col] = mean(cur_lifespan, na.rm=T)
  allFitPar["median", col] = median(cur_lifespan, na.rm=T)
  allFitPar["stddev", col] = sqrt( var(cur_lifespan, na.rm=T))
  allFitPar["CV", col] = allFitPar["stddev", col] / allFitPar["mean", col]
  tmp = wilcox.test( tb[,'wt'], tb[,col] )  
  allFitPar['p_wilcox', col] = tmp$p.value
  tmp2 = ks.test( tb[,'wt'], tb[,col] )  
  allFitPar['p_ks', col] = tmp2$p.value
  
}
allFitPar


#################################
##### LRT for wt, wtCR
### 1) model H0, same G and same I
### 2) model H1i, same G, different I
### 3) model H1g, different G, same I
### 4) model H2, different G, different I
                                  #    I1       G1       I2         G2
 H0  <- function( rawIG ) { IG <- c(rawIG[1], rawIG[2], rawIG[1], rawIG[2] ) }  #all the same
 H1i <- function( rawIG ) { IG <- c(rawIG[1], rawIG[2], rawIG[3], rawIG[2] ) }  #different I
 H1g <- function( rawIG ) { IG <- c(rawIG[1], rawIG[2], rawIG[1], rawIG[4] ) }  # different G
 H2  <- function( rawIG ) { IG <- c(rawIG[1], rawIG[2], rawIG[3], rawIG[4] ) }  # all different

 Hx.llh.gompertz.single.run <- function( rawIG, model, lifespan1, lifespan2 ) {
   IG = model(rawIG); 
   I1 = IG[1]; G1 = IG[2]; I2 = IG[3]; G2 = IG[4];
   my.data1 = lifespan1[!is.na(lifespan1)];
   my.data2 = lifespan2[!is.na(lifespan2)];
   log_s1 = (I1/G1) *(1 - exp(G1* my.data1))
   log_s2 = (I2/G2) *(1 - exp(G2* my.data2))
   log_m1 = log(I1) +  G1 * my.data1 ; 
   log_m2 = log(I2) +  G2 * my.data2 ; 
   my.lh = sum(log_s1) + sum(log_m1) + sum(log_s2) + sum(log_m2)
   #print (IG ); #trace the convergence
   ret = - my.lh # because optim seems to maximize 
 }

## LRT to exam whether wt, wtCR share the same G, I 
llh.H0   = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H0,   lifespan1=tb$wt, lifespan2=tb$wtCR );
llh.H1i  = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H1i,  lifespan1=tb$wt, lifespan2=tb$wtCR );
llh.H1g  = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H1g,  lifespan1=tb$wt, lifespan2=tb$wtCR );
llh.H2  = optim( c(0.01,0.2,0.01,0.1),   Hx.llh.gompertz.single.run, model=H2,  lifespan1=tb$wt, lifespan2=tb$wtCR );

cbind(llh.H0$par, llh.H1i$par, llh.H1g$par, llh.H2$par);

LH = c(llh.H0$value, llh.H1i$value, llh.H1g$value, llh.H2$value);
LH
deltaLH =  - LH + llh.H0$value; 
deltaLH
1 - pchisq( 2*deltaLH, df =c(1,1,1,2) );
#p=1 not significant for all models. 

#################################
###### End of LRT    ############
#################################

## LRT to exam whether ybr053cD, and CR share the same G, I 
llh.H0   = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H0,   lifespan1=tb[,3], lifespan2=tb[,4] );
llh.H1i  = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H1i,  lifespan1=tb[,3], lifespan2=tb[,4] );
llh.H1g  = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H1g,  lifespan1=tb[,3], lifespan2=tb[,4] );
llh.H2   = optim( c(0.01,0.2,0.01,0.1),   Hx.llh.gompertz.single.run, model=H2,  lifespan1=tb[,3], lifespan2=tb[,4] );

cbind(llh.H0$par, llh.H1i$par, llh.H1g$par, llh.H2$par);

LH = c(llh.H0$value, llh.H1i$value, llh.H1g$value, llh.H2$value);
LH
deltaLH =  - LH + llh.H0$value; 
deltaLH
1 - pchisq( 2*deltaLH, df =c(1,1,1,2) );
#[1] 1.00000000 0.01498437 0.01368829 0.04515191

#################################
###### End of LRT    ############
#################################

## LRT to exam whether izh4 and CR share the same G, I 
llh.H0   = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H0,   lifespan1=tb[,5], lifespan2=tb[,6] );
llh.H1i  = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H1i,  lifespan1=tb[,5], lifespan2=tb[,6] );
llh.H1g  = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H1g,  lifespan1=tb[,5], lifespan2=tb[,6] );
llh.H2   = optim( c(0.01,0.2,0.01,0.1),   Hx.llh.gompertz.single.run, model=H2,  lifespan1=tb[,5], lifespan2=tb[,6] );
cbind(llh.H0$par, llh.H1i$par, llh.H1g$par, llh.H2$par);

LH = c(llh.H0$value, llh.H1i$value, llh.H1g$value, llh.H2$value);
LH
# [1] 233.5728 231.6154 233.1406 230.0525
deltaLH =  - LH + llh.H0$value; 
deltaLH
1 - pchisq( 2*deltaLH, df =c(1,1,1,2) );
#[1] 1.00000000 0.04786349 0.35252235 0.02958972

#################################
###### End of LRT    ############
#################################

## LRT to exam whether Pvac8ERG2 and CR share the same G, I 
llh.H0   = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H0,   lifespan1=tb[,7], lifespan2=tb[,8],
                  lower=c(1E-10, 1E-5, 1E-5, 1E-5), method="L-BFGS-B" );
llh.H1i  = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H1i,  lifespan1=tb[,7], lifespan2=tb[,8] );
                  lower=c(1E-10, 1E-5, 1E-5, 1E-5), method="L-BFGS-B" );
llh.H1g  = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H1g,  lifespan1=tb[,7], lifespan2=tb[,8] );
                  lower=c(1E-10, 1E-5, 1E-5, 1E-5), method="L-BFGS-B" );
llh.H2   = optim( c(0.01,0.2,0.01,0.1),   Hx.llh.gompertz.single.run, model=H2,  lifespan1=tb[,7], lifespan2=tb[,8] );
                  lower=c(1E-10, 1E-5, 1E-5, 1E-5), method="L-BFGS-B" );
cbind(llh.H0$par, llh.H1i$par, llh.H1g$par, llh.H2$par);
LH = c(llh.H0$value, llh.H1i$value, llh.H1g$value, llh.H2$value);
LH
#[1] 213.8032 209.7650 209.6608 209.5800
deltaLH =  - LH + llh.H0$value; 
deltaLH
1 - pchisq( 2*deltaLH, df =c(1,1,1,2) );
#[1] 1.000000000 0.004484780 0.003998046 0.014652552

#################################
###### End of LRT    ############
#################################


###############LRT compare mutants with wt (1), 3, 5, 7
llh.H0   = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H0,   lifespan1=tb[,1], lifespan2=tb[,3] );
llh.H1i  = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H1i,  lifespan1=tb[,1], lifespan2=tb[,3] );
llh.H1g  = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H1g,  lifespan1=tb[,1], lifespan2=tb[,3] );
llh.H2   = optim( c(0.01,0.2,0.01,0.1),   Hx.llh.gompertz.single.run, model=H2,  lifespan1=tb[,1], lifespan2=tb[,3] );
cbind(llh.H0$par, llh.H1i$par, llh.H1g$par, llh.H2$par);
LH = c(llh.H0$value, llh.H1i$value, llh.H1g$value, llh.H2$value);
LH
deltaLH =  - LH + llh.H0$value; 
deltaLH
1 - pchisq( 2*deltaLH, df =c(1,1,1,2) );

###############LRT compare mutants with wt (1), 3, 5, 7
llh.H0   = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H0,   lifespan1=tb[,1], lifespan2=tb[,5] );
llh.H1i  = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H1i,  lifespan1=tb[,1], lifespan2=tb[,5] );
llh.H1g  = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H1g,  lifespan1=tb[,1], lifespan2=tb[,5] );
llh.H2   = optim( c(0.01,0.2,0.01,0.1),   Hx.llh.gompertz.single.run, model=H2,  lifespan1=tb[,1], lifespan2=tb[,5] );
cbind(llh.H0$par, llh.H1i$par, llh.H1g$par, llh.H2$par);
LH = c(llh.H0$value, llh.H1i$value, llh.H1g$value, llh.H2$value);
LH
deltaLH =  - LH + llh.H0$value; 
deltaLH
1 - pchisq( 2*deltaLH, df =c(1,1,1,2) );

###############LRT compare mutants with wt (1), 7, significant R and G
llh.H0   = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H0,   lifespan1=tb[,1], lifespan2=tb[,7] );
llh.H1i  = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H1i,  lifespan1=tb[,1], lifespan2=tb[,7] );
llh.H1g  = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H1g,  lifespan1=tb[,1], lifespan2=tb[,7] );
llh.H2   = optim( c(0.01,0.2,0.01,0.1),   Hx.llh.gompertz.single.run, model=H2,  lifespan1=tb[,1], lifespan2=tb[,7] );
cbind(llh.H0$par, llh.H1i$par, llh.H1g$par, llh.H2$par);
LH = c(llh.H0$value, llh.H1i$value, llh.H1g$value, llh.H2$value);
LH
deltaLH =  - LH + llh.H0$value; 
deltaLH
1 - pchisq( 2*deltaLH, df =c(1,1,1,2) );


###############LRT compare mutants with wtCR (2), 4, significant R and G
llh.H0   = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H0,   lifespan1=tb[,2], lifespan2=tb[,4] );
llh.H1i  = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H1i,  lifespan1=tb[,2], lifespan2=tb[,4] );
llh.H1g  = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H1g,  lifespan1=tb[,2], lifespan2=tb[,4] );
llh.H2   = optim( c(0.01,0.2,0.01,0.1),   Hx.llh.gompertz.single.run, model=H2,  lifespan1=tb[,2], lifespan2=tb[,4] );
cbind(llh.H0$par, llh.H1i$par, llh.H1g$par, llh.H2$par);
LH = c(llh.H0$value, llh.H1i$value, llh.H1g$value, llh.H2$value);
LH
deltaLH =  - LH + llh.H0$value; 
deltaLH
1 - pchisq( 2*deltaLH, df =c(1,1,1,2) );

###############LRT compare mutants with wtCR (2), 6, significant R and G
llh.H0   = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H0,   lifespan1=tb[,2], lifespan2=tb[,6] );
llh.H1i  = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H1i,  lifespan1=tb[,2], lifespan2=tb[,6] );
llh.H1g  = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H1g,  lifespan1=tb[,2], lifespan2=tb[,6] );
llh.H2   = optim( c(0.01,0.2,0.01,0.1),   Hx.llh.gompertz.single.run, model=H2,  lifespan1=tb[,2], lifespan2=tb[,6] );
cbind(llh.H0$par, llh.H1i$par, llh.H1g$par, llh.H2$par);
LH = c(llh.H0$value, llh.H1i$value, llh.H1g$value, llh.H2$value);
LH
deltaLH =  - LH + llh.H0$value; 
deltaLH
1 - pchisq( 2*deltaLH, df =c(1,1,1,2) );


###############LRT compare mutants with wtCR (2), 8, significant R or G?
llh.H0   = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H0,   lifespan1=tb[,2], lifespan2=tb[,8] );
llh.H1i  = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H1i,  lifespan1=tb[,2], lifespan2=tb[,8] );
llh.H1g  = optim( c(0.01,0.2,0.01,0.1),  Hx.llh.gompertz.single.run, model=H1g,  lifespan1=tb[,2], lifespan2=tb[,8] );
llh.H2   = optim( c(0.01,0.2,0.01,0.1),   Hx.llh.gompertz.single.run, model=H2,  lifespan1=tb[,2], lifespan2=tb[,8] );
cbind(llh.H0$par, llh.H1i$par, llh.H1g$par, llh.H2$par);
LH = c(llh.H0$value, llh.H1i$value, llh.H1g$value, llh.H2$value);
LH
deltaLH =  - LH + llh.H0$value; 
deltaLH
1 - pchisq( 2*deltaLH, df =c(1,1,1,2) );
