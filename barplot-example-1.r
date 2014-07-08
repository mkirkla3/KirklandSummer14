#092204

 lab <- read.table( "lab.tab", header=T );
 nat <- read.table( "092104.rls.natural-isolates.results.summary.tab", header=T );

# my.breaks = c(13,15,17,19,21,23,25,27,29,31,33,35,37, 39, 41);
 my.breaks = seq( 12.5,36.5,by=3);
 
 h.nat <- hist(nat$rmean, br= my.breaks, xlab = "replicative lifespan", ylab = "relative density", freq=F ) ;

 h.lab <- hist(lab$rmean, br= my.breaks, xlab = "replicative lifespan",  ylab = "relative density", freq=F ) ;
   
 #generate the comparison table
 tb <-  data.frame( rbind(h.nat$density,h.lab$density) )  ;
 names( tb ) <- h.nat$mids;
 row.names(tb) <- c( "natural", "lab" )
 tb
 
 write.table( tb, file="binned-comparison-matrix.csv", sep="\t", col.names=T, quote=T, row.names=T );

 
 postscript( "092104.natural.vs.lab.barplot.ps" );
 #barplot( rbind(h.nat$density,h.lab$density), beside=T, col=c("red","blue"));
 barplot( as.matrix(tb), beside=T, col=c("red","blue"), ylab="Relative Density", xlab="Replicative Lifespan",
     legend= c( "Natural isolates", "Lab strains" ),
 );
 title(main="Comparison of lifespan distributions" )
 dev.off(); 
 
 postscript( "092104.natural.vs.lab.overlay.ps" );
 hist(lab$rmean, br= my.breaks, xlab = "replicative lifespan", ylab = "relative density", freq=F,
  main = "Comparison of natural and lab strains lifespan distributions", 
 );
 lines( h.nat$mids, h.nat$density, col="red"); 
 dev.off();
 
 

 
 
# q("no");



