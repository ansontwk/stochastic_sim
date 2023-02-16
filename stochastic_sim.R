#============================================

#variables
c_transcription <- 0.005
c_translation <- 0.167
c_mrnaloss <- log(2)/120
c_proteinloss <- log(2)/3600

cycle_count <- 30000

#initial values
t <- 0 #time
mRNA <- 0
Protein <- 0
DNA <- 1

#Edit here
#=============================================
#storage
time_store <- c(0)
mRNA_store <- c(mRNA)
protein_store <- c(Protein)

#main code
counter <- 1
while (counter <= cycle_count){
  
  a1 <- c_transcription * DNA
  a2 <- c_translation * mRNA
  a3 <- c_mrnaloss * mRNA
  a4 <- c_proteinloss * Protein
  
  a_total <- a1 + a2 + a3 + a4
  r1 <- runif(1)
  tau <- (1/a_total) * log(1/r1)
  t <- t + tau
  
  r2 <- runif(1)
  comparison = a_total * r2
  
  if (comparison < a1) {
    mRNA <- mRNA +1
    
  } else if (comparison < a1+a2) {
    Protein <- Protein + 1
    
  } else if (comparison < a1+a2+a3) {
    mRNA <- mRNA - 1
    
  } else if (comparison < a_total) {
    Protein <- Protein - 1
    
  }
  
  mRNA_store <- c(mRNA_store, mRNA)
  protein_store <-c(protein_store, Protein)
  time_store <- c(time_store, t)
  counter <- counter + 1}


#Plotting
plot(time_store, protein_store, xlab = 'Time(s)', ylab = 'Protein number', pch = 20)
points(time_store, mRNA_store, pch = 20 )
axis(side=4)
mtext('mRNA number', side = 4)
title(main = 'Stochastic Simulation of Protein and mRNA')

