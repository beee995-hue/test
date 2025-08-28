##Q1
pop1<-as.matrix(rbind(c(1,2),c(1,1),c(2,1),c(2,2),c(2,2)))

pop1[3,]

pop1[4,1]

pop1[4,2]

nrow(pop1)

##Q2

indivX<-pop1[3,]

test <- if (indivX[1] == indivX[2]) {
  print("TRUE")
} else {
  print("FALSE")
}


indivX<-pop1[4,]

test <- if (indivX[1] == indivX[2]) {
  print("TRUE")
} else {
  print("FALSE")
}

isHomozygote <- function(indivX) {
  if (indivX[1] == indivX[2]) {
    return("TRUE")
  } else {
    return("FALSE")
  }
}


indiv1 <- c("A", "A")
indiv2 <- c("A", "T")

isHomozygote(indiv1)  # returns "TRUE"
isHomozygote(indiv2)  # returns "FALSE"


##Q3


results <- apply(pop1, 1, isHomozygote)
results



# Count HZ
nbHomozygoteIndividuals <- sum(results)

nbHomozygoteIndividuals


nbHomozygotes <- function(pop) {
  # Apply isHomozygote across rows, sum the TRUEs
  sum(apply(pop, 1, function(indiv) indiv[1] == indiv[2]))
}

nbHomozygotes(pop1)




##Q4 


countGenotypes <- function(pop) {
  # Check which rows are homozygotes
  isHomo <- apply(pop, 1, function(indiv) indiv[1] == indiv[2])
  
  # Homozygotes: collapse alleles (e.g., "AA", "GG")
  homozygotes <- apply(pop[isHomo, , drop = FALSE], 1, paste, collapse = "")
  
  # Count unique homozygote types
  homozygoteCounts <- table(homozygotes)
  
  # Heterozygotes: just count them all together
  heterozygoteCount <- sum(!isHomo)
  
  return(list(
    homozygotes = homozygoteCounts,
    heterozygotes = heterozygoteCount
  ))
}


countGenotypes(pop1)



##Q5


countAlleles <- function(pop) {
  # Flatten the matrix into a single vector of alleles
  alleles <- as.vector(pop)
  
  # Count occurrences of each allele
  alleleCounts <- table(alleles)
  
  
  return(alleleCounts)
}

countAlleles(pop1)



##Q6


calcHe <- function(pop) {
  # Flatten alleles
  alleles <- as.vector(pop)
  
  # Count alleles
  alleleCounts <- table(alleles)
  
  # Frequencies
  freqs <- alleleCounts / sum(alleleCounts)
  
  # Expected heterozygosity
  He <- 1 - sum(freqs^2)
  
  return(He)
}

calcHe(pop1)


##Q7

calcFis <- function(pop) {
  # Number of individuals
  nInd <- nrow(pop)
  
  # Observed heterozygosity (Ho)
  Ho <- sum(pop[,1] != pop[,2]) / nInd
  
  # Expected heterozygosity (He)
  alleles <- as.vector(pop)
  alleleCounts <- table(alleles)
  freqs <- alleleCounts / sum(alleleCounts)
  He <- 1 - sum(freqs^2)
  
  # Inbreeding coefficient Fis
  Fis <- (He - Ho) / He
  
  return(Fis)
}


calcFis(pop1)


##Q8
nextGen <- function(pop) {
  nOffspring <- nrow(pop)  # same number as parent population
  
  # Flatten the population to get all alleles
  alleles <- as.vector(pop)
  
  # Compute allele frequencies
  alleleFreqs <- table(alleles) / length(alleles)
  alleleNames <- names(alleleFreqs)
  
  # Draw alleles for each offspring
  offspring <- matrix(nrow = nOffspring, ncol = 2)
  
  for (i in 1:nOffspring) {
    offspring[i, 1] <- sample(alleleNames, size = 1, prob = alleleFreqs)
    offspring[i, 2] <- sample(alleleNames, size = 1, prob = alleleFreqs)
  }
  
  return(offspring)
}


nextGen(pop1)


##Q9


popEvol<- function(pop, nGen) {
  currentPop <- pop
  nInd <- nrow(pop)
  
  for (gen in 1:nGen) {
    alleles <- as.vector(currentPop)
    alleleFreqs <- table(alleles) / length(alleles)
    alleleNames <- names(alleleFreqs)
    
    nextPop <- matrix(nrow = nInd, ncol = 2)
    for (i in 1:nInd) {
      nextPop[i, 1] <- sample(alleleNames, size = 1, prob = alleleFreqs)
      nextPop[i, 2] <- sample(alleleNames, size = 1, prob = alleleFreqs)
    }
    
    currentPop <- nextPop
  }
  
  # Compute allele frequencies in last generation
  finalAlleles <- as.vector(currentPop)
  finalAlleleFreqs <- table(finalAlleles) / length(finalAlleles)
  
  # Return as list
  return(list(
    finalPop = currentPop,
    alleleFreqs = finalAlleleFreqs
  ))
}

##Q10


popEvol(pop1,30)




# Number of simulations and generations
nSim <- 40
nGen <- 30
A1_freqs <- numeric(nSim)

set.seed(123)
for (i in 1:nSim) {
  result <- popEvol(pop1, nGen)
  finalAlleleFreqs <- result$alleleFreqs
  # Store frequency of allele 1
  A1_freqs[i] <- ifelse("1" %in% names(finalAlleleFreqs),
                        finalAlleleFreqs["1"],
                        0)
}

A1_freqs

mean(A1_freqs)
var(A1_freqs)


##Q11

nextGenA1 <- function(freqs, nInd) {
  # freqs: named vector of allele frequencies, e.g., c("1" = 0.6, "2" = 0.4)
  # nInd: number of individuals in the next generation
  
  alleles <- names(freqs)
  
  # Generate offspring: two alleles per individual
  offspring <- matrix(nrow = nInd, ncol = 2)
  for (i in 1:nInd) {
    offspring[i, 1] <- sample(alleles, size = 1, prob = freqs)
    offspring[i, 2] <- sample(alleles, size = 1, prob = freqs)
  }
  
  # Compute frequency of allele 1 in the new generation
  allAlleles <- as.vector(offspring)
  freqA1 <- sum(allAlleles == "1") / length(allAlleles)
  
  return(freqA1)
}



nextGenA1 <- function(freqA1, nInd) {
  # freqA1: current frequency of allele 1
  # nInd: population size (number of individuals)
  
  # Number of allele 1 copies in next generation
  nA1 <- rbinom(1, size = 2 * nInd, prob = freqA1)
  
  # Frequency of A1 in next generation
  freqNext <- nA1 / (2 * nInd)
  
  return(freqNext)
}



## Q12

# Parameters
nSim <- 40        # number of replicates
nGen <- 200       # number of generations
nInd <- 50        # population size
initFreq <- 0.4   # initial frequency of allele 1

# Vector to store final frequencies
finalFreqs <- numeric(nSim)

set.seed(123)  # for reproducibility

for (sim in 1:nSim) {
  freq <- initFreq
  for (gen in 1:nGen) {
    freq <- nextGenA1(freq, nInd)  # update allele frequency each generation
  }
  finalFreqs[sim] <- freq
}

# Inspect the resulting vector
finalFreqs


mean(finalFreqs)
var(finalFreqs)

##500

# Parameters
nSim <- 40        # number of replicates
nGen <- 200       # number of generations
nInd <- 500        # population size
initFreq <- 0.4   # initial frequency of allele 1

# Vector to store final frequencies
finalFreqs <- numeric(nSim)

set.seed(123)  # for reproducibility

for (sim in 1:nSim) {
  freq <- initFreq
  for (gen in 1:nGen) {
    freq <- nextGenA1(freq, nInd)  # update allele frequency each generation
  }
  finalFreqs[sim] <- freq
}

# Inspect the resulting vector
finalFreqs


mean(finalFreqs)
var(finalFreqs)

##5000

# Parameters
nSim <- 40        # number of replicates
nGen <- 200       # number of generations
nInd <- 5000       # population size
initFreq <- 0.4   # initial frequency of allele 1

# Vector to store final frequencies
finalFreqs <- numeric(nSim)

set.seed(123)  # for reproducibility

for (sim in 1:nSim) {
  freq <- initFreq
  for (gen in 1:nGen) {
    freq <- nextGenA1(freq, nInd)  # update allele frequency each generation
  }
  finalFreqs[sim] <- freq
}

# Inspect the resulting vector
finalFreqs


mean(finalFreqs)
var(finalFreqs)









##50000

# Parameters
nSim <- 40        # number of replicates
nGen <- 200       # number of generations
nInd <- 50000        # population size
initFreq <- 0.4   # initial frequency of allele 1

# Vector to store final frequencies
finalFreqs <- numeric(nSim)

set.seed(123)  # for reproducibility

for (sim in 1:nSim) {
  freq <- initFreq
  for (gen in 1:nGen) {
    freq <- nextGenA1(freq, nInd)  # update allele frequency each generation
  }
  finalFreqs[sim] <- freq
}

# Inspect the resulting vector
finalFreqs


mean(finalFreqs)
var(finalFreqs)






