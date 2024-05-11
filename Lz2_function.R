## ------------------------------------------------------------------------
## This function is based on multilevel 2PL model
## ------------------------------------------------------------------------ 

## Dat: response pattern
   # Dimension: number of test-takers * number of items
      
## Group: group index matrix
   # Dimension: number of test takers * 1
    
## Ability: latent abilities
   # Ability[, 1]: group level ability (fix)
   # Dimension: number of test-takers * 1
    
## IP: item parameters (multilevel 2PL model)
   # IP[, 1]: item discrimination / slope
   # IP[, 2]: item intercept
   # Dimension: number of items * 2

Lz2 <- function(Dat, Group, Ability, IP){
  Group <- sort(Group)
  
  nquad <- 49L
  qnd <- seq(-6, 6, ,nquad)
  qwt <- dnorm(qnd)
  qwt <- qwt / sum(qwt)
  
  Pjk <- function(Dat, Ability, IP)
  { 
    n <- nrow(Dat)      # Number of individual test taker
    m <- ncol(Dat)      # Test Length
    pjk <- matrix(qwt,nrow = n ,ncol = length(qwt),byrow=TRUE)
    
    for (j in 1L:n){
      for (i in 1L:m)
      {
        # multiply trace lines
        if (Dat[j, i] == 1) {
          pjk[j, ] <- pjk[j, ] * plogis(IP[i, 1] * (qnd + Ability[j, 1]) + IP[i, 2]) }
        else {
          pjk[j, ] <- pjk[j, ] * ( 1 - plogis(IP[i, 1] * (qnd + Ability[j, 1]) + IP[i, 2])) }
      }
    }
    
    return(rowSums(pjk)) 
  }
  
  # Equation(10) L_02
  
  L02 <- function(Dat, Group, Ability, IP){
    aggregate(list(l02 = log(Pjk(Dat, Ability, IP))), list(Group = Group),  sum)
  }
  
  # Brute Force Computation of Expectation 
  EV.L02 <- function(Dat, Group, Ability, IP){
    m <- ncol(Dat)
    rp <- expand.grid( rep(list(0L:1L), m) )
    xi <- aggregate(list(xi = Ability[, 1]), list(Group = Group), mean)
    gsize <- data.frame(table(Group))
    
    EVl02 <- data.frame(Group = unique(Group),
                        El02 = NA,
                        Vl02 = NA)
    
    for (j in 1:nrow(xi)){
      pjk <- Pjk(Dat = rp, Ability = matrix(rep(xi[j, 2], nrow(rp))), IP = IP)
      plogp <- sum(pjk*(log(pjk)))
      p.logp2 <- sum(pjk*((log(pjk))^2))
      plogp.square <- plogp^2
      EVl02[j, 2] <- gsize[j, 2] * plogp
      EVl02[j, 3] <- gsize[j, 2] * (p.logp2 - plogp.square)
    }
    return(EVl02)
  }
  
  l02 <- L02(Dat, Group, Ability, IP)
  evl02 <- EV.L02(Dat, Group, Ability, IP)
  
  Lz2 <- (l02$l02 - evl02$El02)/sqrt(evl02$Vl02) 
  return(Lz2)
}

## -----------------------------------
## Output information
## -----------------------------------

## The function returns a dataframe:
    # Column 1: Group index
    # Column 2: l02  (function 10)
    # Column 3: El02 (function 11)
    # Column 4: Vl02 (function 12)
    # Column 5: Lz2  (function 13)

