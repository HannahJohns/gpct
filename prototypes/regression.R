library(tidyverse)

# Prototype for regression using CLES

# Given a pairwise comparison matrix of all-to-all,
# assume data comes from a normal distribution modelled by
# Y ~ B1 X1 + B2 X2 + ... + Bn Xn + 0 + epsilon

# We dont' need an intercept term because everything is relative
# Note that we're assuming homoskedasticity here as well
df <- data.frame(gp = c(rep(0,300),rep(1,300)))
b <- 2
df$val <- sapply(df$gp,function(x){rnorm(1,b*x,1)})
ggplot(df,aes(x=val,fill=factor(gp)))+geom_density(alpha=0.7)


# Now that we have data, get a pairwise comparison matrix for all-to-all data.
# Does row i beat column j?
outer(df$val,df$val,">=") -> pairwise

# For each of these comparisons, assuming the model above we can get the probability of seeing the result.



# THIS IS WRONG, THESE PROBABILITIES ARE NOT INDEPENDENT. THE PROBABILITY MODEL HAS TO BE DEPENDENT ON THE RANKS OF EVERYTHING
# MEANING IT'S A NESTED INTEGRAL FOR DAYS



out <- {}
for(b in seq(0,3,length.out = 20))
{
for(e in seq(0.1,3,length.out = 20))
{

  cat(sprintf("b %f e %f\n", b,e))
    pars <- c(b,e)

    L <- 0
    for(i in 2:(nrow(pairwise)-1))
    {
    for(j in (i+1):nrow(pairwise))
    {
      result <- pairwise[i,j]

      prob <- ifelse(result,
                     pnorm(0,pars[1],sd=sqrt(pars[2]^2 + pars[2]^2),lower.tail = FALSE),
                     pnorm(0,pars[1],sd=sqrt(pars[2]^2 + pars[2]^2),lower.tail = TRUE)
      )

      L <- L + log(prob)

    }
    }

    out <- rbind(out,c(b,e,L))
}
}

out <- as.data.frame(out)

colnames(out) <- c("b","e","L")


# Note that we have drastically reduced epsilon values,
# This doesn't work
ggplot(out,aes(x=b,y=e,fill=L))+geom_raster()






# mean(pairwise)
#
# guess <- 2
#
# pnorm(0,guess,sd=sqrt(2),lower.tail = FALSE)
#
# table(c(pairwise),
# ifelse(c(pairwise),pnorm(0,guess,sd=sqrt(2),lower.tail = FALSE),1-pnorm(0,guess,sd=sqrt(2),lower.tail = FALSE))
# )


guesses <- seq(0,3,length.out = 60)
sapply(guesses,function(bhat){
sum(log(ifelse(c(pairwise),pnorm(0,bhat,sd=sqrt(2),lower.tail = FALSE),
       1-pnorm(0,bhat,sd=sqrt(2),lower.tail = FALSE))))
}) -> L

guesses[which(L==max(L))]






