popDynamics <- function(alpha,beta,Nadults,model="Beverton-Holt")
{
  if(model=="Beverton-Holt")
  {
    # population model is R ~ (CR * S)/(1+(CR-1)/K * S): reparameterized Beverton-Holt
    recruits <- (alpha*Nadults)/(1+((alpha-1)/beta)*Nadults)
  }
  if(model=="Ricker")
  {
    # population model is R ~ CR * S * e(-log(CR)/K * S): reparameterized Ricker
    recruits <- alpha*Nadults*exp(beta*Nadults)
  }
  return(list("recruits"=recruits))
}