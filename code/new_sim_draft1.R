
#Parameters
b=0.3 #biting rate
Tvh = 0.5 #vector transmission-probability susceptible vector being infected after contact with infectious host
Nh = 10000 #number of hosts in system
Nv = 100000 #number of vectors in system
muv = 0.1 #vector mortality rate
gammah = 0.1 #host recovery rate
alphab = 1 #turn on noise for biting rate
alphat = 1 #turn on noise for transmission rate
alpham = 1 #turn on noise for mortality rate

Thv = 0 #host transmission-probability susceptible host being infected after contact with infectious vector
#this one is our level for R0-code later

sigma = 0 #environmental noise level
#this is varied from 0 to 0.3

#initial values
H0 = 0
V0 = 10

#set up vectors to store results with value at time=0
H <- c(H0)
V <- c(V0)
time <- c(0)

H_t <- H0
V_t <- V0
deltat <- 1


for (t in (1:3650)){
  #following equations in paper draft-double check these
  
  H_t_1 <- H_t + ((b*Thv/Nh)*V_t*(Nh-H_t)-gammah*H_t)*deltat+
    sqrt(((b*Thv/Nh)*V_t*(Nh-H_t)+gammah*H_t))*rnorm(1, mean=0, sd=sqrt(deltat))+
    sigma*alphab*(Thv/Nh)*V_t*(Nh-H_t)*rnorm(1, mean=0, sd=sqrt(deltat))
  
  V_t_1 <- V_t + ((b*Tvh/Nh)*H_t*(Nv-V_t)-muv*V_t)*deltat+
    sqrt(((b*Tvh/Nh)*H_t*(Nv-V_t)+muv*V_t))*rnorm(1, mean=0, sd=sqrt(deltat))+
    sigma*(alphab*(Tvh/Nh)*H_t*(Nv-V_t)+alphat*(b/Nh)*(Nv-V_t)+alpham*V_t)*
    rnorm(1, mean=0, sd=sqrt(deltat))
  
  H <- append(H, H_t_1)
  V <- append(V, V_t_1)
}