

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

#initial values
H0 = 0
V0 = 10

#R0 values to test
R0s <- c(1, 2)
  #c(0.75, 0.95, 1, 1.05, 1.25, 2, 4, 6.5)

#Sigma values (environmental noise level in 0-0.3) to test
sigmas <- c(0, 0.1)
  #c(0, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3)

#set up end dataframe
#track probability of endemic disease-lasting until end of simulation
#probability of outbreaks greater than 10 and 100 hosts
#mean max number of cases
#mean outbreak duration
results <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(results) <- c('R0', 'sigma', 'prob_e', 'prob_o_10', 'prob_o_100', 'max_cases', 'duration')


for (q in (1:length(R0s))) {
  R0=R0s[q]

  Thv <- (R0^2) / ((b^2 * Tvh * Nv)/(Nh * gammah * muv)) #calculate Thv based on R0-lever
  #host transmission-probability susceptible host being infected after contact with infectious vector

for (p in (1:length(sigmas))){
  sigma = sigmas[p] #environmental noise level


#set up vectors to store results with value at time=0
H <- c(H0)
V <- c(V0)
time <- c(0)

H_t <- H0
V_t <- V0
deltat <- 1

#set up a dataframe to store summary results of each R0, sigma combination
temp_results <- data.frame(matrix(ncol = 8, nrow = 0))
colnames(temp_results) <- c('trial', 'R0', 'sigma', 'prob_e', 'prob_o_10', 'prob_o_100', 'max_cases', 'duration')

for (w in (1:10)){

for (t in (1:3650)){
  #following equations in paper draft-double check these ****
  
  H_t_1 <- H_t + ((b*Thv/Nh)*V_t*(Nh-H_t)-gammah*H_t)*deltat+
    sqrt(((b*Thv/Nh)*V_t*(Nh-H_t)+gammah*H_t))*rnorm(1, mean=0, sd=sqrt(deltat))+
    sigma*alphab*(Thv/Nh)*V_t*(Nh-H_t)*rnorm(1, mean=0, sd=sqrt(deltat))
  
  V_t_1 <- V_t + ((b*Tvh/Nh)*H_t*(Nv-V_t)-muv*V_t)*deltat+
    sqrt(((b*Tvh/Nh)*H_t*(Nv-V_t)+muv*V_t))*rnorm(1, mean=0, sd=sqrt(deltat))+
    sigma*(alphab*(Tvh/Nh)*H_t*(Nv-V_t)+alphat*(b/Nh)*(Nv-V_t)+alpham*V_t)*
    rnorm(1, mean=0, sd=sqrt(deltat))
  
  #I am not sure if we need this
  #Ensure that H and V do not go negative
  if (H_t_1<0){
    H_t_1 = 0
  }
  if (V_t_1<0){
    V_t_1 = 0
  }
  
  #Add the values to vector
  H <- append(H, H_t_1)
  V <- append(V, V_t_1)
  time <- append(time, t)
  
  #advance one timestep
  H_t <- H_t_1
  V_t <- V_t_1
  
}
  #make dataframe of results from one simulation
  sim_results <- data.frame(time, V, H)
  endemic <- tail(sim_results$H, n=1) #if disease was endemic in simulation H>0 at end
  max_cases <- max(sim_results$H)
  
  
  duration <- max(sim_results[sim_results$H>1, 'time'])-min(sim_results[sim_results$H>1, 'time']) 
    #determine duration by the time H>1. 
  #***Here I am assuming there is one outbreak per simulation, not multiple
  #I do not think this assumption is valid, but I cannot think of another way####
  
  #set up dataframe to store this simulation's results
  sim_summary<-data.frame(w,R0, sigma, endemic, max_cases, max_cases, max_cases, duration)
  colnames(sim_summary) <- c('trial', 'R0', 'sigma', 'prob_e', 'prob_o_10', 'prob_o_100', 'max_cases', 'duration')
  
  sim_summary$prob_e[sim_summary$prob_e>0] <- 1 #assign value of 1 if there are infected hosts at the end
  sim_summary$prob_o_10[sim_summary$prob_o_10<11] <- 0 #assign value of 0 if there is never more than 10 infected hosts at any point
  sim_summary$prob_o_10[sim_summary$prob_o_10>10] <- 1 #assign value of 1 if there are more than 10 infected hosts at any point
  sim_summary$prob_o_100[sim_summary$prob_o_100<101] <- 0 #assign value of 0 if there is never more than 100 infected hosts at any point
  sim_summary$prob_o_100[sim_summary$prob_o_100>100] <- 1 #assign value of 1 if there are more than 100 infected hosts at any point
  
  #add the results from this run to the temp_results dataframe
  temp_results <- rbind(temp_results, sim_summary)
}
  #set up dataframe storing the results of this R0-sigma combination
  R0_sigma_summary <-data.frame(mean(temp_results$R0), mean(temp_results$R0), 
                                mean(temp_results$prob_e), mean(temp_results$prob_o_10), 
                                mean(temp_results$prob_o_100), mean(temp_results$max_cases), mean(temp_results$duration))
  colnames(R0_sigma_summary) <- c('R0', 'sigma', 'prob_e', 'prob_o_10', 'prob_o_100', 'max_cases', 'duration')
  
  #add the results from this r0 sigma combination to the results dataframe
  results <- rbind(results, R0_sigma_summary)
  
}
}
