require(deSolve)
require(tidyr)
require(ggplot2)

stoc.RM = function(t.vec, y, parms){
  
  with(as.list(c(parms, y)), {
    
    # b0.stoc = rpois(1, b0)
    b0.stoc = rnorm(1, sd=size, mean=b0)
    
    dSh = (Kh - Sh)*mh - b0.stoc/(Sh + Ih)*thv*Sh*Iv + gh*Ih
    dIh = b0.stoc/(Sh + Ih)*thv*Sh*Iv - (gh+mh)*Ih
    
    dSv = (Kv - Sv)*mv - b0.stoc/(Sh + Ih)*tvh*Ih*Sv
    dIv = b0.stoc/(Sh + Ih)*tvh*Ih*Sv - mv*Iv
    
    list(c(dSh, dIh, dSv, dIv))
  })
}


simulate.noise = function(parms.grid, y0.vec, t.vec, nsim){
  noise.out = NULL
  int = 1
  for(i in 1:dim(parms.grid)[1]){
    parms.use = parms.grid[i, ]
    
    for(j in 1:nsim){
      print(int)
      
      ode.temp = data.frame(ode(y=y0, times=t.vec, func=stoc.RM, parms=parms.use, method="rk4"))
      # ode.temp = data.frame(ode(y=y0, times=t.vec, func=stoc.RM, parms=parms.grid))
      
      ode.temp2 = cbind(ode.temp, size = rep(parms.use$size, dim(ode.temp)[1]), int = as.factor(rep(int, dim(ode.temp)[1])))
      noise.out = rbind(noise.out, ode.temp2)
      
      int = int + 1
    }
  }
  return(noise.out)
}

estimate_mode <- function(x,from=min(x), to=max(x)) {
  d <- density(x, from=from, to=to)
  d$x[which.max(d$y)]
}

stoc.stats = function(sim.df){
  mode.Ih.out = NULL
  mode.Iv.out = NULL
  mode.Sh.out = NULL
  mode.Sv.out = NULL
  size.out = NULL
  for(i in 1:max(as.numeric(sim.df$int))){
    sim.data.use = subset(sim.df, int %in% i)
    mode.Ih.out = append(mode.Ih.out, estimate_mode(sim.data.use$Ih))
    mode.Iv.out = append(mode.Iv.out, estimate_mode(sim.data.use$Iv))
    mode.Sh.out = append(mode.Sh.out, estimate_mode(sim.data.use$Sh))
    mode.Sv.out = append(mode.Sv.out, estimate_mode(sim.data.use$Sv))
    size.out = append(size.out, sim.data.use$size[1])
  }
  stats.out = data.frame(mode.Ih.out, mode.Iv.out, mode.Sh.out, mode.Sv.out, size.out)
  
  return(stats.out)
}


###################################
#Implementation
###################################
t0 = 0
tmax = 100
dt = 0.01
t.vec = seq(t0, tmax, dt)

#host params 
Kh = 10
mh = 1.0

thv = 1.0
gh = 3.0

#vector parms 
Kv = 100
mv = 2.0

tvh = 1.0
b0 = 1.0

size = seq(0, 100, 10)

R0 = sqrt((b0^2*thv*tvh*Kv)/(Kh*(mh+gh)*mv))

# parms.grid = c(Kh=Kh, mh=mh, thv=thv, gh=gh,
#               Kv=Kv, mv=mv, tvh=tvh, b0=b0,
#               size=size)

parms.grid = expand.grid(Kh=Kh, mh=mh, thv=thv, gh=gh,
                         Kv=Kv, mv=mv, tvh=tvh, b0=b0,
                         size=size)

#intial conditions
Sh0 = Kh
Ih0 = 0

Sv0 = Kv
Iv0 = 1

y0 = c(Sh=Sh0, Ih=Ih0, Sv=Sv0, Iv=Iv0)

# ode.out = data.frame(ode(y=y0, times=t.vec, func=stoc.RM, parms=parms.grid))
# ode.out = data.frame(ode(y=y0, times=t.vec, func=stoc.RM, parms=parms.grid, method="rk4"))
# 
# ode.out2 = gather(ode.out, key="state", value="number", Sh:Iv)
# 
# ggplot(ode.out2) + geom_line(aes(x=time, y=number, color=state))

sim.test = simulate.noise(parms.grid, y0, t.vec, 10)

sim.test2 = gather(sim.test, key="state", value="number", Sh:Iv)

ggplot(subset(sim.test2, state %in% "Iv")) + geom_line(aes(x=time, y=number, color=int)) + facet_wrap(~size) + theme(legend.position = "none")

sim.stats.test = stoc.stats(sim.test) 

sim.stats.test2 =  gather(sim.stats.test, key="stat", value="value", mode.Ih.out:mode.Sv.out)

ggplot(sim.stats.test2) + geom_point(aes(x=size.out, y=value)) + facet_wrap(~stat, scales="free") + theme(legend.position = "none")
