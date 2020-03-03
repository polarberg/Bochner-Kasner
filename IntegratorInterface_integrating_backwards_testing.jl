#DEPOT_PATH
#DEPOT_PATH[1]="D:\\.julia"
#DEPOT_PATH[2]="D:\\Julia-1.2.0\\local\\share\\julia"
#DEPOT_PATH[3]="D:\\Julia-1.2.0\\share\\julia"
#DEPOT_PATH

#import Pkg; using Pkg
#Pkg.add("DifferentialEquations"); Pkg.add("Plots")
using DifferentialEquations, Plots 

function eq1(dv,v,y,p,t)
    #dv[1] = -y[1]
    dv[1] = v[1]/(3*t) - ( (1)/t^(4/3) + 1*t^(2/3) ) * y[1] 
end

tspan = (2.0,1e-16); #tspan=big.(tspan)
#u0 = big.(2.0)
#du0= big.(1.0)

u0=2.0
du0=1.0

ode = SecondOrderODEProblem(eq1, [du0], [u0], tspan)
sol=solve(ode,DPRKN12(),dense=true,abstol=1e-21,reltol=1e-21, dt=-1.0, maxiters=1e50)


plot(sol)

z=zeros(length(sol.t),1)
length(z)
for i in 1:1:length(z)
    z[i]=sol.u[i][2]
end

z 
plot(sol.t,z)


integrator = init(ode,DPRKN12(),dense=true,abstol=1e-15,reltol=1e-15, dt=-1.0, maxiters=1e50)
dt=-1.0
step!(integrator,dt)
integrator
integrator.sol
sol