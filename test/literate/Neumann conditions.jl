using ComputationalHeatTransfer

using Plots
using LaTeXStrings

ρ = 1;
c = 1;
k = 0.1;
d = 0.1;
params = HeatConductionParameters(ρ,c,k,thickness=d)

xlim = (-2.1,2.1)
ylim = (-2.1,2.1)
Δx, Δt = setstepsizes(params.α,gridPe=0.4,fourier=2.0)

#=
Here we set up a rectangular boundary
=#
#bdry = Circle(2.0,1.5Δx)
bdry = Rectangle(2.0,2.0,1.5Δx,shifted=true)

#eb = Circle(1.0,1.5*Δx)
eb = Rectangle(0.5,1.0,1.5*Δx)
Te = RigidTransform((0.0,0.0),0.0)
Te(eb)
cb1 = Rectangle(0.5,1.0,1.5*Δx)
cb2 = Rectangle(0.5,1.0,1.5*Δx)
Tc1 = RigidTransform((1.0,0.0),0.0)
Tc2 = RigidTransform((-1.0,0.0),0.0)

Tc1(cb1)
Tc2(cb2)


qe = 5.0
hc = 10.0
Tc = -10.0

eparams = PrescribedHeatFluxRegion(qe,eb);
cparams1 = PrescribedHeatModelRegion(hc,Tc,cb1);
cparams2 = PrescribedHeatModelRegion(hc,Tc,cb2);

#=
Note the use of the flag `bctype=ComputationalHeatTransfer.AdiabaticBC` here.
We also specify the boundary of the domain in `bodies`.
=#
sys = HeatConduction(params,Δx,xlim,ylim,Δt,bodies=BodyList([bdry]),qflux=eparams,qmodel=[cparams1,cparams2],bctype=ComputationalHeatTransfer.AdiabaticBC)

u0 = newstate(sys)
tspan = (0.0,10.0)
integrator = init(u0,tspan,sys)

step!(integrator,2.0)

plot(temperature(integrator),sys.grid,legend=true,color=cgrad(:RdBu,rev=true))
plot!(bdry,fillalpha=0.0)

xc, yc = coordinates(temperature(integrator),sys.grid);

plot(xc,temperature(integrator)[:,floor(Int,size(sys.grid,2)/2)])
