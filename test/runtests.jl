using ComputationalHeatTransfer
using LaTeXStrings
using Revise
using LinearAlgebra
using DifferentialEquations
using Interpolations
using JLD2

divider = Sys.iswindows() ? "\\" : "/"
cd(dirname(pwd()))
cd("src")
includet(pwd()*divider*"OneDOHP.jl")
using .OneDOHP

using Test


@testset "liquid heat transfer" begin

tube_d = 1e-3; # tube diameter

ρ = 2700; # density
c = 8.97e02; # specific heat
k = 1.67e02; # heat conductivity
# d = 2e-3;
plate_d = 1e-3;  # plate thickness
params = HeatConductionParameters(ρ,c,k,thickness=plate_d)

fluid_type = "butane"
Tᵥ = 295.0

Cpₗ = CoolProp.PropsSI("CPMASS","T",Tᵥ,"Q",0.0,fluid_type)
ρₗ  = CoolProp.PropsSI("D","T",Tᵥ,"Q",0.0,fluid_type)
μₗ  = CoolProp.PropsSI("V","T",Tᵥ,"Q",0.0,fluid_type)
hₗ = CoolProp.PropsSI("H","T",Tᵥ,"Q",0.0,fluid_type)
kₗ = CoolProp.PropsSI("CONDUCTIVITY","T",Tᵥ,"Q",0.0,fluid_type)
Prₗ = CoolProp.PropsSI("PRANDTL","T",Tᵥ,"Q",0.0,fluid_type)

Cpᵥ = CoolProp.PropsSI("CPMASS","T",Tᵥ,"Q",1.0,fluid_type)
ρᵥ  = CoolProp.PropsSI("D","T",Tᵥ,"Q",1.0,fluid_type)
μᵥ  = CoolProp.PropsSI("V","T",Tᵥ,"Q",1.0,fluid_type);
hᵥ = CoolProp.PropsSI("H","T",Tᵥ,"Q",1.0,fluid_type)
kᵥ = CoolProp.PropsSI("CONDUCTIVITY","T",Tᵥ,"Q",1.0,fluid_type)
Prᵥ = CoolProp.PropsSI("PRANDTL","T",Tᵥ,"Q",1.0,fluid_type)

σ = CoolProp.PropsSI("I","T",Tᵥ,"Q",0.0,fluid_type)
P = CoolProp.PropsSI("P","T",Tᵥ,"Q",0.0,fluid_type)
R = CoolProp.PropsSI("GAS_CONSTANT","T",Tᵥ,"Q",1.0,fluid_type)
M = CoolProp.PropsSI("M","T",Tᵥ,"Q",1.0,fluid_type)
Rkg = R/M

νₗ = μₗ/ρₗ
νᵥ = μᵥ/ρᵥ;
hₗᵥ = hᵥ-hₗ;

Lx = 0.1524; # plate size x
Ly = 0.0648; # plate size y

xlim = (-Lx/2,Lx/2) .*1.0
ylim = (-Ly/2,Ly/2) .*1.0
Δx, Δt = setstepsizes(params.α,gridPe=8.0,fourier=0.3)


inch = 2.54e-2;

power = 40 #watts
total_heater_area = 2.0inch*2.0inch;

qe = power/total_heater_area
hc = 1500.0
Tc = Tᵥ

eb1 = Rectangle(0.5inch,1.0inch,1.5*Δx)
Tfe = RigidTransform((0.7inch,-0.0),0.0)
Tfe(eb1)

eb2 = Rectangle(0.5inch,1.0inch,1.5*Δx)
Tfe = RigidTransform((-0.7inch,-0.0),0.0)
Tfe(eb2)

cb1 = Rectangle(0.5inch*0.9,Ly*0.9/2,1.5*Δx)
Tfc = RigidTransform((-2.5inch,-0.0),0.0)
Tfc(cb1)

cb2 = Rectangle(0.5inch*0.9,Ly*0.9/2,1.5*Δx)
Tfc = RigidTransform((2.5inch,-0.0),0.0)
Tfc(cb2)

eparams1 = PrescribedHeatFluxRegion(qe,eb1);
eparams2 = PrescribedHeatFluxRegion(qe,eb2);
cparams1 = PrescribedHeatModelRegion(hc,Tc,cb1);
cparams2 = PrescribedHeatModelRegion(hc,Tc,cb2);

ds = 1.5Δx
nturn = 16
width_ohp = 46.25*1e-3
length_ohp = 133.83*1e-3
gap = plate_d
pitch = width_ohp/(2*nturn+1)
x0, y0 = length_ohp/2, width_ohp/2

# x, y, xf, yf = ComputationalHeatTransfer.construct_ohp_curve(nturn,pitch,length_ohp,gap,ds,x0,y0,false,false,3pi/2)
one_loop_gap = 10*plate_d
x, y, xf, yf = construct_oneloop_curve(0,0,ds,length_ohp,one_loop_gap,pi/2)


ohp = BasicBody(x,y)

# plot(ohp,fillalpha=0,linecolor=:blue,xlims=xlim,ylims=ylim,framestyle = :box)

ohpgeom = ComputationalHeatTransfer.LineSourceParams(ohp)

sys = HeatConduction(params,Δx,xlim,ylim,Δt,qline=ohpgeom,qflux=[eparams1,eparams2],qmodel=[cparams1,cparams2])

t_to_nondi_t = μₗ/ρₗ/(tube_d^2)

hl_to_nondihl = 4*tube_d/(μₗ*Cpₗ)

hv_to_nondihv =  4*Rkg*Tᵥ^2*tube_d/(P*hₗᵥ*νₗ)

nondihv_tonondihl = hl_to_nondihl / hv_to_nondihv

nondiqv_to_qv = Tᵥ *(pi*tube_d) / hv_to_nondihv

nondiql_to_ql = Tᵥ * (pi*tube_d)/ hl_to_nondihl

nondi_Q_to_Q = nondiqv_to_qv
function di_T_to_nondi_T(di_T;T0=295.0)

    di_T/T0

end

    ω = sqrt(2.45E+03);
    ℘L = 2.05E+05;
    k = 0.106024683

    θinitial=1.0
    θc = 1.0; # useless

    Hwc = 0.0; #not useful later on. H between wall and condenser,

    htₗ = 4.62e2 # Nusselt number 4.36
    Htₗ_real = 4*htₗ*tube_d/(μₗ*Cpₗ)
    Htₗ_OHP  = Htₗ_real/nondihv_tonondihl

    Hδ = kₗ/(tube_d/2) * hv_to_nondihv


    nondi_diameter = 6.56e-3 # dimensionless diameter
    L = 2.0 + one_loop_gap*pi/length_ohp  # total length of the pipe when streched to a 1D pipe
#     L = 34.35  # total length of the pipe when streched to a 1D pipe
    L2D = 1.0 # the actual length of the bended pipe in the real world
    angle = 0*pi/2 # inclination angle
    ΔTthres = 0.3/295.0 *100; # use a small threshold
    closedornot = true

    nucleatenum = 256
    Xstations = sort(rand(nucleatenum).*L);
#     Xstations = LinRange(0.0,L-L/nucleatenum,nucleatenum)
#     boilingΔt = 0.05

tube = Tube(nondi_diameter,L,L2D,angle,ΔTthres,closedornot);

# useless
Xe = map(tuple, [1.0], [3.0])
We = [0.0]
evap = Evaporator(Xe,We);

# useless
Xc = map(tuple, [0.0,3.0], [1.0,4.0])
cond = Condenser(θc,Xc,Hwc);

# X0,realratio = randomXp(L,numofslugs=60,chargeratio=0.45)
X0 = [(L/2,L/2 + 1.0)]

dXdt0_l = zeros(length(X0))
dXdt0_r = zeros(length(X0))
dXdt0 = map(tuple,dXdt0_l,dXdt0_r)


# construct liquids
N=numpts(ohp)

ρ = 102.593344 # density ratio
Xarrays,θarrays = constructXarrays(X0,N,θinitial,L)
liquids=Liquid(Htₗ_OHP,ρ,ω,℘L,X0,dXdt0,Xarrays,θarrays);

#
# realratio

# P = [1.0,1.0,1.0,1.0,1.0]; # closed end

γ = 1.4;
# Hδ = 4.64E+02
P = 0*zeros(length(X0)) .+ 1.0;
δ = 0*zeros(length(X0)) .+ 2.0E-01;
# δ = 0*zeros(length(X0)) .+ 2.78E-02;
# δ = 0*zeros(length(X0)) .+ 2.78E-02 .* 10; # use a random film thickness
vapors=Vapor(γ,Hδ,P,δ);

di_T_to_nondi_T(Tᵥ)

Xwallarray,θwallarray = constructXarrays(sys.qline[1].arccoord,L,θinitial)

α = 1.154e-5 # nondimensional thermal diffusivity = (αₐ d^2) / (νₗ height^2)
# Xwallarray,θwallarray = constructXarrays(L,N,θinitial)
Xwallarray,θwallarray = constructXarrays(sys.qline[1].arccoord,L,θinitial)
θwallarray .= di_T_to_nondi_T(Tᵥ);

wall = Wall(α,Δt,Xstations,Xwallarray,θwallarray);

sys0_nomapping = PHPSystem_nomapping(tube,evap,cond,liquids,vapors,wall);
θ_interp_walltoliquid, θ_interp_liquidtowall, H_interp_liquidtowall, P_interp_liquidtowall = sys_interpolation(sys0_nomapping)
mapping = Mapping(θ_interp_walltoliquid, θ_interp_liquidtowall, H_interp_liquidtowall, P_interp_liquidtowall);

sys0 = PHPSystem(tube,evap,cond,liquids,vapors,wall,mapping);

Lvaporplug = XptoLvaporplug(X0,sys0.tube.L,sys0.tube.closedornot)
M = nondi_PtoD(P) .* Lvaporplug
# M = P.^(1/γ).* Lvaporplug


u=[XMδtovec(X0,dXdt0,M,δ); liquidθtovec(sys0.liquid.θarrays)];

cb_boiling =  DiscreteCallback(boiling_condition,boiling_affect!)
# cb_boiling =  PeriodicCallback(boiling_affect!,0.01*t_to_nondi_t)
cb_merging =  DiscreteCallback(merging_condition,merging_affect!)

cbst = CallbackSet(cb_boiling,cb_merging);

tspan = (0.0, 1.0);
dt_record = tspan[2] /100;

tstep=1e-3

N_iter = 1
# tstep_plate = tstep/N_iter
# dt_record = tstep

u0 = newstate(sys)
integrator_plate = init(u0,tspan,sys)
# Tplate = temperature(integrator_plate);
temperature(integrator_plate) .= Tᵥ  .+ 1.0;

p = sys0
u=[XMδtovec(X0,dXdt0,M,δ); liquidθtovec(sys0.liquid.θarrays)];
prob = ODEProblem(ODE_innertube, u, tspan, p)
integrator_tube = init(prob, RK4(),save_everystep=false, dt=tstep*t_to_nondi_t, callback=cbst);

plate_hist = []
tube_hist  = []

currentsys = integrator_tube.p
currentsys = getcurrentsys(integrator_tube.u,currentsys);


@time for t in tspan[1]:tstep:tspan[2]


#   for j in 1:N_iter
        currentsys.wall.θarray = di_T_to_nondi_T(temperature_linesource(integrator_plate))
        currentsys = getcurrentsys(integrator_tube.u,currentsys)


        nondi_qtmp = sys_to_heatflux(currentsys)
        qtmp = nondi_Q_to_Q*nondi_qtmp
        set_linesource_strength!(sys,qtmp)


#      ADI_timemarching!(temperature(integrator_plate),sys,tstep)
#     end
    integrator_plate.t += tstep

    # reinitialize the integrator_tube to avoid some mysterious problems
    prob = ODEProblem(ODE_innertube, deepcopy(integrator_tube.u), (integrator_plate.t*t_to_nondi_t-tstep*t_to_nondi_t,integrator_plate.t*t_to_nondi_t), currentsys)
    integrator_tube = init(prob, RK4(), callback=cbst, dt=tstep*t_to_nondi_t);
#     step!(integrator_tube);
    step!(integrator_tube,tstep*t_to_nondi_t,true);

    if (mod(integrator_plate.t,dt_record) < 1e-6) || (mod(-integrator_plate.t,dt_record) < 1e-6)
        push!(plate_hist,deepcopy(integrator_plate));
        push!(tube_hist,deepcopy(integrator_tube));
        println(integrator_plate.t)
    end

end

        push!(plate_hist,deepcopy(integrator_plate));
        push!(tube_hist,deepcopy(integrator_tube));
        println(integrator_plate.t)

sysfinal = []
for i = 1:length(tube_hist)
    push!(sysfinal, deepcopy(getcurrentsys(tube_hist[i].u,tube_hist[i].p)))
end

# plot(ohp,fillalpha=0,line_z=temperature_linesource(integrator_plate))

# maximum(temperature(plate_hist[end])[:])
#
# Tmax = maximum(temperature(plate_hist[end-1])[:])
# Tmin = 0.0
# @gif for i = 1:1:length(plate_hist)
# # @gif     for i = 1:1:1
#
# plot(temperature(plate_hist[i]),sys.grid,legend=true,color=cgrad(:RdBu,rev=true),clim=(Tmin,Tmax),line_z=0,xlabel="x [m]",ylabel="y [m]",title=string("time = ", round(plate_hist[i].t, digits=2), "[s] \n",  "T - T0 [K]"))
# end
#
# plot(sys.qline[1].arccoord,sys_to_heatflux(sysfinal[end]))

# @gif for i=1:1:1
#     Htmp = sys_to_Harray(sysfinal[i])
#     plot(ohp,fillalpha=0,linecolor=cgrad([:gold, :blue],rev=true),line_z=Htmp,xlabel="x ",ylabel="y ",border=:none,axis=nothing)
# end

# @gif for i=1:1:length(sysfinal)
#     Htmp = sys_to_Harray(sysfinal[i])
#     plot(ohp,fillalpha=0,linecolor=cgrad([:gold, :blue],rev=true),line_z=Htmp,xlabel="x ",ylabel="y ",border=:none,axis=nothing)
# end

# @gif for ii=1:length(sysfinal)
#      plot(sysfinal[ii],plottype="T",xlim=(0.0,2.0))
# end

# plate_hist, integrator_plate = load("F:\OHPresults\1mm-correctliquidH\plate_OHP_0.2thickness.jld2", "plate_hist",  "integrator_plate")

# tube_hist, integrator_tube = load("F:\OHPresults\1mm-correctliquidH\tube_OHP_0.2thickness.jld2", "tube_hist",  "integrator_tube")

# save("plate_OHP_mediumlowH.jld2", "plate_hist", plate_hist, "integrator_plate", integrator_plate)

# save("tube_OHP_mediumlowH.jld2", "tube_hist", tube_hist, "integrator_tube",integrator_tube)

# inch = 2.54e-2;

# # x = [-2.7inch,0.0,2.7inch];
# # y = [0.0,0.0,0.0];

# x = [0.0]
# y = [0.0]

# X =  VectorData(x,y);

# H = Regularize(X,cellsize(sys),I0=origin(sys.grid))
# g = ScalarData(X);

# ghist = []
# thist = []
# for i = 1:length(plate_hist)
#     H(g,temperature(plate_hist[i]))
#     append!(ghist,deepcopy(g))
#     append!(thist,plate_hist[i].t)
# end

# RTDx,RTDy = load("RTD4.jld2","RTDx","RTDy")

# plot(thist,ghist)
# scatter!(RTDx .- RTDx[1],RTDy .- RTDy[1] .+ 295.0)

paraₗ = 4*htₗ/(tube_d*ρₗ*Cpₗ)
t = LinRange(tspan[1],tspan[2],length(sysfinal))

Thist_refₗ = (Tᵥ + 1.0) .- exp.(-paraₗ .* (t .- t[1]));


Verification_Thist = zeros(length(sysfinal))
Verification_Thist[1] = Tᵥ
for i = 1:length(sysfinal)-1
    Verification_Thist[i+1] = (maximum(sysfinal[i].liquid.θarrays[1]) * Tᵥ)
end

@test norm(Verification_Thist - Thist_refₗ,Inf) < 1e-4


end;

# ##using TestSetExtensions
# using Literate
#
# const GROUP = get(ENV, "GROUP", "All")
#
# notebookdir = "../examples"
# docdir = "../docs/src/manual"
# litdir = "./literate"
#
# if GROUP == "All" || GROUP == "Auxiliary"
#   #include("pointforce.jl")
# end
#
#
# if GROUP == "All" || GROUP == "Notebooks"
#   for (root, dirs, files) in walkdir(litdir)
#     for file in files
#       #endswith(file,".jl") && startswith(file,"6") && Literate.notebook(joinpath(root, file),notebookdir)
#       endswith(file,".jl") && Literate.notebook(joinpath(root, file),notebookdir)
#     end
#   end
# end
#
# if GROUP == "Documentation"
#   for (root, dirs, files) in walkdir(litdir)
#     for file in files
#       endswith(file,".jl") && Literate.markdown(joinpath(root, file),docdir)
#     end
#   end
# end
