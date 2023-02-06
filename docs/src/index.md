# ComputationalHeatTransfer.jl

Documentation for ComputationalHeatTransfer.jl

```julia
using Pkg
```


```julia
Pkg.activate(dirname(pwd()))
```

    [32m[1m  Activating[22m[39m project at `~/Documents/GitHub/ComputationalHeatTransfer.jl`



```julia
using ComputationalHeatTransfer
using LaTeXStrings
using JLD2
using Plots
gr()  

```

    ┌ Info: Precompiling ComputationalHeatTransfer [5fc296c8-2eb5-40dc-a46d-98a68011a900]
    └ @ Base loading.jl:1662





    Plots.GRBackend()



## Control Console


```julia
OHPtype = "ASETS-II OHP 2 LARGE HEATER"
Rn = 3e-6 # nucleation site radius
δfilm = 2e-5 # initial film thickness
ad_fac = 1.3 # film thickness factor
plate_d = 1.5e-3; # plate thickness
Eratio_plus = 0.15 + 0.5 # η+
Eratio_minus = 0.15 # η-
```




    0.15



# Properies

### Solid Physical parameters


```julia
ρₛ = 2730; # density
cₛ  = 8.93e02; # specific heat
kₛ  = 1.93e02; # heat conductivity
params = HeatConductionParameters(ρₛ ,cₛ ,kₛ ,thickness=plate_d)
```




    HeatConductionParameters(2730.0, 893.0, 193.0, 7.916682048820907e-5, 0.0015)



### Fluid Physical parameters


```julia
Tref = 291.2 # reference temperature
fluid_type = "Butane"
p_fluid = SaturationFluidProperty(fluid_type,Tref);
```

# Plate Conduction

### Geometry parameters


```julia
Lx = 0.1524; # plate size x
Ly = 0.0648; # plate size y
xlim = (-Lx/2,Lx/2) # plate x limits
ylim = (-Ly/2,Ly/2) # plate y limits
```




    (-0.0324, 0.0324)



### Set mesh size and maximum time step


```julia
Δx,Δt_max = setstepsizes(params.α,gridPe=8.0,fourier=0.3) 
```




    (0.0006333345639056725, 0.001520002953373614)



### Set up the evaporators and condensers


```julia
power = 40 #watts
hc = 3000.0
Tc = Tref;
eparams,cparams = OHPConfiguration(OHPtype,power,Tc,hc,Δx,hc2ratio=1/30);
```

### Set up OHP channels


```julia
x, y, xf, yf = construct_ohp_curve("ASETS",Δx)
ohp = BasicBody(x,y)
ohpgeom = ComputationalHeatTransfer.LineSourceParams(ohp);
```

### Create HeatConduction system


```julia
sys = HeatConduction(params,Δx,xlim,ylim,Δt_max,qline=ohpgeom,qflux=eparams,qmodel=cparams)
```




    Unbounded Heat conduction system on a grid of size 256 x 120 and 0 static immersed points




# OHP inner part


```julia
sys_tube,u,cbst = initialize_ohpsys(OHPtype,fluid_type,sys,p_fluid,Tref,δfilm,Eratio_plus,Eratio_minus,Rn);
```

### set time step


```julia
tspan = (0.0, 1.0);
dt_record = 0.01
num_data = tspan[2] / dt_record

tstep = 5e-4
```




    0.0005



### combine inner tube and plate together


```julia
u_plate = newstate(sys) .+ Tref # initialize T field
integrator_plate = init(u_plate,tspan,sys)
```




    t: 0.0
    u: (Primal nodes in a (nx = 256, ny = 120) cell grid of type Float64 data
      Number of Primal nodes: (nx = 255, ny = 119), Float64[])




```julia
prob = ODEProblem(ODE_innertube, u, tspan, sys_tube)
integrator_tube = init(prob, RK4(),save_on=false, dt=tstep, callback=cbst);
```

## Resume


```julia
# integrator_plate = integrator_plate_temp;
```


```julia
# integrator_tube = tube_hist[end];
```

## Start


```julia
boil_hist=[]
plate_T_hist = []
tube_hist_u  = []
tube_hist_t = []
tube_hist_θwall = []
```




    Any[]




```julia
using Distributed
using ProgressMeter
```


```julia
@showprogress for t in tspan[1]:tstep:tspan[2]

    timemarching!(integrator_tube,integrator_plate,tstep)
    
    if (mod(integrator_plate.t,dt_record) < 1e-6) || (mod(-integrator_plate.t,dt_record) < 1e-6)
        push!(plate_T_hist,deepcopy(temperature(integrator_plate))); 
        push!(tube_hist_θwall,deepcopy(integrator_tube.p.wall.θarray))
        push!(tube_hist_u,deepcopy(integrator_tube.u));
        push!(tube_hist_t,deepcopy(integrator_tube.t));
        integrator_plate_temp = deepcopy(integrator_plate)
#         println(sys.qline[1].q[1:5])
    end
    
end
# integrator_plate.t
```
