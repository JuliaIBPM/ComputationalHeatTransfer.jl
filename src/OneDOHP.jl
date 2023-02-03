export inch,gravity,timemarching!

const inch = 2.54e-2; 
const gravity = 9.8;


include("intube/Systems.jl")
include("intube/Thermomodel.jl")
include("intube/Tools.jl")
include("intube/Mapping.jl")
include("intube/TimeMarching.jl")
include("intube/PostProcessing.jl")
include("intube/CoolProp.jl")
include("intube/Plotrecipe.jl")
include("intube/DrawingOHP.jl")
include("intube/PreProcessing.jl")
# include("intube/FluidProperty.jl")
include("intube/HeaterCondenser.jl")
include("intube/callback/boiling.jl")
include("intube/callback/vapormerging.jl")
include("intube/callback/liquidmerging.jl")
include("intube/callback/fixdx.jl")

# weakly coupled alternate time marching
function timemarching!(integrator_tube,integrator_plate,tstep)
    
    step!(integrator_tube,tstep,true);

    currentsys = deepcopy(getcurrentsys(integrator_tube.u,integrator_tube.p))
    currentsys.wall.θarray = deepcopy(temperature_linesource(integrator_plate))
    integrator_tube.p = deepcopy(currentsys)
    qtmp = deepcopy(sys_to_heatflux(currentsys))
    sys_plate = integrator_plate.p
    set_linesource_strength!(sys_plate,qtmp)

    ADI_timemarching!(temperature(integrator_plate),sys_plate,tstep)
    integrator_plate.t += tstep
    
    integrator_tube,integrator_plate
end