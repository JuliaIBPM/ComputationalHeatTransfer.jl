# module Systems

export PHPSystem_nomapping,PHPSystem,Tube,Evaporator,Condenser,Liquid,Vapor,Wall,Mapping
# ,PHPResult

# using ..Tools

"""
Tube is a struct containing tube geometries
    d           ::Float64   tube characteristic diameter
    peri        ::Float64   tube perimeter (may be an arbitrary shape so cannot be derived from d)
    Ac          ::Float64   tube cross-sectional area
    L           ::Float64   tube one dimensional length
    L2D         ::Float64   tube length for each section of the tube
    angle       ::Float64   inclination angle of the OHP
    g           ::Float64   gravity
    closedornot ::Bool      if the tube is closed loop or not (open loop)
"""

mutable struct Tube
    d::Float64
    peri::Float64
    Ac::Float64
    L::Float64
    L2D::Float64
    angle::Float64
    g::Float64
    closedornot::Bool
    N::Int64  
    fluid_type::String
end
"""
Liquid is a struct containing liquid properties at a ref temperature
    Hₗ::Float64                              heat transfer coefficient between the wall and the pure liquid slug
    ρ::Float64                              liquid density
    Cp::Float64                             liquid specific heat capacity
    α::Float64                              liquid heat diffusivity
    μₗ::Float64                              liquid dynamic viscosity
    σ::Float64                              liquid surface tension
    Xp::Array{Tuple{Float64,Float64},1}     interface locations for each liquid slug
    dXdt::Array{Tuple{Float64,Float64},1}   interface velocity for each liquid slug
    Xarrays::Array{Array{Float64,1},1}      finite difference location points within each liquid slug
    θarrays::Array{Array{Float64,1},1}      finite difference temperature points within each liquid slug
"""

mutable struct Liquid
    Hₗ::Float64
    ρ::Float64
    Cp::Float64
    α::Float64
    μₗ::Float64
    σ::Float64
    Xp::Array{Tuple{Float64,Float64},1}
    dXdt::Array{Tuple{Float64,Float64},1}
    Xarrays::Array{Array{Float64,1},1}
    θarrays::Array{Array{Float64,1},1}

end

"""
Vapor is a struct containing vapor properties at a ref temperature
    Hᵥ::Float64             heat transfer coefficient between the wall and the pure vapor bubble
    k::Float64              heat conductivity of liquid in the film
    δmin::Float64           the delta with maximum heat transfer coefficient in the H interpolation
    P::Array{Float64,1}     pressure in each vapor
    δ::Array{Float64,1}     film thickness in each vapor
end
"""

mutable struct Vapor
    Hᵥ::Float64
    k::Float64
    δmin::Float64
    Eratio::Float64
    P::Array{Float64,1}
    # δ::Array{Float64,1}
    δfilm_deposit::Float64
    δstart::Array{Float64,1}
    δend::Array{Float64,1}
    Lfilm_start::Array{Float64,1}
    Lfilm_end::Array{Float64,1}
end

"""
Wall is a struct containing wall properties
    ΔTthres::Float64                superheat threshold to trigger boiling
    Xstations::Array{Float64,1}     locations of boiling stations on the wall
    Xarray::Array{Float64,1}        wall discrete location points for immersed boundary method
    θarray::Array{Float64,1}        wall discrete temperature points for immersed boundary method
end
"""

mutable struct Wall
    boil_type::String
    boil_interval::Float64
    Rn::Float64
    Xstations::Array{Float64,1}
    boiltime_stations::Array{Float64,1}
    Xarray::Array{Float64,1}
    θarray::Array{Float64,1}
end

"""
Mapping is a struct containing interpolation data
    θ_interp_walltoliquid   temperature interpolation from wall to OHP
    θ_interp_liquidtowall   temperature interpolation from OHP to wall
    H_interp_liquidtowall   heat transfer coefficient interpolation from OHP to wall
    P_interp_liquidtowall   pressure interpolation from OHP to wall
"""

# mutable struct Mapping
#     walltoliquid::Array{Tuple{Int64,Int64},1}
#     liquidtowall::Array{Array{Int64,1},1}
# end

mutable struct Mapping
    θ_interp_walltoliquid
    θ_interp_liquidtowall
    H_interp_liquidtowall
    P_interp_liquidtowall
end

"""
PHPSystem is a struct containing
    tube    ::Tube
    liquid  ::Liquid
    vapor   ::Vapor
    wall    ::Wall
    mapping ::Mapping
"""

mutable struct PHPSystem
    tube    ::Tube
    liquid  ::Liquid
    vapor   ::Vapor
    wall    ::Wall
    mapping ::Mapping
end

"""
PHPSystem_nomapping is a struct containing
    Tube
    Liquid
    Vapor
    Wall
It is used to construct PHPSystem, no other use.
"""
mutable struct PHPSystem_nomapping
    tube    ::Tube
    liquid  ::Liquid
    vapor   ::Vapor
    wall    ::Wall
end
