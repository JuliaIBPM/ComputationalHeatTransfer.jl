# module Systems

export PHPSystem,Tube,Evaporator,Condenser,Liquid,Vapor,Wall,Mapping
# ,PHPResult

# using ..Tools

"""
PHPSystem is a struct containing
    γ
    Hc
    He
    θc
    θe
    ω
    ζ
    L dimensionless pipe total length
    Xc dimensionless condenser range
    Xe dimensionless evaporater range
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
end
"""
PHPSystem is a struct containing
    γ
    Hc
    He
    θc
    θe
    ω
    ζ
    L dimensionless pipe total length
    Xc dimensionless condenser range
    Xe dimensionless evaporater range
"""

mutable struct Liquid
    Hₗ::Float64
    ρ::Float64
    Cp::Float64
    α::Float64
    μₗ::Float64
    Xp::Array{Tuple{Float64,Float64},1}
    dXdt::Array{Tuple{Float64,Float64},1}
    Xarrays::Array{Array{Float64,1},1}
    θarrays::Array{Array{Float64,1},1}

end

"""
PHPSystem is a struct containing
    γ
    Hc
    He
    θc
    θe
    ω
    ζ
    L dimensionless pipe total length
    Xc dimensionless condenser range
    Xe dimensionless evaporater range
"""

mutable struct Vapor
    # γ::Float64
    Hᵥ::Float64
    k::Float64
    δmin::Float64
    P::Array{Float64,1}
    δ::Array{Float64,1}
end

"""
PHPSystem is a struct containing
    γ
    Hc
    He
    θc
    θe
    ω
    ζ
    L dimensionless pipe total length
    Xc dimensionless condenser range
    Xe dimensionless evaporater range
"""

mutable struct Wall
    # α::Float64
    ΔTthres::Float64
    Xstations::Array{Float64,1}
    Xarray::Array{Float64,1}
    θarray::Array{Float64,1}
end

"""
PHPSystem is a struct containing
    γ
    Hc
    He
    θc
    θe
    ω
    ζ
    L dimensionless pipe total length
    Xc dimensionless condenser range
    Xe dimensionless evaporater range
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
    γ
    Hc
    He
    θc
    θe
    ω
    ζ
    L dimensionless pipe total length
    Xc dimensionless condenser range
    Xe dimensionless evaporater range
"""
# add type names!
mutable struct PHPSystem
    tube::Tube
    # evaporator::Evaporator
    # condenser::Condenser
    liquid::Liquid
    vapor::Vapor
    wall::Wall
    mapping::Mapping
end


mutable struct PHPSystem_nomapping
    tube::Tube
    # evaporator::Evaporator
    # condenser::Condenser
    liquid::Liquid
    vapor::Vapor
    wall::Wall
end
# """
# PHPSystem is a struct containing
#     γ
#     Hc
#     He
#     θc
#     θe
#     ω
#     ζ
#     L dimensionless pipe total length
#     Xc dimensionless condenser range
#     Xe dimensionless evaporater range
# """
#
# struct PHPResult
#     t
#     Xp
#     dXdt
#     P
#     θ
#     M
# end
#

# end
