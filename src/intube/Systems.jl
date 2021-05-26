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
    L::Float64
    L2D::Float64
    angle::Float64
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

mutable struct Evaporator
    # He::Float64
    # θe::Float64
    Xe::Array{Tuple{Float64,Float64},1}
    We::Array{Float64,1}
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

mutable struct Condenser
    # Hc::Float64
    θc::Float64
    Xc::Array{Tuple{Float64,Float64},1}
    Hwc::Float64
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
    γ::Float64
    Hₗ::Float64
    ρ::Float64
    ω::Float64
    ℘L::Float64
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
    γ::Float64
    k::Float64
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
    α::Float64
    Δt::Float64
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

mutable struct Mapping
    walltoliquid::Array{Tuple{Int64,Int64},1}
    liquidtowall::Array{Array{Int64,1},1}
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
    evaporator::Evaporator
    condenser::Condenser
    liquid::Liquid
    vapor::Vapor
    wall::Wall
    mapping::Mapping
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
