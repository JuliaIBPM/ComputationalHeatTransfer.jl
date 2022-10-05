export SaturationFluidProperty

using CoolProp

struct SaturationFluidProperty
    Cpₗ::Float64
    ρₗ::Float64
    μₗ::Float64
    hₗ::Float64
    kₗ::Float64
    Prₗ::Float64

    Cpᵥ::Float64
    ρᵥ::Float64
    μᵥ::Float64
    hᵥ::Float64
    kᵥ::Float64
    Prᵥ::Float64

    σ::Float64
    P::Float64
    R::Float64
    M::Float64
    Rkg::Float64

    αₗ::Float64
    νₗ::Float64
    νᵥ::Float64
    hₗᵥ::Float64
end

function SaturationFluidProperty(Cpₗ,ρₗ,μₗ,hₗ,kₗ,Prₗ,Cpᵥ,ρᵥ,μᵥ,hᵥ,kᵥ,Prᵥ,σ,P,R,M)

    Rkg = R/M
    αₗ = kₗ/ρₗ/Cpₗ
    νₗ = μₗ/ρₗ
    νᵥ = μᵥ/ρᵥ;
    hₗᵥ = hᵥ-hₗ;

    SaturationFluidProperty(Cpₗ,ρₗ,μₗ,hₗ,kₗ,Prₗ,Cpᵥ,ρᵥ,μᵥ,hᵥ,kᵥ,Prᵥ,σ,P,R,M,Rkg,αₗ,νₗ,νᵥ,hₗᵥ)
end

function SaturationFluidProperty(fluid_type::String,Tᵥ)
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

    SaturationFluidProperty(Cpₗ,ρₗ,μₗ,hₗ,kₗ,Prₗ,Cpᵥ,ρᵥ,μᵥ,hᵥ,kᵥ,Prᵥ,σ,P,R,M)
end
