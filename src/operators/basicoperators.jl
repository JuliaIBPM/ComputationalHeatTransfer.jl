### Basic operators for any heat conduction system ###
# Note that the input vector `T` is temperature, and
# we expect the rhs of the equations to have units of dT/dt

export heatconduction_rhs!


# RHS of heat conduction equations
function heatconduction_rhs!(dT::Nodes{Primal,NX,NY},T::Nodes{Primal,NX,NY},sys::HeatConduction{NX,NY},t::Real) where {NX,NY}
  dT .= 0.0
  #_heatconduction_rhs_convectivederivative!(sys.Sc,T,sys,t)
  #_heatconduction_rhs_double_layer!(dT,sys,t)
  _heatconduction_rhs_double_layer_adiabatic!(dT,T,sys,sys.bctype)
  _heatconduction_rhs_forcing!(dT,sys,t)
  _heatconduction_rhs_qflux!(dT,sys)
  _heatconduction_rhs_qmodel!(dT,T,sys)
  _heatconduction_rhs_qline!(dT,sys)
  return dT
end

_heatconduction_rhs_qflux!(dT::Nodes{Primal,NX,NY},sys::HeatConduction{NX,NY}) where {NX,NY} =
      _heatconduction_rhs_qflux!(dT,sys,sys.qflux)


function _heatconduction_rhs_qflux!(dT::Nodes{Primal,NX,NY},sys::HeatConduction{NX,NY},qflux) where {NX,NY}
  @unpack params = sys
  @unpack ρ, c, d = params
  fact = 1.0/(ρ*c)
  fact = isnothing(d) ? fact : fact/d
  for q in qflux
    dT .+= fact.*q()
  end
  dT
end

_heatconduction_rhs_qflux!(dT::Nodes{Primal,NX,NY},sys::HeatConduction{NX,NY},::Nothing) where {NX,NY} = dT

_heatconduction_rhs_qmodel!(dT::Nodes{Primal,NX,NY},T::Nodes{Primal,NX,NY},sys::HeatConduction{NX,NY}) where {NX,NY} =
    _heatconduction_rhs_qmodel!(dT,T,sys,sys.qhdT)

function _heatconduction_rhs_qmodel!(dT::Nodes{Primal,NX,NY},T::Nodes{Primal,NX,NY},sys::HeatConduction{NX,NY},qhdT) where {NX,NY}
  @unpack params = sys
  @unpack ρ, c, d = params
  fact = 1.0/(ρ*c)
  fact = isnothing(d) ? fact : fact/d
  for q in qhdT
    dT .+= fact.*q(T)
  end
  dT
end

_heatconduction_rhs_qmodel!(dT::Nodes{Primal,NX,NY},T::Nodes{Primal,NX,NY},sys::HeatConduction{NX,NY},::Nothing) where {NX,NY} = dT


_heatconduction_rhs_qline!(dT::Nodes{Primal,NX,NY},sys::HeatConduction{NX,NY}) where {NX,NY} =
    _heatconduction_rhs_qline!(dT,sys,sys.qline)

function _heatconduction_rhs_qline!(dT::Nodes{Primal,NX,NY},sys::HeatConduction{NX,NY},qline) where {NX,NY}
  @unpack qline, params = sys
  @unpack ρ, c, d = params
  fact = 1.0/(ρ*c)
  fact = isnothing(d) ? fact : fact/d
  for ql in qline
    ql.cache1 .= ql.R*ql.q
    dT .-= fact.*ql.cache1
  end
  dT
end

_heatconduction_rhs_qline!(dT::Nodes{Primal,NX,NY},sys::HeatConduction{NX,NY},::Nothing) where {NX,NY} = dT



_heatconduction_rhs_forcing!(dT::Nodes{Primal,NX,NY},sys::HeatConduction{NX,NY},t) where {NX,NY} = _heatconduction_rhs_forcing!(dT,sys.qforce,cellsize(sys),t)

_heatconduction_rhs_forcing!(dT,::Nothing,Δx,t) = dT

function _heatconduction_rhs_forcing!(dT,qforce::Vector{<:ModulatedField},Δx,t)
  for q in qforce
    dT .+= q(t) #*Δx
  end
  dT
end


@inline _heatconduction_rhs_double_layer!(dT::Nodes{Primal,NX,NY},sys::HeatConduction{NX,NY,N,MT,
                                  ExternalInternalProblem},t::Real) where {NX,NY,N,MT} = dT

@inline _heatconduction_rhs_double_layer!(dT::Nodes{Primal,NX,NY},sys::HeatConduction{NX,NY,0},t::Real) where {NX,NY} = dT

function _heatconduction_rhs_double_layer!(dT::Nodes{Primal,NX,NY},sys::HeatConduction{NX,NY,N,MT,SD},t::Real) where {NX,NY,N,MT,SD}
    @unpack params = sys
    @unpack α = params
    Δx⁻¹ = 1/cellsize(sys)
    fact = α*Δx⁻¹
    surface_temperature_jump!(sys.ΔTs,sys,t)
    fill!(sys.Sc,0.0)
    sys.dlc(sys.Sc,sys.ΔTs)
    sys.Sc .*= fact
    dT .-= sys.Sc
end

@inline _heatconduction_rhs_double_layer_adiabatic!(dT::Nodes{Primal,NX,NY},T::Nodes{Primal,NX,NY},sys::HeatConduction{NX,NY,0},bctype) where {NX,NY} = dT

@inline _heatconduction_rhs_double_layer_adiabatic!(dT::Nodes{Primal,NX,NY},T::Nodes{Primal,NX,NY},sys::HeatConduction{NX,NY,0},::Type{AdiabaticBC}) where {NX,NY} = dT


@inline _heatconduction_rhs_double_layer_adiabatic!(dT::Nodes{Primal,NX,NY},T::Nodes{Primal,NX,NY},sys::HeatConduction{NX,NY,N},bctype) where {NX,NY,N} = dT


function _heatconduction_rhs_double_layer_adiabatic!(dT::Nodes{Primal,NX,NY},T::Nodes{Primal,NX,NY},sys::HeatConduction{NX,NY,N},::Type{AdiabaticBC}) where {NX,NY,N}
    @unpack params, ΔTs, Sc, Cc = sys
    @unpack k = params
    fill!(ΔTs,0.0)
    heatflux_constraint_op!(ΔTs,T,sys)
    ΔTs .= -(Cc\ΔTs)
    ΔTs ./= k
    fill!(Sc,0.0)
    heatflux_constraint_force!(Sc,ΔTs,sys)
    dT .-= Sc
end



#=
function _heatconduction_rhs_convectivederivative!(u::Edges{Primal,NX,NY},w::Nodes{Dual,NX,NY},sys::HeatConduction{NX,NY},t) where {NX,NY}
    Δx⁻¹ = 1/cellsize(sys)
    fill!(sys.Vv,0.0)
    velocity!(sys.Vv,w,sys,t)
    _unscaled_convective_derivative!(sys.Vv,sys)
    sys.Vv .*= Δx⁻¹
    u .-= sys.Vv
end

_unscaled_convective_derivative!(u::Edges{Primal,NX,NY},sys::HeatConduction{NX,NY}) where {NX,NY} =
      _unscaled_convective_derivative!(u,sys.Vtf,sys.DVf,sys.VDVf)

# Operates in-place on `u`, which comes in with the velocity field and
# returns the unscaled convective derivative
function _unscaled_convective_derivative!(u::Edges{Primal,NX,NY},
                                          Vtf::EdgeGradient{Primal,Dual,NX,NY},
                                          DVf::EdgeGradient{Primal,Dual,NX,NY},
                                          VDVf::EdgeGradient{Primal,Dual,NX,NY}) where {NX,NY}
    transpose!(Vtf,grid_interpolate!(DVf,u))
    DVf .= 0.0
    grad!(DVf,u)
    product!(VDVf,Vtf,DVf)
    u .= 0.0
    grid_interpolate!(u,VDVf)
    u
end
=#
