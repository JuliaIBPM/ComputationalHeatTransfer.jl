#### Operators for a system with a body and stationary points ####

export bc_constraint_rhs!

function bc_constraint_rhs!(Ts::ScalarData{N},sys::HeatConduction{NX,NY,N},t::Real) where {NX,NY,N}
    surface_temperature!(sys.ΔTs,sys,t)
    Ts .= sys.ΔTs
    return Ts
end


# Constraint operators, using stored regularization and interpolation operators
# B₁ᵀ = Eᵀ
function heatconduction_op_constraint_force!(out::Nodes{Primal,NX,NY},σ::ScalarData{N},sys::HeatConduction{NX,NY,N}) where {NX,NY,N}
    out .= sys.Rc*σ
    return out
end

#=
function _vel_ns_op_constraint_force!(u::Edges{Primal,NX,NY},τ,sys::HeatConduction{NX,NY,0}) where {NX,NY}
    return u
end

function _vel_ns_op_constraint_force!(u::Edges{Primal,NX,NY},τ::VectorData{N},sys::HeatConduction{NX,NY,N}) where {NX,NY,N}
    u .= sys.Rf*τ
    return u
end
=#

# B₂ = E
function heatconduction_constraint_op!(out::ScalarData{N},T::Nodes{Primal,NX,NY},sys::HeatConduction{NX,NY,N}) where {NX,NY,N}
    out .= sys.Ec*T
    return out
end
