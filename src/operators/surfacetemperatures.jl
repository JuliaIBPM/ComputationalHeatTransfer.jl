# Routines to compute surface temperatures and their jumps

@inline function surface_temperature_jump!(Δus::VectorData{N},sys::HeatConduction{NX,NY,N},t::Real) where {NX,NY,N}
    #assign_temperature!(sys.Sb,sys.bodies,sys.bodytemps,t)
    fill!(sys.Sb,0.0)
    surface_temperature_jump!(Δus,sys.Sb,sys)
    return Δus
  end

@inline surface_temperature_jump!(ΔTs::ScalarData{N},Ts::ScalarData{N},
                               ::HeatConduction{NX,NY,N,MT,InternalProblem}) where {NX,NY,N,MT} =
                               (ΔTs .= -Ts; ΔTs)

@inline surface_temperature_jump!(ΔTs::ScalarData{N},Ts::ScalarData{N},
                               ::HeatConduction{NX,NY,N,MT,ExternalProblem}) where {NX,NY,N,MT} =
                               (ΔTs .= Ts; ΔTs)

@inline surface_temperature_jump!(ΔTs::ScalarData{N},Ts::ScalarData{N},
                               ::HeatConduction{NX,NY,N,MT,ExternalInternalProblem}) where {NX,NY,N,MT} =
                               (ΔTs .= 0.0; ΔTs)


@inline function surface_temperature!(T̄s::ScalarData{N},sys::HeatConduction{NX,NY,N},t::Real) where {NX,NY,N}
    #assign_temperature!(sys.Sb,sys.bodies,sys.bodytemps,t)
    fill!(sys.Sb,0.0)
    surface_temperature!(T̄s,sys.Sb,sys)
    return T̄s
  end

@inline surface_temperature!(T̄s::ScalarData{N},Ts::ScalarData{N},
                               ::HeatConduction{NX,NY,N,MT,ExternalInternalProblem}) where {NX,NY,N,MT} =
                               (T̄s .= Ts; T̄s)

@inline surface_temperature!(T̄s::ScalarData{N},Ts::ScalarData{N},
                               ::HeatConduction{NX,NY,N,MT,SD}) where {NX,NY,N,MT,SD} =
                               (T̄s .= 0.5*Ts; T̄s)

@inline function relative_surface_temperature!(T̄sr::VectorData{N},sys::HeatConduction{NX,NY,N},t::Real) where {NX,NY,N}
  assign_temperature!(sys.Sb,sys.bodies,sys.bodytemps,t)
  surface_temperature!(T̄sr,sys.Sb,sys)
  T̄sr .-= sys.Sb
end
