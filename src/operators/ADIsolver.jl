### ADI operators for any heat conduction system ###

export ADI_laplacian!,ADI_newT!,ADI_timemarching!

# RHS of heat conduction equations
function ADI_laplacian!(NX,NY)

    N = NX*NY
    dlx = ones(N-1)
    dx  = -2ones(N)
    drx = ones(N-1)

    dx[1:NX:end] .= -1.0
    dx[NX:NX:end] .= -1.0
    dlx[NX:NX:end] .= 0.0
    drx = dlx;

    dly = ones(N-1)
    dy  = -2ones(N)
    dry = deepcopy(dly);

    dy[1:NY:end] .= -1.0
    dy[NY:NY:end] .= -1.0
    dly[NY:NY:end] .= 0.0
    dry = dly;

    DDx = Tridiagonal(dlx,dx,drx)
    DDy = Tridiagonal(dly,dy,dry)

  return DDx,DDy
end

# RHS of heat conduction equations
function ADI_newT!(T::Nodes{Primal,NX,NY},sys::HeatConduction{NX,NY,N},Δt::Real) where {NX,NY,N}
        @unpack DDx,DDy,params = sys
        @unpack α = params

        Fo = α*Δt/cellsize(sys)/cellsize(sys)
        fact = Fo/2

        Told = transpose(T)[:]

        T2Dold = reshape((I + fact*DDy)*Told,NY-1,NX-1)
        Tvecold = transpose(T2Dold)[:]
        Ttemp = (I - fact*DDx)\Tvecold

        T2Dtemp= reshape((I + fact*DDx)*Ttemp,NX-1,NY-1)
        Tvectemp = transpose(T2Dtemp)[:];
        Tnew  = (I - fact*DDy)\Tvectemp

        T .= transpose(reshape(Tnew,NY-1,NX-1))

  return T
end

# RHS of heat conduction equations
function ADI_timemarching!(T::Nodes{Primal,NX,NY},sys::HeatConduction{NX,NY,N},Δt::Real) where {NX,NY,N}

        dT = Nodes(Primal,NX,NY)
        Trhs = Δt * heatconduction_rhs!(dT::Nodes{Primal,NX,NY},T::Nodes{Primal,NX,NY},sys::HeatConduction{NX,NY},0.0)

        T = ADI_newT!(T::Nodes{Primal,NX,NY},sys::HeatConduction{NX,NY,N},Δt::Real)

        T .= T + Trhs

  return T
end
