# module Thermomodel

export dMdtdynamicsmodel,wallmodel,liquidmodel,dynamicsmodel,sys_to_heatflux,sys_to_Harray,integrator_to_heatflux,integrator_to_Harray


function dynamicsmodel(u::Array{Float64,1},p::PHPSystem)


    du = zero(deepcopy(u))

    (Xp,dXdt0,M,δ)=vectoXMδ(u)


    numofliquidslug = p.tube.closedornot ? Integer( (length(u))/6 ) : Integer( (length(u) - 2)/6 )
    sys = deepcopy(p)

    # γ = sys.vapor.γ
    # ω = sys.liquid.ω
    # ℘L = sys.liquid.℘L
    d = sys.tube.d
    peri = sys.tube.peri
    Ac = sys.tube.Ac
    angle = sys.tube.angle
    g = sys.tube.g

    δ = sys.vapor.δ
    δarea = Ac .* (1 .- ((d .- 2*δ ) ./ d) .^ 2);
    dXdt_factor =  Ac ./ (Ac .- δarea)

    # println(sys.liquid.Xp)
    # println(sys.liquid.dXdt)

    δdeposit = δfilm
    Adeposit = getAdeposit(sys,δdeposit)
    # println(Adeposit)

    μₗ = sys.liquid.μₗ

    Lvaporplug = XptoLvaporplug(Xp,sys.tube.L,sys.tube.closedornot)
    Lliquidslug = XptoLliquidslug(Xp,sys.tube.L)
    height = getheight(Xp,sys.tube.L2D,sys.tube.angle)
    Xpvapor = getXpvapor(Xp,sys.tube.L,sys.tube.closedornot)


    ρ = M ./ Lvaporplug ./ (Ac .* ((d .- 2δ) ./ d).^2 )
    P = DtoP.(ρ)
    θ = PtoT.(P)
    #
    # println((Ac .* ((d .- 2δ[23]) ./ d).^2 ))
    # println(ρ[23])
    # println(δ[23])
    # println(M)
    # println(ρ)
    # println(P)



    ρₗ = p.liquid.ρ
# get differential equation factors
    lhs = ρₗ*Ac .* Lliquidslug

    rhs_press = Ac ./ lhs

# modify f !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dXdt_to_stress = -8*μₗ/d
    rhs_dXdt = peri .* Lliquidslug .* dXdt_to_stress ./ lhs


    rhs_g = Ac*ρₗ*g*cos(angle) ./ lhs



    # P = nondi_DtoP.(M./Lvaporplug)
    # θ = nondi_PtoT.(P)
    # P = real.((M./Lvaporplug .+ 0im).^(γ))
    # θ = real.((P .+ 0im).^((γ-1)/γ))

    # Prescribed pressure fields

    if p.tube.closedornot == false

            numofvaporbubble = numofliquidslug + 1

        for i = 1:numofliquidslug
            du[2*i-1] = u[2*numofliquidslug+2*i-1]
            du[2*i] = du[2*i-1]
            du[2*numofliquidslug + 2*i-1] = rhs_dXdt[i]*u[2*numofliquidslug + 2*i-1] +rhs_g[i]*(height[i][1]-height[i][end]) + rhs_press[i] * (P[i]-P[i+1])
            du[2*numofliquidslug + 2*i] = du[2*numofliquidslug + 2*i-1]

        end


        dMdt = dMdtdynamicsmodel(Xpvapor,sys)
        dδdt = (dMdt) .* d^2 ./ ( (-4) .* (d .- 2δ) .* ρₗ .* Lvaporplug .* Ac)
        du[4*numofliquidslug+1:5*numofliquidslug+1] .= dMdt
        du[5*numofliquidslug+2:end] .= dδdt

        return du
        end

    if p.tube.closedornot == true

        numofvaporbubble = numofliquidslug

        v_liquid_left=zeros(numofliquidslug)
        v_liquid_right=zeros(numofliquidslug)
        v_vapor_left=zeros(numofvaporbubble)
        v_vapor_right=zeros(numofvaporbubble)

            for i = 1:numofliquidslug

    # temperary v, need a new momentum equation to solve for it
                v_momentum = u[2*numofliquidslug+2*i-1]

                v_liquid_left[i] = v_momentum + v_momentum*Adeposit[i][1]/(Ac-Adeposit[i][1])
                v_liquid_right[i] = v_momentum + v_momentum*Adeposit[i][end]/(Ac-Adeposit[i][end])

                # loop_index = i == numofliquidslug ? 1 : i+1
                # v_liquid_left[i] = v_momentum + v_momentum*δarea[i][1]/(Ac-δarea[i][1])
                # v_liquid_right[i] = v_momentum + v_momentum*δarea[loop_index]/(Ac-δarea[loop_index])

                # v_liquid_left[i] = u[2*numofliquidslug+2*i-1]
                # v_liquid_right[i] = u[2*numofliquidslug+2*i]

                rhs_dLdt = -v_momentum*(v_liquid_right[i]-v_liquid_left[i])/Lliquidslug[i]

                du[2*i-1] = v_liquid_left[i]
                du[2*i] = v_liquid_right[i]

                if i != numofliquidslug
                    du[2*numofliquidslug + 2*i-1] = rhs_dXdt[i]*u[2*numofliquidslug + 2*i-1] + rhs_g[i]*(height[i][1]-height[i][end]) + rhs_press[i] * (P[i]-P[i+1]) + rhs_dLdt
                else
                    du[2*numofliquidslug + 2*i-1] = rhs_dXdt[i]*u[2*numofliquidslug + 2*i-1] + rhs_g[i]*(height[i][1]-height[i][end]) + rhs_press[i] * (P[i]-P[1]) + rhs_dLdt
                end

                du[2*numofliquidslug + 2*i] = du[2*numofliquidslug + 2*i-1]

            end

            # anti_loop_index = (i != 1) ? i-1 : numofvaporbubble
            # v_momentum_leftslug = u[2*numofliquidslug+2*loop_index]
            # v_vapor_left[2:end] = v_momentum_leftslug + v_momentum_leftslug*Adeposit[loop_index][end]/(Ac-Adeposit[loop_index][end])
            v_vapor_left[2:end] = v_liquid_right[1:end-1]
            v_vapor_left[1] = v_liquid_right[end]
            v_vapor_right = v_liquid_left



            A_dδdt_right_vapor = [elem[1] for elem in Adeposit]

            A_dδdt_right_liquid = [elem[2] for elem in Adeposit]
            A_dδdt_left_vapor = zeros(size(A_dδdt_right_vapor))
            A_dδdt_left_vapor[2:end] = A_dδdt_right_liquid[1:end-1]
            A_dδdt_left_vapor[1] = A_dδdt_right_liquid[end]


            # dδdt = (dMdt - ρₗ.*v_vapor_left.*Adeposit[loop_index][end] + ρₗ.*v_vapor_right.*Adeposit[i][1]) .* d^2 ./ ( (-4) .* (d .- 2δ) .* ρₗ .* Lvaporplug .* Ac)
            dMdt_sensible,dMdt_latent = dMdtdynamicsmodel(Xpvapor,sys)
            # println(dMdt_latent[21:26])
            # println(dMdt_sensible[21:26])
            # println(ρₗ .* A_dδdt_right_vapor .* v_vapor_right)
            # println(ρₗ .* A_dδdt_left_vapor  .* v_vapor_left)


            E = ρₗ .* Ac .* 4 .* δ .* (d .- δ) ./ (d^2)
            C = ρₗ .* Ac .* 4 .* (d .- 2δ) ./ (d^2)
            # dδdt = (dMdt_latent) .* d^2 ./ ( (-4) .* (d .- 2δ) .* ρₗ .* Lvaporplug .* Ac)
            dδdt = (-dMdt_latent + ρₗ .* A_dδdt_right_vapor .* v_vapor_right - ρₗ .* A_dδdt_left_vapor  .* v_vapor_left - E .* (v_vapor_right-v_vapor_left)) ./ (C .* Lvaporplug)
            # println( v_vapor_right)
            # println( v_vapor_left)
            # println(dMdt_latent)
            du[4*numofliquidslug+1:5*numofliquidslug] .= dMdt_latent+dMdt_sensible
            du[5*numofliquidslug+1:end] .= dδdt # equals to 0 for now

            return du
    end


end

function dMdtdynamicsmodel(Xpvapor::Array{Tuple{Float64,Float64},1},sys::PHPSystem)

    # dMdt=zeros(length(Xpvapor))
    dMdt_sensible=zeros(length(Xpvapor))
    dMdt_latent=zeros(length(Xpvapor))

    # get Hvapor
    peri = sys.tube.peri
    Ac = sys.tube.Ac
    k = sys.vapor.k
    δ = sys.vapor.δ
    # Hvapor = k ./ δ
    Hvapor = Hfilm.(δ,[sys])

    #get θ
    P = sys.vapor.P
    # γ = sys.vapor.γ

    θarrays = sys.liquid.θarrays
    Xarrays = sys.liquid.Xarrays

    θ = PtoT.(P)
    # θ = nondi_PtoT.(P)
    # θ = real.((P .+ 0im).^((γ-1)/γ)) # isentropic

    Hfg = PtoHfg.(P)

    # modified here
    # Hfg = Hfg*1e-3
    # println(P)

    dx_wall = sys.wall.Xarray[2]-sys.wall.Xarray[1]

    for i = 1:length(Xpvapor)
        a, b = Xpvapor[i][1], Xpvapor[i][2];
        L_temp    = mod(b-a,sys.tube.L)		   ## note n=10
        n = Int64(div(L_temp,dx_wall) == 0 ? 1 : div(L_temp,dx_wall))


        dx_vapor = L_temp/n
        xs = mod.(a .+ (0:n) * dx_vapor,[sys.tube.L]);          ## n, right is 1:n * delta

        θ_wall_inter = sys.mapping.θ_interp_walltoliquid

        fx = map(θ_wall_inter, xs) .- θ[i]
        # println(maximum(fx))

        slope_r = getslope(θarrays[i][2],θarrays[i][1],Xarrays[i][2],Xarrays[i][1])
        slope_l = (i == 1) ? getslope(θarrays[1][end],θarrays[1][end-1],Xarrays[1][end],Xarrays[1][end-1]) : getslope(θarrays[i-1][end],θarrays[i-1][end-1],Xarrays[i-1][end],Xarrays[i-1][end-1])

        # slope_r = 21.95530408719068
        # slope_l = -21.95530408719068

        axial_rhs = Ac*k*(slope_r-slope_l) /Hfg[i]

        # simulate dryout
                δ_dry = 1e-5
                if δ[i] > δ_dry
                    dMdt_sensible[i] = 0.0
                    dMdt_latent[i] = (sum(fx) * dx_vapor) * (Hvapor[i] *peri/Hfg[i]) + axial_rhs
                else
                    dMdt_sensible[i] = (sum(fx) * dx_vapor) * (Hfilm(0.0,sys) *peri/Hfg[i])
                    dMdt_latent[i] = axial_rhs
                end
                # println(dMdt_latent)
                # println((sum(fx) * dx_vapor) * (Hvapor[i] *peri/Hfg[i]))
            end

            # for i = 1:length(Xpvapor)
            #     indexes = findall( x -> (mod(x[1],length(Xpvapor)) == mod(i,length(Xpvapor)) && (x[end] == -1)), sys.mapping.walltoliquid)
            #
            #         for j in length(indexes)
            #             dMdt[i] += Hvapor[i]*dx*(sys.wall.θarray[indexes[j]] - θ[i])
            #         end
            # end
            return dMdt_sensible,dMdt_latent
end


function liquidmodel(p::PHPSystem)
    sys = deepcopy(p)
    θarrays = sys.liquid.θarrays
    # nondihv_tonondihl = 0.0046206704347650325 # temperary variable to fix different nondimensionlaization

    du = zero.(deepcopy(θarrays))

    # γ = sys.vapor.γ
    Hₗ = sys.liquid.Hₗ
    peri = sys.tube.peri
    α = sys.liquid.α
    Cpₗ = sys.liquid.Cp
    ρₗ = sys.liquid.ρ

    H_rhs = peri / (ρₗ*Cpₗ*Ac)


    θarray_temp_wall = zero.(deepcopy(θarrays))


    for i = 1:length(θarrays)

        Ttuple = getadjacentT(p,i)

        xs = sys.liquid.Xarrays[i];
        dx = mod(xs[2] - xs[1], sys.tube.L)

        θ_wall_inter = sys.mapping.θ_interp_walltoliquid

        fx = map(θ_wall_inter, xs) - θarrays[i]
        du[i] = α .* laplacian(θarrays[i]) ./ dx ./ dx + Hₗ .* fx .* H_rhs

        # same temperature B.C.
        du[i][1] = (Ttuple[1]-θarrays[i][1])/tstep
        du[i][end] = (Ttuple[end]-θarrays[i][end])/tstep

        # # adiabatic temperature B.C.
        # du[i][1] += (θarrays[i][1]-θarrays[i][1])/tstep
        # du[i][end] += (Ttuple[end]-θarrays[i][end])/tstep

        # println(du[i][1])
    end
    return du
end

function getadjacentT(p::PHPSystem,i::Int64)
    Tfirst = PtoT(p.vapor.P[i])

    if p.tube.closedornot == true
        Tlast = i >= length(p.vapor.P) ? PtoT(p.vapor.P[1]) : PtoT(p.vapor.P[i+1])
    else
        Tlast = PtoT(p.vapor.P[i+1])
    end

    (Tfirst,Tlast)
end

# """
#     (Depreciated)
#     get the array of evaporator's heat flux along the wall if the 1D tube wall model is used.
# """
#
# function getwallWearray(Xarray,p::PHPSystem)
#
#     Wearray = zero(deepcopy(Xarray))
#
#
#     for i = 1:length(Xarray)
#         for j = 1:length(p.evaporator.Xe)
#             if ifamongone(Xarray[i],p.evaporator.Xe[j])
#                 Wearray[i] = p.evaporator.We[j]
#             end
#         end
#     end
#
#     return Wearray
# end


"""
    This is a function to get the laplacian of a vector field u

    For now zero-gradient boundary condition is used.

    u    ::  an array
"""


function laplacian(u)
    unew = deepcopy(u)

    dl = ones(length(u)-1)
    dr = dl
    d  = -2*ones(length(u))
    d[1] = -1
    d[end]= -1

    A = Tridiagonal(dl, d, dr)

    unew = A*u

    # #periodic B.C.
    # if periodic
    #
    #     unew[1]   = u[2] - 2*u[1] + u[end]
    #
    #     unew[end] = u[1] - 2*u[end] + u[end-1]
    #
    # else

    # zero gradient B.C.
    # unew[1]   = 0.0
    # unew[end] = 0.0
    # end

    return (unew)
end
#
# function laplacian(u,Ttuple)
#     unew = deepcopy(u)
#
#     dl = ones(length(u)-1)
#     dr = dl
#     d  = -2*ones(length(u))
#
#     A = Tridiagonal(dl, d, dr)
#
#     unew = A*u
#     # zero gradient B.C.
#     unew[1]   = (Ttuple[1]-u[1])*tstep
#     unew[end] = (Ttuple[end]-u[end])*tstep
#     # end
#
#     return (unew)
# end



# q'
function sys_to_heatflux(p::PHPSystem)

    sys = deepcopy(p)

    θarray = sys.wall.θarray
    # γ = sys.vapor.γ
    Hₗ = sys.liquid.Hₗ
    # He = sys.evaporator.He
    k = sys.vapor.k
    δ = sys.vapor.δ
    Hvapor = k ./ δ

    peri = p.tube.peri

    # dx = sys.wall.Xarray[2]-sys.wall.Xarray[1]

    # Xarray = sys.wall.Xarray
    θ_interp_liquidtowall = sys.mapping.θ_interp_liquidtowall
    H_interp_liquidtowall = sys.mapping.H_interp_liquidtowall

    xs =  sys.wall.Xarray

    dθarray = map(θ_interp_liquidtowall, xs) .- θarray
    Harray  = map(H_interp_liquidtowall, xs)

    # qwallarray = -Harray.*dθarray
    qwallarray = -Harray.*dθarray*peri
end


function integrator_to_heatflux(inte)
    sys_to_heatflux(inte.p)
end

function sys_to_Harray(p::PHPSystem)

    sys = deepcopy(p)

    θarray = sys.wall.θarray
    # γ = sys.vapor.γ
    Hₗ = sys.liquid.Hₗ
    # He = sys.evaporator.He
    k = sys.vapor.k
    δ = sys.vapor.δ
    Hvapor = k ./ δ


    # dx = sys.wall.Xarray[2]-sys.wall.Xarray[1]

    # Xarray = sys.wall.Xarray
    # θ_interp_liquidtowall = sys.mapping.θ_interp_liquidtowall
    H_interp_liquidtowall = sys.mapping.H_interp_liquidtowall

    xs =  sys.wall.Xarray

    # dθarray = map(θ_interp_liquidtowall, xs) .- θarray
    Harray  = map(H_interp_liquidtowall, xs)

    Harray
end

function integrator_to_Harray(inte)
    sys_to_Harray(inte.p)
end

getslope(y2,y1,x2,x1) = (y2-y1)/(x2-x1)
