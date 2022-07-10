# module Thermomodel

export dMdtdynamicsmodel,wallmodel,liquidmodel,dynamicsmodel,sys_to_heatflux,sys_to_Harray,integrator_to_heatflux,integrator_to_Harray


function dynamicsmodel(u::Array{Float64,1},p::PHPSystem)


    du = zero(deepcopy(u))

    sys = deepcopy(p)

    Xp = sys.liquid.Xp
    numofliquidslug = length(Xp)

    d = sys.tube.d
    peri = sys.tube.peri
    Ac = sys.tube.Ac
    angle = sys.tube.angle
    g = sys.tube.g


    P = sys.vapor.P
    Eratio = sys.vapor.Eratio
    δstart = sys.vapor.δstart
    δend = sys.vapor.δend
    Lfilm_start = sys.vapor.Lfilm_start
    Lfilm_end = sys.vapor.Lfilm_end

    δarea_start = Ac .* (1 .- ((d .- 2*δstart) ./ d) .^ 2);
    δarea_end = Ac .* (1 .- ((d .- 2*δend) ./ d) .^ 2);

    δdeposit = δfilm
    Adeposit = getAdeposit(sys,δdeposit)


    μₗ = sys.liquid.μₗ

    Lvaporplug = XptoLvaporplug(Xp,sys.tube.L,sys.tube.closedornot)
    Lliquidslug = XptoLliquidslug(Xp,sys.tube.L)
    height = getheight(Xp,sys.tube.L2D,sys.tube.angle)
    Xpvapor = getXpvapor(Xp,sys.tube.L,sys.tube.closedornot)

    # volume_vapor = Lvaporplug .* Ac - Lfilm_start .* δarea_start - Lfilm_end .* δarea_end



    ρₗ = p.liquid.ρ
# get differential equation factors
    lhs = ρₗ*Ac .* Lliquidslug

    rhs_press = Ac ./ lhs

  dXdt_to_stress = -8*μₗ/d
    rhs_dXdt = peri .* Lliquidslug .* dXdt_to_stress ./ lhs


    rhs_g = Ac*ρₗ*g*cos(angle) ./ lhs

    if p.tube.closedornot == false
        println("open loop not supported!")
         return "open loop not supported!"
        end

    if p.tube.closedornot == true

        numofvaporbubble = numofliquidslug

        v_liquid_left=zeros(numofliquidslug)
        v_liquid_right=zeros(numofliquidslug)
        v_vapor_left=zeros(numofvaporbubble)
        v_vapor_right=zeros(numofvaporbubble)

            for i = 1:numofliquidslug

    # temperary v, need a new momentum equation to solve for it
    #  momentum for liquid
                v_momentum = u[2*numofliquidslug+2*i-1]

                # v_liquid_left[i] = v_momentum + v_momentum*Adeposit[i][1]/(Ac-Adeposit[i][1])
                # v_liquid_right[i] = v_momentum + v_momentum*Adeposit[i][end]/(Ac-Adeposit[i][end])

                v_liquid_left[i] = v_momentum
                v_liquid_right[i] = v_momentum

                rhs_dLdt = -v_momentum*(v_liquid_right[i]-v_liquid_left[i])/Lliquidslug[i]

                du[2*i-1] = v_liquid_left[i]
                du[2*i] = v_liquid_right[i]

                # println(i,",",numofliquidslug)

                if i != numofliquidslug
                    du[2*numofliquidslug + 2*i-1] = rhs_dXdt[i]*u[2*numofliquidslug + 2*i-1] + rhs_g[i]*(height[i][1]-height[i][end]) + rhs_press[i] * (P[i]-P[i+1]) + rhs_dLdt
                else
                    du[2*numofliquidslug + 2*i-1] = rhs_dXdt[i]*u[2*numofliquidslug + 2*i-1] + rhs_g[i]*(height[i][1]-height[i][end]) + rhs_press[i] * (P[i]-P[1]) + rhs_dLdt
                end

                du[2*numofliquidslug + 2*i] = du[2*numofliquidslug + 2*i-1]

            end

            # vapor δ
            v_vapor_left[2:end] = v_liquid_right[1:end-1]
            v_vapor_left[1] = v_liquid_right[end]
            v_vapor_right = v_liquid_left

            # println(δstart)
            # println(δend)
            # println(Adeposit)
            # println(Ac)

            A_dδdt_right_vapor = [elem[1] for elem in Adeposit]

            A_dδdt_right_liquid = [elem[2] for elem in Adeposit]
            A_dδdt_left_vapor = zeros(size(A_dδdt_right_vapor))
            A_dδdt_left_vapor[2:end] = A_dδdt_right_liquid[1:end-1]
            A_dδdt_left_vapor[1] = A_dδdt_right_liquid[end]

            dMdt_latent_start,dMdt_sensible,dMdt_latent_end = dMdtdynamicsmodel(Xpvapor,sys)


            F_start = ρₗ .* Ac .* 4 .* δstart .* (d .- δstart) ./ (d^2)
            C_start = ρₗ .* Ac .* 4 .* (d .- 2δstart) ./ (d^2)

            F_end = ρₗ .* Ac .* 4 .* δend .* (d .- δend) ./ (d^2)
            C_end = ρₗ .* Ac .* 4 .* (d .- 2δend) ./ (d^2)

            L0threshold_film = 3e-4
            L0threshold_pure_vapor = 3e-4

            dLdt_start_normal = (-dMdt_latent_start .* (1 .- Eratio) - 0*ρₗ .* A_dδdt_left_vapor  .* v_vapor_left) ./ F_start - v_vapor_left
            dLdt_end_normal = (-dMdt_latent_end .* (1 .- Eratio) + 0*ρₗ .* A_dδdt_right_vapor .* v_vapor_right) ./ F_end + v_vapor_right

            # dLdt_start_normal = (-dMdt_latent_start .* (1 .- Eratio) - ρₗ .* A_dδdt_left_vapor  .* v_vapor_left) ./ F_start
            # dLdt_end_normal = (-dMdt_latent_end .* (1 .- Eratio) + ρₗ .* A_dδdt_right_vapor .* v_vapor_right) ./ F_end 

            he_start_short = Bool.(heaviside.(-Lfilm_start .+ L0threshold_film))
            he_end_short = Bool.(heaviside.(-Lfilm_end .+ L0threshold_film))
            he_start_positive = Bool.(heaviside.(dLdt_start_normal))
            he_end_positive = Bool.(heaviside.(dLdt_end_normal))
            he_meet= Bool.(heaviside.(-Lvaporplug .+ Lfilm_start .+ Lfilm_end .+ L0threshold_pure_vapor))

            dLdt_start = zeros(numofvaporbubble)
            dLdt_end = zeros(numofvaporbubble)
        

            #  println(length(Eratio))
           
            # dLdt_start_case1 = 0
            # dLdt_start_case2 = dLdt_start_normal
            # dLdt_start_case3 = -v_vapor_left
            # dLdt_start_case4 = v_vapor_right .- v_vapor_left

            for i = 1:numofvaporbubble
                if he_meet[i]
                    if he_start_short[i]
                        dLdt_start[i] = 0
                    elseif he_start_positive[i]
                        dLdt_start[i] = he_end_short[i] ? v_vapor_right[i] .- v_vapor_left[i] : -v_vapor_left[i]
                        if i == 1 
                            # println(-v_vapor_left[i],v_vapor_right[i] .- v_vapor_left[i],dLdt_start[i])

                        end
                    else
                        dLdt_start[i] = minimum([dLdt_start_normal[i],v_vapor_right[i] .- v_vapor_left[i], -v_vapor_left[i]])
                    end
                elseif he_start_short[i]
                    dLdt_start[i] = he_start_positive[i] ? dLdt_start_normal[i] : 0
                else
                    dLdt_start[i] = dLdt_start_normal[i]
                end
            end

            # println(F_start )

            # dLdt_start = he_start_short .* dLdt_start_case1 + ((1 .- he_start_short) .* (1 .- he_meet) .+  he_start_short .* (1 .- he_meet) .* he_start_positive .+ (1 .- he_start_short) .* he_meet .* (1 .- he_start_positive)) .* dLdt_start_case2 + (1 .- he_start_short) .* he_meet .* he_start_positive .* (1 .- he_end_short) .* dLdt_start_case3 + (1 .- he_start_short) .* he_meet .* he_start_positive .* he_end_short .* dLdt_start_case4

            # dLdt_end_case1 = 0
            # dLdt_end_case2 = dLdt_end_normal
            # dLdt_end_case3 = v_vapor_right
            # dLdt_end_case4 = -(v_vapor_left .- v_vapor_right)

            for i = 1:numofvaporbubble
                if he_meet[i]
                    if he_end_short[i]
                        dLdt_end[i] = 0
                    elseif he_end_positive[i]
                        dLdt_end[i] = he_start_short[i] ? v_vapor_right[i] .- v_vapor_left[i] : v_vapor_right[i]
                    else
                        dLdt_end[i] = minimum([dLdt_end_normal[i],v_vapor_right[i] .- v_vapor_left[i],v_vapor_right[i]])
                    end
                elseif he_end_short[i]
                    dLdt_end[i] = he_end_positive[i] ? dLdt_end_normal[i] : 0
                else
                    dLdt_end[i] = dLdt_end_normal[i]
                end
            end

            # dLdt_end = he_end_short .* dLdt_end_case1 + ((1 .- he_end_short) .* (1 .- he_meet) .+ he_end_short .* (1 .- he_meet) .* he_end_positive .+ (1 .- he_end_short) .* he_meet .* (1 .- he_end_positive)).* dLdt_end_case2 + (1 .- he_end_short) .* he_meet .* he_end_positive .* (1 .- he_start_short) .* dLdt_end_case3 + (1 .- he_end_short) .* he_meet .* he_end_positive .* he_start_short .* dLdt_end_case4

            # dLdt_start = (-dMdt_latent_start .* (1 .- Eratio) - ρₗ .* A_dδdt_left_vapor  .* v_vapor_left) ./ F_start .* heaviside_L_start .+  (v_vapor_right .- v_vapor_left) .* heaviside_L_total_start
            # dLdt_end = (-dMdt_latent_end .* (1 .- Eratio) + ρₗ .* A_dδdt_right_vapor .* v_vapor_right) ./ F_end .* heaviside_L_end .+ (v_vapor_right .- v_vapor_left) .* heaviside_L_total_end
            
            # dLdt_start = (-dMdt_latent_start .* (1 .- Eratio) - ρₗ .* A_dδdt_left_vapor  .* v_vapor_left) ./ F_start .* heaviside_L_start
            # dLdt_end = (-dMdt_latent_end .* (1 .- Eratio) + ρₗ .* A_dδdt_right_vapor .* v_vapor_right) ./ F_end .* heaviside_L_end

            dδdt_start_normal = (-dMdt_latent_start .* Eratio) ./ (C_start .* Lfilm_start)
            dδdt_end_normal = (-dMdt_latent_end .* Eratio) ./ (C_end .* Lfilm_end)

            he_dδdt_start_positive = dδdt_start_normal .> 0
            he_dδdt_end_positive = dδdt_end_normal .> 0
            he_dδdt_start_toobig = δstart .> δfilm*3
            he_dδdt_start_toosmall = δstart .< δfilm/3
            he_dδdt_end_toobig = δend .> δfilm*3
            he_dδdt_end_toosmall = δend .< δfilm/3

            dδdt_start = (1 .- (he_dδdt_start_toosmall .* (1 .- he_dδdt_start_positive) .+ he_dδdt_start_toobig .* he_dδdt_start_positive)) .* dδdt_start_normal
            # println((1 .- div.((he_dδdt_end_toosmall .* (1 .- he_dδdt_end_positive) .+ he_dδdt_end_toobig .* he_dδdt_end_positive),2)))
            dδdt_end = (1 .- ((he_dδdt_end_toosmall .* (1 .- he_dδdt_end_positive) .+ he_dδdt_end_toobig .* he_dδdt_end_positive))) .* dδdt_end_normal

            # dδdt_start = (δstart < [δfilm/3] && dδdt_start_normal < [0]) || (δstart > [δfilm*3] && dδdt_start_normal > [0]) ? 0.0 : dδdt_start_normal
            # dδdt_end = (δend < [δfilm/3] && dδdt_end_normal < [0]) || (δend > [δfilm*3] && dδdt_end_normal > [0]) ? 0.0 : dδdt_end_normal



            du[4*numofliquidslug+1:5*numofliquidslug] .= dMdt_latent_start+dMdt_sensible+dMdt_latent_end
            du[5*numofliquidslug+1:6*numofliquidslug] .= dδdt_start # equals to 0 for now
            du[6*numofliquidslug+1:7*numofliquidslug] .= dδdt_end # equals to 0 for now
            du[7*numofliquidslug+1:8*numofliquidslug] .= dLdt_start # equals to 0 for now
            du[8*numofliquidslug+1:9*numofliquidslug] .= dLdt_end # equals to 0 for now

            # println(δend[1])
            # testindex = 10
            # println(he_start_short[testindex])
            # println(he_end_short[testindex])
            # println(he_start_positive[testindex])
            # println(he_end_positive[testindex])
            # println(he_meet[testindex])
            # println(Xpvapor)
            # println(δstart)
            # println(Xpvapor[testindex])
            # println(v_vapor_right[testindex])
            # println(v_vapor_left[testindex])
            # println(dMdt_latent_start[testindex])
            # # println(Lfilm_start[testindex])
            # println(\delta     film_end[testindex])
            # # println(Lvaporplug[testindex])
            # # # # println(he_start_short)
            # println(dLdt_start[testindex])
            # println(dLdt_end[testindex])
            # println(dLdt_end[5])
            # # println(dLdt_start_case2[67])
            # # println(dLdt_end[67])
            # # println(dLdt_end_case2[67])
            # # println(Lfilm_end)
            # # println(dLdt_end)
            # # # # println(dLdt_start_case1[21])
            # # println(dLdt_end[6])
            # # # println(dLdt_end_case4[21])
            # # println(v_vapor_right[6])
            # # println(v_vapor_left[6])

            return du
    end


end

function dMdtdynamicsmodel(Xpvapor::Array{Tuple{Float64,Float64},1},sys::PHPSystem)

    dMdt_sensible=zeros(length(Xpvapor))
    dMdt_latent_start=zeros(length(Xpvapor))
    dMdt_latent_end=zeros(length(Xpvapor))

    L = sys.tube.L

    peri = sys.tube.peri
    Ac = sys.tube.Ac
    k = sys.vapor.k

    Lfilm_start = sys.vapor.Lfilm_start
    Lfilm_end = sys.vapor.Lfilm_end

    P = sys.vapor.P
    θarrays = sys.liquid.θarrays
    Xarrays = sys.liquid.Xarrays

    θ = PtoT.(P)

    Hfg = PtoHfg.(P)

    
    H_interp = sys.mapping.H_interp_liquidtowall
    θ_wall_interp = sys.mapping.θ_interp_walltoliquid

    dx_wall = sys.wall.Xarray[2]-sys.wall.Xarray[1]

    Lvapor = XptoLvaporplug(sys.liquid.Xp,sys.tube.L,sys.tube.closedornot)
    Lvapor_pure = max.(Lvapor - Lfilm_start - Lfilm_end,0.0)

    for i = 1:length(Xpvapor)
        Nstart = Int64(max(2 , div(Lfilm_start[i],dx_wall)))
        heatflux_start = quad_trap(H_interp,θ_wall_interp,θ[i],Xpvapor[i][1],mod(Xpvapor[i][1]+Lfilm_start[i],L),L,Nstart)

        Nvapor_pure = Int64(max(2 , div(Lvapor_pure[i],dx_wall)))
        heatflux_pure_vapor = quad_trap(H_interp,θ_wall_interp,θ[i],mod(Xpvapor[i][1]+Lfilm_start[i],L),mod(Xpvapor[i][2]-Lfilm_end[i],L),L,Nvapor_pure)

        Nend = Int64(max(2 , div(Lfilm_end[i],dx_wall)))
        heatflux_end = quad_trap(H_interp,θ_wall_interp,θ[i],mod(Xpvapor[i][2]-Lfilm_end[i],L),Xpvapor[i][2],L,Nend)

        slope_r = getslope(θarrays[i][2],θarrays[i][1],Xarrays[i][2],Xarrays[i][1])
        slope_l = (i == 1) ? getslope(θarrays[1][end],θarrays[1][end-1],Xarrays[1][end],Xarrays[1][end-1]) : getslope(θarrays[i-1][end],θarrays[i-1][end-1],Xarrays[i-1][end],Xarrays[i-1][end-1])


        axial_rhs_end = Ac*k*slope_r /Hfg[i]
        axial_rhs_start = Ac*k*(-slope_l) /Hfg[i]

        dMdt_latent_start[i] = heatflux_start*peri/Hfg[i] + axial_rhs_start

        dMdt_sensible[i] = heatflux_pure_vapor*peri/Hfg[i]

        dMdt_latent_end[i] = heatflux_end*peri/Hfg[i] + axial_rhs_end
        
        end
            return dMdt_latent_start,dMdt_sensible,dMdt_latent_end
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
    # k = sys.vapor.k
    # δ = sys.vapor.δ
    # Hvapor = k ./ δ

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

    # println(Harray)
end


function integrator_to_heatflux(inte)
    sys_to_heatflux(inte.p)
end

function sys_to_Harray(p::PHPSystem)

    sys = deepcopy(p)

    # θarray = sys.wall.θarray
    # γ = sys.vapor.γ
    # Hₗ = sys.liquid.Hₗ
    # He = sys.evaporator.He
    # k = sys.vapor.k
    # δ = sys.vapor.δ
    # Hvapor = k ./ δ


    # dx = sys.wall.Xarray[2]-sys.wall.Xarray[1]

    # Xarray = sys.wall.Xarray
    # θ_interp_liquidtowall = sys.mapping.θ_interp_liquidtowall
    H_interp_liquidtowall = sys.mapping.H_interp_liquidtowall

    xs =  sys.wall.Xarray

    # dθarray = map(θ_interp_liquidtowall, xs) .- θarray
    Harray  = map(H_interp_liquidtowall, xs)

    Harray
end


function quad_trap(H_interp,θ_interp,θvapor_one, a,b,L,N) 
    h = mod((b-a),L)/N
    
    # println(h)

    H_interp(a)*(θ_interp(a)-θvapor_one)
    # println(b)

    int = h * ( H_interp(a)*(θ_interp(a)-θvapor_one) + H_interp(b)*(θ_interp(b)-θvapor_one) ) / 2
    for k=1:N-1
        xk = mod((b-a),L) * k/N + a
        int = int + h*H_interp(mod(xk,L))*(θ_interp(mod(xk,L))-θvapor_one)
    end
    return int
end

function integrator_to_Harray(inte)
    sys_to_Harray(inte.p)
end

getslope(y2,y1,x2,x1) = (y2-y1)/(x2-x1)

heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.5)))
