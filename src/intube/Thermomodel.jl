using Statistics

export dMdtdynamicsmodel,dynamicsmodel_steadyfilm,wallmodel,liquidmodel,dynamicsmodel,sys_to_heatflux,sys_to_Harray,integrator_to_heatflux,integrator_to_Harray


function dynamicsmodel(u::Array{Float64,1},p::PHPSystem)


    du = zero(deepcopy(u))

    sys = deepcopy(p)

    σ = sys.liquid.σ
    μₗ = sys.liquid.μₗ
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

    V = [elem[2] for elem in sys.liquid.dXdt]
    Vavg = mean(abs.(V))
    Ca = getCa.(μₗ,σ,Vavg)
    
    ad_fac = Main.ad_fac
    δdeposit = Catoδ(d,Ca,adjust_factor=ad_fac)
    # δdeposit = δfilm
    Adeposit = getAdeposit(sys,δdeposit)


    Lvaporplug = XptoLvaporplug(Xp,sys.tube.L,sys.tube.closedornot)
    Lliquidslug = XptoLliquidslug(Xp,sys.tube.L)
    height = getheight(Xp,sys.tube.L2D,sys.tube.angle)
    Xpvapor = getXpvapor(Xp,sys.tube.L,sys.tube.closedornot)


    ρₗ = p.liquid.ρ
# get differential equation factors
    lhs = ρₗ*Ac .* Lliquidslug

    rhs_press = Ac ./ lhs

    Re_list = ρₗ .* abs.(V) .* d ./ μₗ
    # println(V)
    f_coefficient = f_churchill.(Re_list .+ 1e-4)
    dXdt_to_stress = -0.125 .* f_coefficient .* ρₗ .* V .* abs.(V)
    # dXdt_to_stress = -8*μₗ/d

    rhs_dXdt = peri .* Lliquidslug .* dXdt_to_stress ./ lhs

    # println(dXdt_to_stress)
    # println(-8*μₗ/d .* V)

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
            v_momentum = u[2*numofliquidslug+2*i-1]

            v_liquid_left[i] = v_momentum #approximation
            v_liquid_right[i] = v_momentum #approximation
        end

            # vapor δ
            v_vapor_left[2:end] = v_liquid_right[1:end-1]
            v_vapor_left[1] = v_liquid_right[end]
            v_vapor_right = v_liquid_left

            A_dδdt_right_vapor = [elem[1] for elem in Adeposit]
            A_dδdt_left_liquid = A_dδdt_right_vapor
            A_dδdt_right_liquid = [elem[2] for elem in Adeposit]
            A_dδdt_left_vapor = zeros(size(A_dδdt_right_vapor))
            A_dδdt_left_vapor[2:end] = A_dδdt_right_liquid[1:end-1]
            A_dδdt_left_vapor[1] = A_dδdt_right_liquid[end]

            dMdt_latent_start,dMdt_sensible,dMdt_latent_end = dMdtdynamicsmodel(Xpvapor,sys)
            dMdt_latent_start_positive,dMdt_sensible_positive,dMdt_latent_end_positive = dMdtdynamicsmodel_positive(Xpvapor,sys)

            # println(dMdt_latent_start[4:6])
            # println(dMdt_sensible[4:6])
            # println(dMdt_latent_end[1:10])
            # println(dMdt_latent_start_positive[4:6])
            # println(dMdt_latent_end_positive[4:6])

            F_start = ρₗ .* Ac .* 4 .* δstart .* (d .- δstart) ./ (d^2)
            C_start = ρₗ .* Ac .* 4 .* (d .- 2δstart) ./ (d^2)

            F_end = ρₗ .* Ac .* 4 .* δend .* (d .- δend) ./ (d^2)
            C_end = ρₗ .* Ac .* 4 .* (d .- 2δend) ./ (d^2)

            L0threshold_film = 4e-4
            L0threshold_pure_vapor = 1e-3

            dLdt_start_normal = (-dMdt_latent_start_positive .* Eratio) ./ F_start - v_vapor_left
            dLdt_end_normal = (-dMdt_latent_end_positive .* Eratio) ./ F_end + v_vapor_right

            # println(maximum(abs.(dLdt_start_normal + v_vapor_left)))
            # println(maximum(abs.(dLdt_end_normal - v_vapor_right)))
            # println(v_vapor_left)


            he_start_short = Bool.(heaviside.(-Lfilm_start .+ L0threshold_film))
            he_end_short = Bool.(heaviside.(-Lfilm_end .+ L0threshold_film))
            he_start_positive = Bool.(heaviside.(dLdt_start_normal))
            he_end_positive = Bool.(heaviside.(dLdt_end_normal))
            he_meet= Bool.(heaviside.(-Lvaporplug .+ Lfilm_start .+ Lfilm_end .+ L0threshold_pure_vapor))

            dLdt_start = zeros(numofvaporbubble)
            dLdt_end = zeros(numofvaporbubble)
        
            for i = 1:numofvaporbubble
                if he_meet[i]
                    if he_start_short[i]
                        dLdt_start[i] = 0
                    elseif he_start_positive[i]
                        dLdt_start[i] = he_end_short[i] ? v_vapor_right[i] .- v_vapor_left[i] : -v_vapor_left[i]
                    else
                        dLdt_start[i] = minimum([dLdt_start_normal[i],v_vapor_right[i] .- v_vapor_left[i], -v_vapor_left[i]])
                    end
                elseif he_start_short[i]
                    dLdt_start[i] = he_start_positive[i] ? dLdt_start_normal[i] : 0
                else
                    dLdt_start[i] = dLdt_start_normal[i]
                end
            end

            # aa = maximum(abs.((dLdt_start + v_vapor_left) .* (1 .- he_meet) .* (1 .- he_start_short)))
            # if aa > 1e-8 
            #     println(aa)
            # end
     
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

            dδdt_start_normal = (-dMdt_latent_start .- (-dMdt_latent_start_positive .* Eratio)) ./ (C_start .* Lfilm_start)  + (- ρₗ .* A_dδdt_left_vapor .* v_vapor_left + F_start .* v_vapor_left) ./ (C_start .* Lfilm_start) 
            dδdt_end_normal = (-dMdt_latent_end .- (-dMdt_latent_end_positive .* Eratio)) ./ (C_end .* Lfilm_end) - (- ρₗ .* A_dδdt_right_vapor .* v_vapor_right + F_end .* v_vapor_right) ./ (C_end .* Lfilm_end)

            he_dδdt_start_positive = dδdt_start_normal .> 0
            he_dδdt_end_positive = dδdt_end_normal .> 0
            he_dδdt_start_toobig = δstart .> 1e-4
            he_dδdt_start_toosmall = δstart .< 3e-6
            he_dδdt_end_toobig = δend .> 1e-4
            he_dδdt_end_toosmall = δend .< 3e-6

            dδdt_start = (1 .- (he_dδdt_start_toosmall .* (1 .- he_dδdt_start_positive) .+ he_dδdt_start_toobig .* he_dδdt_start_positive)) .* dδdt_start_normal
            # println((1 .- div.((he_dδdt_end_toosmall .* (1 .- he_dδdt_end_positive) .+ he_dδdt_end_toobig .* he_dδdt_end_positive),2)))
            dδdt_end = (1 .- ((he_dδdt_end_toosmall .* (1 .- he_dδdt_end_positive) .+ he_dδdt_end_toobig .* he_dδdt_end_positive))) .* dδdt_end_normal

            Astart = getδarea(Ac,d,δstart)
            Aend = getδarea(Ac,d,δend)

            dMdt_film_start = ρₗ .* Lfilm_start .* 4Ac ./ d .* (1 .- 2 .* δstart ./ d) .* dδdt_start .+ ρₗ .* Astart .* dLdt_start .- dMdt_latent_start
            dMdt_film_end = ρₗ .* Lfilm_end .* 4Ac ./ d .* (1 .- 2 .* δend ./ d) .* dδdt_end .+ ρₗ .* Aend .* dLdt_end .- dMdt_latent_end

            Nvapor = length(P)
            loop_plus_index = [2:Nvapor;1]
            for i = 1:numofliquidslug
                du[2*i-1] = v_liquid_left[i] + dMdt_film_end[i] ./ ρₗ ./ (Ac .- 2Aend[i])
                du[2*i] = v_liquid_right[i] + dMdt_film_start[loop_plus_index[i]] ./ ρₗ ./ (Ac .- 2Astart[loop_plus_index[i]])

                v_momentum = v_liquid_right[i]
                rhs_dLdt = -v_momentum*(du[2*i]-du[2*i-1])/Lliquidslug[i]
            
                if i != numofliquidslug
                    du[2*numofliquidslug + 2*i-1] = rhs_dXdt[i] + rhs_g[i]*(height[i][1]-height[i][end]) + rhs_press[i] * (P[i]-P[i+1]) + rhs_dLdt
                else
                    du[2*numofliquidslug + 2*i-1] = rhs_dXdt[i] + rhs_g[i]*(height[i][1]-height[i][end]) + rhs_press[i] * (P[i]-P[1]) + rhs_dLdt
                end                
            
                du[2*numofliquidslug + 2*i] = du[2*numofliquidslug + 2*i-1]
            end
            

            du[4*numofliquidslug+1:5*numofliquidslug] .= dMdt_latent_start+dMdt_latent_end
            # du[4*numofliquidslug+1:5*numofliquidslug] .= dMdt_latent_start+dMdt_sensible+dMdt_latent_end
            du[5*numofliquidslug+1:6*numofliquidslug] .= dδdt_start # equals to 0 for now
            du[6*numofliquidslug+1:7*numofliquidslug] .= dδdt_end # equals to 0 for now
            du[7*numofliquidslug+1:8*numofliquidslug] .= dLdt_start # equals to 0 for now
            du[8*numofliquidslug+1:9*numofliquidslug] .= dLdt_end # equals to 0 for now

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

function dMdtdynamicsmodel_positive(Xpvapor::Array{Tuple{Float64,Float64},1},sys::PHPSystem)

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
        heatflux_start = quad_trap_positive(H_interp,θ_wall_interp,θ[i],Xpvapor[i][1],mod(Xpvapor[i][1]+Lfilm_start[i],L),L,Nstart)

        Nvapor_pure = Int64(max(2 , div(Lvapor_pure[i],dx_wall)))
        heatflux_pure_vapor = quad_trap_positive(H_interp,θ_wall_interp,θ[i],mod(Xpvapor[i][1]+Lfilm_start[i],L),mod(Xpvapor[i][2]-Lfilm_end[i],L),L,Nvapor_pure)

        Nend = Int64(max(2 , div(Lfilm_end[i],dx_wall)))
        heatflux_end = quad_trap_positive(H_interp,θ_wall_interp,θ[i],mod(Xpvapor[i][2]-Lfilm_end[i],L),Xpvapor[i][2],L,Nend)

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
    tstep = Main.tstep
    sys = deepcopy(p)
    θarrays = sys.liquid.θarrays
    # nondihv_tonondihl = 0.0046206704347650325 # temperary variable to fix different nondimensionlaization

    du = zero.(deepcopy(θarrays))

    # γ = sys.vapor.γ
    Ac = p.tube.Ac

    Hₗ = sys.liquid.Hₗ
    peri = sys.tube.peri
    α = sys.liquid.α
    Cpₗ = sys.liquid.Cp
    ρₗ = sys.liquid.ρ

    H_rhs = peri / (ρₗ*Cpₗ*Ac)


    θarray_temp_wall = zero.(deepcopy(θarrays))


    for i in eachindex(θarrays)

        Ttuple = getadjacentT(p,i)

        xs = sys.liquid.Xarrays[i];
        dx = mod(xs[2] - xs[1], sys.tube.L)

        θ_wall_inter = sys.mapping.θ_interp_walltoliquid

        fx = map(θ_wall_inter, xs) - θarrays[i]
        du[i] = α .* laplacian(θarrays[i]) ./ dx ./ dx + Hₗ .* fx .* H_rhs
        # du[i] = Hₗ .* fx .* H_rhs

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

    return (unew)
end


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

    H_interp_liquidtowall = sys.mapping.H_interp_liquidtowall

    xs =  sys.wall.Xarray

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
    return (b>a) ? int : 0
end

function quad_trap_positive(H_interp,θ_interp,θvapor_one, a,b,L,N) 
    h = mod((b-a),L)/N
    
    # println(h)

    H_interp(a)*(θ_interp(a)-θvapor_one)
    # println(b)

    int = maximum([h * ( H_interp(a)*(θ_interp(a)-θvapor_one) + H_interp(b)*(θ_interp(b)-θvapor_one) ) / 2, 0.0])
    for k=1:N-1
        xk = mod((b-a),L) * k/N + a
        int = int + maximum([h*H_interp(mod(xk,L))*(θ_interp(mod(xk,L))-θvapor_one),0.0])
    end
    return (b>a) ? int : 0
end

function integrator_to_Harray(inte)
    sys_to_Harray(inte.p)
end

getslope(y2,y1,x2,x1) = (y2-y1)/(x2-x1)

heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.0)))

function dynamicsmodel_steadyfilm(u::Array{Float64,1},p::PHPSystem)


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

    # temporary!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    δfilm = Main.δfilm

    δdeposit = δfilm
    Adeposit = getAdeposit(sys,δdeposit)


    μₗ = sys.liquid.μₗ

    Lvaporplug = XptoLvaporplug(Xp,sys.tube.L,sys.tube.closedornot)
    Lliquidslug = XptoLliquidslug(Xp,sys.tube.L)
    height = getheight(Xp,sys.tube.L2D,sys.tube.angle)
    Xpvapor = getXpvapor(Xp,sys.tube.L,sys.tube.closedornot)


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

    # temperary v, need a new momentum equation to solve for it, but now I have to let it be, for liquid mass conservation
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
                    du[2*numofliquidslug + 2*i-1] = (rhs_dXdt[i]*u[2*numofliquidslug + 2*i-1] + rhs_g[i]*(height[i][1]-height[i][end]) + rhs_press[i] * (P[i]-P[i+1]) + rhs_dLdt) .* 0.0
                else
                    du[2*numofliquidslug + 2*i-1] = (rhs_dXdt[i]*u[2*numofliquidslug + 2*i-1] + rhs_g[i]*(height[i][1]-height[i][end]) + rhs_press[i] * (P[i]-P[1]) + rhs_dLdt) .* 0.0
                end

                du[2*numofliquidslug + 2*i] = du[2*numofliquidslug + 2*i-1]

            end

            # vapor δ
            v_vapor_left[2:end] = v_liquid_right[1:end-1]
            v_vapor_left[1] = v_liquid_right[end]
            v_vapor_right = v_liquid_left

            A_dδdt_right_vapor = [elem[1] for elem in Adeposit]

            A_dδdt_right_liquid = [elem[2] for elem in Adeposit]
            A_dδdt_left_vapor = zeros(size(A_dδdt_right_vapor))
            A_dδdt_left_vapor[2:end] = A_dδdt_right_liquid[1:end-1]
            A_dδdt_left_vapor[1] = A_dδdt_right_liquid[end]

            dMdt_latent_start,dMdt_sensible,dMdt_latent_end = dMdtdynamicsmodel(Xpvapor,sys)
            dMdt_latent_start_positive,dMdt_sensible_positive,dMdt_latent_end_positive = dMdtdynamicsmodel_positive(Xpvapor,sys)

            # println(dMdt_latent_start[4:6])
            # println(dMdt_sensible[4:6])
            # println(dMdt_latent_end[1:10])
            # println(dMdt_latent_start_positive[4:6])
            # println(dMdt_latent_end_positive[4:6])

            F_start = ρₗ .* Ac .* 4 .* δstart .* (d .- δstart) ./ (d^2)
            C_start = ρₗ .* Ac .* 4 .* (d .- 2δstart) ./ (d^2)

            F_end = ρₗ .* Ac .* 4 .* δend .* (d .- δend) ./ (d^2)
            C_end = ρₗ .* Ac .* 4 .* (d .- 2δend) ./ (d^2)

            L0threshold_film = 4e-4
            L0threshold_pure_vapor = 1e-3

            dLdt_start_normal = (-dMdt_latent_start_positive .* Eratio) ./ F_start - v_vapor_left
            dLdt_end_normal = (-dMdt_latent_end_positive .* Eratio) ./ F_end + v_vapor_right

            # println(maximum(abs.(dLdt_start_normal + v_vapor_left)))
            # println(maximum(abs.(dLdt_end_normal - v_vapor_right)))
            # println(v_vapor_left)


            he_start_short = Bool.(heaviside.(-Lfilm_start .+ L0threshold_film))
            he_end_short = Bool.(heaviside.(-Lfilm_end .+ L0threshold_film))
            he_start_positive = Bool.(heaviside.(dLdt_start_normal))
            he_end_positive = Bool.(heaviside.(dLdt_end_normal))
            he_meet= Bool.(heaviside.(-Lvaporplug .+ Lfilm_start .+ Lfilm_end .+ L0threshold_pure_vapor))

            dLdt_start = zeros(numofvaporbubble)
            dLdt_end = zeros(numofvaporbubble)
        
            for i = 1:numofvaporbubble
                if he_meet[i]
                    if he_start_short[i]
                        dLdt_start[i] = 0
                    elseif he_start_positive[i]
                        dLdt_start[i] = he_end_short[i] ? v_vapor_right[i] .- v_vapor_left[i] : -v_vapor_left[i]
                    else
                        dLdt_start[i] = minimum([dLdt_start_normal[i],v_vapor_right[i] .- v_vapor_left[i], -v_vapor_left[i]])
                    end
                elseif he_start_short[i]
                    dLdt_start[i] = he_start_positive[i] ? dLdt_start_normal[i] : 0
                else
                    dLdt_start[i] = dLdt_start_normal[i]
                end
            end

            # aa = maximum(abs.((dLdt_start + v_vapor_left) .* (1 .- he_meet) .* (1 .- he_start_short)))
            # if aa > 1e-8 
            #     println(aa)
            # end
     
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

            # println(he_end_short[5])
            # println(he_end_positive[5])
            # println(he_meet[5])
            # println(Lfilm_end[4:6])
            # println(dLdt_end[4:6])
            # println(dLdt_end_normal[4:6])
            # println(δend[4:6])

            # dδdt_start_normal = (-dMdt_latent_start .- (-dMdt_latent_start_positive .* Eratio)) ./ (C_start .* Lfilm_start) 
            # dδdt_end_normal = (-dMdt_latent_end .- (-dMdt_latent_end_positive .* Eratio)) ./ (C_end .* Lfilm_end)

            dδdt_start_normal = (-dMdt_latent_start .- (-dMdt_latent_start_positive .* Eratio)) ./ (C_start .* Lfilm_start)  + (- ρₗ .* A_dδdt_left_vapor .* v_vapor_left + F_start .* v_vapor_left) ./ (C_start .* Lfilm_start) 
            dδdt_end_normal = (-dMdt_latent_end .- (-dMdt_latent_end_positive .* Eratio)) ./ (C_end .* Lfilm_end) - (- ρₗ .* A_dδdt_right_vapor .* v_vapor_right + F_end .* v_vapor_right) ./ (C_end .* Lfilm_end)

            he_dδdt_start_positive = dδdt_start_normal .> 0
            he_dδdt_end_positive = dδdt_end_normal .> 0
            he_dδdt_start_toobig = δstart .> 1e-4
            he_dδdt_start_toosmall = δstart .< 3e-6
            he_dδdt_end_toobig = δend .> 1e-4
            he_dδdt_end_toosmall = δend .< 3e-6

            dδdt_start = (1 .- (he_dδdt_start_toosmall .* (1 .- he_dδdt_start_positive) .+ he_dδdt_start_toobig .* he_dδdt_start_positive)) .* dδdt_start_normal
            # println((1 .- div.((he_dδdt_end_toosmall .* (1 .- he_dδdt_end_positive) .+ he_dδdt_end_toobig .* he_dδdt_end_positive),2)))
            dδdt_end = (1 .- ((he_dδdt_end_toosmall .* (1 .- he_dδdt_end_positive) .+ he_dδdt_end_toobig .* he_dδdt_end_positive))) .* dδdt_end_normal

            volume_vapor = Lvaporplug .* Ac - Lfilm_start .* δarea_start - Lfilm_end .* δarea_end;

            du[4*numofliquidslug+1:5*numofliquidslug] .= (dMdt_latent_start+dMdt_sensible+dMdt_latent_end).*0.0
            du[5*numofliquidslug+1:6*numofliquidslug] .= dδdt_start # equals to 0 for now
            du[6*numofliquidslug+1:7*numofliquidslug] .= dδdt_end # equals to 0 for now
            du[7*numofliquidslug+1:8*numofliquidslug] .= dLdt_start # equals to 0 for now
            du[8*numofliquidslug+1:9*numofliquidslug] .= dLdt_end # equals to 0 for now

            return du
    end


end