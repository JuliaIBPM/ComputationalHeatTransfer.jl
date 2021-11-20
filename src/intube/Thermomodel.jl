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

    μₗ = sys.liquid.μₗ

    Lvaporplug = XptoLvaporplug(Xp,sys.tube.L,sys.tube.closedornot)
    Lliquidslug = XptoLliquidslug(Xp,sys.tube.L)
    height = getheight(Xp,sys.tube.L2D,sys.tube.angle)
    Xpvapor = getXpvapor(Xp,sys.tube.L,sys.tube.closedornot)


    ρ = M ./ Lvaporplug ./ Ac
    P = DtoP.(ρ)
    θ = PtoT.(P)




# get differential equation factors
    lhs = ρ*Ac .* Lliquidslug

    rhs_press = Ac ./ lhs

    dXdt_to_stress = -8*μₗ/d
    rhs_dXdt = peri .* Lliquidslug .* dXdt_to_stress ./ lhs

    rhs_g = Ac*ρ*g*cos(angle) ./ lhs


    # P = nondi_DtoP.(M./Lvaporplug)
    # θ = nondi_PtoT.(P)
    # P = real.((M./Lvaporplug .+ 0im).^(γ))
    # θ = real.((P .+ 0im).^((γ-1)/γ))


    if p.tube.closedornot == false

    for i = 1:numofliquidslug
        du[2*i-1] = u[2*numofliquidslug+2*i-1]
        du[2*i] = du[2*i-1]
        du[2*numofliquidslug + 2*i-1] = rhs_dXdt[i]*u[2*numofliquidslug + 2*i-1] +rhs_g[i]*(height[i][1]-height[i][end]) + rhs_press[i] * (P[i]-P[i+1])
        # use ℘L and ω for now, in the future change them to one value rather than a array.
        du[2*numofliquidslug + 2*i] = du[2*numofliquidslug + 2*i-1]

    end

# not sure which one is better
    du[4*numofliquidslug+1:5*numofliquidslug+1] .= dMdtdynamicsmodel(Xpvapor,sys)
    du[5*numofliquidslug+2:end] .= [0.0]

    return du
    end

if p.tube.closedornot == true

        for i = 1:numofliquidslug
            du[2*i-1] = u[2*numofliquidslug+2*i-1]
            du[2*i] = du[2*i-1]

            if i != numofliquidslug
                du[2*numofliquidslug + 2*i-1] = rhs_dXdt[i]*u[2*numofliquidslug + 2*i-1] + rhs_g[i]*(height[i][1]-height[i][end]) + rhs_press[i] * (P[i]-P[i+1])
            else
                du[2*numofliquidslug + 2*i-1] = rhs_dXdt[i]*u[2*numofliquidslug + 2*i-1] + rhs_g[i]*(height[i][1]-height[i][end]) + rhs_press[i] * (P[i]-P[1])
            end
            # use ℘L and ω for now, in the future change them to one value rather than a array.
            du[2*numofliquidslug + 2*i] = du[2*numofliquidslug + 2*i-1]

        end

    # not sure which one is better
        du[4*numofliquidslug+1:5*numofliquidslug] .= dMdtdynamicsmodel(Xpvapor,sys)
        du[5*numofliquidslug+1:end] .= [0.0] # equals to 0 for now

        return du
end


end

function dMdtdynamicsmodel(Xpvapor::Array{Tuple{Float64,Float64},1},sys::PHPSystem)

    dMdt=zeros(length(Xpvapor))

    # get Hvapor
    peri = sys.tube.peri
    k = sys.vapor.k
    δ = sys.vapor.δ
    Hvapor = k ./ δ

    #get θ
    P = sys.vapor.P
    # γ = sys.vapor.γ

    θ = PtoT.(P)
    # θ = nondi_PtoT.(P)
    # θ = real.((P .+ 0im).^((γ-1)/γ)) # isentropic

    Hfg = PtoHfg.(P)

    # println(Hfg)

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
        dMdt[i] = (sum(fx) * dx_vapor) * (Hvapor[i] *peri/Hfg[i])

        # println(dMdt[i])
    end

    # for i = 1:length(Xpvapor)
    #     indexes = findall( x -> (mod(x[1],length(Xpvapor)) == mod(i,length(Xpvapor)) && (x[end] == -1)), sys.mapping.walltoliquid)
    #
    #         for j in length(indexes)
    #             dMdt[i] += Hvapor[i]*dx*(sys.wall.θarray[indexes[j]] - θ[i])
    #         end
    # end
    return dMdt

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

        xs = sys.liquid.Xarrays[i];
        dx = mod(xs[2] - xs[1], sys.tube.L)

        θ_wall_inter = sys.mapping.θ_interp_walltoliquid

        fx = map(θ_wall_inter, xs) - θarrays[i]
        du[i] = α .* laplacian(θarrays[i]) ./ dx ./ dx + Hₗ .* fx .* H_rhs
    end
    return du
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


function laplacian(u,periodic=true)
    unew = deepcopy(u)

    dl = ones(length(u)-1)
    dr = dl
    d  = -2*ones(length(u))

    A = Tridiagonal(dl, d, dr)

    unew = A*u

    #periodic B.C.
    if periodic

        unew[1]   = u[2] - 2*u[1] + u[end]

        unew[end] = u[1] - 2*u[end] + u[end-1]

    else

    # zero gradient B.C.
    unew[1]   = unew[2]
    unew[end] = unew[end-1]
    end

    return (unew)
end


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
