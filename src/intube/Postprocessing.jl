export getcurrentsys


"""
    give a new u and an old system, return a new system sysnew
"""

function getcurrentsys(u,sys0)

        indexes = Int64[]
        θliquidrec = Array[]

        for i = 1:length(u)
            if abs(u[i]+1e10) <= 10^(-1)
                push!(indexes,i)
            end
        end


    Xp,dXdt,M,δ = vectoXMδ(u[1:indexes[1]-1])
    modX!(Xp,sys0.tube.L)
    # θwallrec = u[indexes[1]+1:indexes[2]-1]

    for i = 1:length(indexes)-1
    push!(θliquidrec, u[indexes[i]+1:indexes[i+1]-1])
    end
    push!(θliquidrec, u[indexes[end]+1:end])

    sysnew = deepcopy(sys0)

    sysnew.liquid.Xp = Xp
    sysnew.liquid.dXdt = dXdt
    sysnew.liquid.θarrays = θliquidrec
    sysnew.liquid.Xarrays = updateXarrays(Xp,sysnew.liquid.θarrays,sysnew.tube.L)


    Lvaporplug = XptoLvaporplug(Xp,sys0.tube.L,sys0.tube.closedornot)
    # γ = sys0.vapor.γ

    Ac = sysnew.tube.Ac


    d = sys0.tube.d
    ρ = M ./ Lvaporplug ./ (Ac .* ((d .- 2δ) ./ d).^2 )
    P = DtoP.(ρ)
    # P = DtoP.(M ./ Lvaporplug ./ Ac)

        # println(P)
    # P = nondi_DtoP.(M./Lvaporplug)
    # P = real.((M./Lvaporplug .+ 0im).^(γ))
    sysnew.vapor.P = P
    sysnew.vapor.δ = δ

    # sysnew.wall.θarray = θwall

    θ_interp_walltoliquid, θ_interp_liquidtowall, H_interp_liquidtowall, P_interp_liquidtowall = sys_interpolation(sysnew)
    sysnew.mapping = Mapping(θ_interp_walltoliquid, θ_interp_liquidtowall, H_interp_liquidtowall, P_interp_liquidtowall)

    return sysnew
end

function modX!(Xp,L)
    for i = 1:length(Xp)
         Xp[i] = mod.(Xp[i],L)
    end

    return Xp
end

"""
    When having a new Xp because of dynamics, Xarrays need to be updated, too.
"""

function updateXarrays(Xp,θarrays,L)

    Xarrays = deepcopy(θarrays)

    for i = 1:length(Xarrays)
        if Xp[i][1] < Xp[i][2]
            Xarrays[i] = range(Xp[i][1], Xp[i][2], length=length(Xarrays[i]))
        else
            Xarrays[i] = range(Xp[i][1], Xp[i][2]+L, length=length(Xarrays[i])) .- L
            Xarrays[i] = mod.(Xarrays[i], L)
        end
    end

    # for i = 1:length(Xp)
    #     Xarrays[i] = LinRange(Xp[i][1],Xp[i][2],length(θarrays[i]))
    # end

    return Xarrays
end
