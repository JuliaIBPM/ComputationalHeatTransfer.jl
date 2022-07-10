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


    Xp,dXdt,M,δstart,δend,Lfilm_start,Lfilm_end = vectoXMδL(u[1:indexes[1]-1])

    # println(u[1:300])
    modX!(Xp,sys0.tube.L)
   
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

    Ac = sysnew.tube.Ac


    d = sys0.tube.d
    δarea_start = Ac .* (1 .- ((d .- 2*δstart) ./ d) .^ 2);
    δarea_end = Ac .* (1 .- ((d .- 2*δend) ./ d) .^ 2);

    volume_vapor = Lvaporplug .* Ac - Lfilm_start .* δarea_start - Lfilm_end .* δarea_end
    ρ = M ./ volume_vapor
    P = DtoP.(ρ)
  
    sysnew.vapor.P = P
    sysnew.vapor.δstart = δstart
    sysnew.vapor.δend = δend
    sysnew.vapor.Lfilm_start = Lfilm_start
    sysnew.vapor.Lfilm_end = Lfilm_end

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

    Xarrays = zero.(deepcopy(θarrays))

    for i = 1:length(Xarrays)
        if Xp[i][1] < Xp[i][2]
            Xarrays[i] = range(Xp[i][1], Xp[i][2], length=length(Xarrays[i]))
        else
            Xarrays[i] = mod.(range(Xp[i][1], Xp[i][2]+L, length=length(Xarrays[i])),L)
        end
        Xarrays[i][1] = Xp[i][1]
        Xarrays[i][end] = Xp[i][end]
    end

    return Xarrays
end
