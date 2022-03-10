export merging_affect!,merging_condition,nucleateboiling

function merging_affect!(integrator)
    println("merged!")

    p = deepcopy(getcurrentsys(integrator.u,integrator.p));
    δv = p.tube.d > (integrator.dt*maximum(p.liquid.dXdt)[1]) ? p.tube.d : (integrator.dt*maximum(p.liquid.dXdt)[1])
    # println(δv)
    L = p.tube.L;

    merge_flags = getmerge_flags(δv,p)
    indexmergingsite = sort(findall(x->x == true, merge_flags),rev = true)

    # println(indexmergingsite)

    # println(length(p.liquid.Xarrays))
    # println(length(p.liquid.Xp))

    for i in indexmergingsite
        p = merging(p,i)

    #
    #     Lvaporplug = XptoLvaporplug(p.liquid.Xp,p.tube.L,p.tube.closedornot)
    #     M = p.vapor.P.^(1/p.vapor.γ).* Lvaporplug
    #     unew=[XMδtovec(p.liquid.Xp,p.liquid.dXdt,M,p.vapor.δ); liquidθtovec(p.liquid.θarrays)];
    # #
    # #
    #
    #     p = getcurrentsys(unew,p)
    #
    #     unew=[XMδtovec(p.liquid.Xp,p.liquid.dXdt,M,p.vapor.δ); liquidθtovec(p.liquid.θarrays)];
    end




    Lvaporplug = XptoLvaporplug(p.liquid.Xp,p.tube.L,p.tube.closedornot)
    # M = nondi_PtoD.(p.vapor.P) .* Lvaporplug

    Ac = p.tube.Ac
    δ = p.vapor.δ
    M = PtoD.(p.vapor.P) .* Lvaporplug .* Ac .* ((p.tube.d .- 2 .* δ) ./ p.tube.d) .^2
    # M = p.vapor.P.^(1/p.vapor.γ).* Lvaporplug

    unew=[XMδtovec(p.liquid.Xp,p.liquid.dXdt,M,p.vapor.δ); liquidθtovec(p.liquid.θarrays)];

    resize!(integrator.u,size(unew,1)::Int)
    integrator.u = deepcopy(unew)
end

function merging_condition(u,t,integrator)     # only for closed loop tube

    p = deepcopy(getcurrentsys(integrator.u,integrator.p));
    δv = p.tube.d > (integrator.dt*maximum(p.liquid.dXdt)[1]) ? p.tube.d : (integrator.dt*maximum(p.liquid.dXdt)[1])

    sys = deepcopy(getcurrentsys(integrator.u,integrator.p));

    merge_flags = getmerge_flags(δv,sys)

    return sum(merge_flags) != 0
    # return true
end

function merging(p,i)

        # get the liquid interface velocities and lengthes for merging
    Lliquidslug = XptoLliquidslug(p.liquid.Xp,p.tube.L)
    Lvaporplug =    XptoLvaporplug(p.liquid.Xp,p.tube.L,p.tube.closedornot)


# get compensated L of merged liquid slug for mass conservation
    left_index = i > 1 ? i-1 : length(Lvaporplug)
    right_index = i < length(Lvaporplug) ? i+1 : 1

    Mperlength_left = getMperlength(p,left_index)
    Mperlength_right = getMperlength(p,right_index)

    Mfilm = getMfilm(p)
    Mvapor = getMvapor(p)
    Mmerged = Mfilm[i]+Mvapor[i]
    # Mshortenedfilm = getMshortenedfilm(p,i)
    # Mshortenedvapor = getMshortenedvapor(p,i)

    #

    # println(Mshortenedfilm)
    # println(Mshortenedvapor)
    # println(i)
    # println(Mvapor[i])
    #
    # Lliquid_add = (Mfilm[i]+Mvapor[i])/ρₗ/Ac

    Lliquid_adjust = (Mmerged - Mperlength_left*Lvaporplug[i]/2 - Mperlength_right*Lvaporplug[i]/2) / (ρₗ*Ac - Mperlength_left/2 - Mperlength_right/2)
    # println(Lliquid_modify)
    # println(p.liquid.Xp)
    # println("Lliquid_adjust",Lliquid_adjust)
    # println(Mmerged)
    # println(ρₗ*Ac*Lliquid_adjust +(Mperlength_left+Mperlength_right)*(Lvaporplug[i]/2 - Lliquid_adjust/2) )
    # println((Mperlength_left+Mperlength_right)*(Lvaporplug[i]/2 - Lliquid_adjust/2))
    # println((Mperlength_left+Mperlength_right)*(Lvaporplug[i]/2))
    # println(ρₗ*Ac)
    # println((Mperlength_left+Mperlength_right)*(0.5))

    Xpnewone = (i != 1) ? (mod((p.liquid.Xp[i-1][1] + Lvaporplug[i]/2 - Lliquid_adjust/2),L), mod((p.liquid.Xp[i][end] - Lvaporplug[i]/2 + Lliquid_adjust/2),L)) : (mod((p.liquid.Xp[end][1] + Lvaporplug[i]/2 - Lliquid_adjust/2),L), mod((p.liquid.Xp[i][end] - Lvaporplug[i]/2 + Lliquid_adjust/2),L))
    # Xpnewone = (i != 1) ? (mod((p.liquid.Xp[i-1][1] + Lvaporplug[i]/2),L), mod((p.liquid.Xp[i][end] - Lvaporplug[i]/2),L)) : (mod((p.liquid.Xp[end][1] + Lvaporplug[i]/2),L), mod((p.liquid.Xp[i][end] - Lvaporplug[i]/2),L))
    dXdtnewonevalue = (i != 1) ? (p.liquid.dXdt[i-1][1]*Lliquidslug[i-1] + p.liquid.dXdt[i][end]*Lliquidslug[i])/(Lliquidslug[i-1]+Lliquidslug[i]) : (p.liquid.dXdt[end][1]*Lliquidslug[end] + p.liquid.dXdt[i][end]*Lliquidslug[i])/(Lliquidslug[end]+Lliquidslug[i])


    # println(Xpnewone)


        systemp = deepcopy(p)


    if i != 1
        splice!(systemp.liquid.Xp,i-1:i,[Xpnewone])
    else
        splice!(systemp.liquid.Xp,length(systemp.liquid.Xp));
        splice!(systemp.liquid.Xp,1);
        insert!(systemp.liquid.Xp,1,Xpnewone)
    end

    if i != 1
        splice!(systemp.liquid.dXdt,i-1:i,[(dXdtnewonevalue,dXdtnewonevalue)])
    else
        splice!(systemp.liquid.dXdt,length(systemp.liquid.dXdt));
        splice!(systemp.liquid.dXdt,1);
        insert!(systemp.liquid.dXdt,1,(dXdtnewonevalue,dXdtnewonevalue))
    end

    splice!(systemp.vapor.δ,i)
    splice!(systemp.vapor.P,i)

    # N = (i != 1) ? length([sys.liquid.Xarrays[i-1]; sys.liquid.Xarrays[i]]) : length([sys.liquid.Xarrays[end]; sys.liquid.Xarrays[i]])

    # N = length(p.wall.Xarray)


    Nliquids = (i != 1) ? length([systemp.liquid.Xarrays[i-1]; systemp.liquid.Xarrays[i]]) : length([systemp.liquid.Xarrays[end]; systemp.liquid.Xarrays[i]])
    # Nliquid = Nliquid - 1

    Xarraysnewone = constructoneXarray((i != 1) ? systemp.liquid.Xp[i-1] : systemp.liquid.Xp[end],Nliquids-1,p.tube.L)

    # println(length(Xarraysnewone))
    # println(Nliquid)
    # println(length(systemp.liquid.Xarrays[i]))
    # println( (i != 1) ? length(systemp.liquid.Xarrays[i-1]) : length(systemp.liquid.Xarrays[end]))

    splice!(systemp.liquid.Xarrays,i);
    (i != 1) ? splice!(systemp.liquid.Xarrays,i-1) : splice!(systemp.liquid.Xarrays,length(systemp.liquid.Xarrays));
    (i != 1) ? insert!(systemp.liquid.Xarrays,i-1,Xarraysnewone) : insert!(systemp.liquid.Xarrays,i,Xarraysnewone);

    θarraysnewone = (i != 1) ? [p.liquid.θarrays[i-1][1:end-1]; (p.liquid.θarrays[i-1][end-1]+p.liquid.θarrays[i][2])/2 ;p.liquid.θarrays[i][2:end]] : [p.liquid.θarrays[end][1:end-1]; (p.liquid.θarrays[end][end-1]+p.liquid.θarrays[i][2]) / 2 ;p.liquid.θarrays[i][2:end]]
    splice!(systemp.liquid.θarrays,i);
    (i != 1) ? splice!(systemp.liquid.θarrays,i-1) : splice!(systemp.liquid.θarrays,length(systemp.liquid.θarrays));
    (i != 1) ? insert!(systemp.liquid.θarrays,i-1,θarraysnewone) : insert!(systemp.liquid.θarrays,i,θarraysnewone);

    # println(sum(getMfilm(p)))
    # println(sum(getMfilm(systemp)))
    # println(sum(getMvapor(p)))
    # println(sum(getMvapor(systemp)))
    # println(sum(getMfilm(p))+sum(getMvapor(p)))
    # println(sum(getMfilm(systemp))+sum(getMvapor(systemp)))
    # println(getMliquid(p))
    # println(getMliquid(systemp))
    # println(p.vapor.P)
    # println(systemp.vapor.P)
    # println(XptoLvaporplug(p.liquid.Xp,p.tube.L,p.tube.closedornot))
    # println(XptoLvaporplug(systemp.liquid.Xp,systemp.tube.L,systemp.tube.closedornot))

    return deepcopy(systemp)
end

function getmerge_flags(δv,sys)

    # only for closed loop tube
    numofliquidslug = length(sys.liquid.Xp)
    numofmergingsite = numofliquidslug
    merge_flags = Array{Bool,1}(undef, numofmergingsite)

    Xpvapor = getXpvapor(sys.liquid.Xp,sys.tube.L,sys.tube.closedornot)

    for i in 1:numofmergingsite
     # merging bubble length threshold
        merge_flags[i] = ((Xpvapor[i][2] - Xpvapor[i][1]) < δv && (Xpvapor[i][2] - Xpvapor[i][1]) >= 0) || ((Xpvapor[i][2] - Xpvapor[i][1] + sys.tube.L) < δv && (Xpvapor[i][2] - Xpvapor[i][1]) < 0) ? true : false

    end

    return merge_flags
end
