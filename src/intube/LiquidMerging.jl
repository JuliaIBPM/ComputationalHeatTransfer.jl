export merging_affect!,merging_condition,nucleateboiling,merging

function merging_affect!(integrator)


    p = deepcopy(getcurrentsys(integrator.u,integrator.p));
    δv = 2p.tube.d

    merge_flags = getmerge_flags(δv,p)
    indexmergingsite = sort(findall(x->x == true, merge_flags),rev = true)

    for i in indexmergingsite
        p = merging(p,i)
    end

    Lvaporplug = XptoLvaporplug(p.liquid.Xp,p.tube.L,p.tube.closedornot)

    Ac = p.tube.Ac

    d = p.tube.d
    δstart = p.vapor.δstart
    δend = p.vapor.δend
    Lfilm_start = p.vapor.Lfilm_start
    Lfilm_end = p.vapor.Lfilm_end

    δarea_start = Ac .* (1 .- ((d .- 2*δstart) ./ d) .^ 2);
    δarea_end = Ac .* (1 .- ((d .- 2*δend) ./ d) .^ 2);

    volume_vapor = Lvaporplug .* Ac - Lfilm_start .* δarea_start - Lfilm_end .* δarea_end
    M = PtoD.(p.vapor.P) .* volume_vapor

    unew=[XMδLtovec(p.liquid.Xp,p.liquid.dXdt,M,δstart,δend,Lfilm_start,Lfilm_end); liquidθtovec(p.liquid.θarrays)];

    resize!(integrator.u,size(unew,1)::Int)
    integrator.u = deepcopy(unew)
end

function merging_condition(u,t,integrator)     # only for closed loop tube

    p = deepcopy(getcurrentsys(integrator.u,integrator.p));
    # δv = p.tube.d > (integrator.dt*maximum(p.liquid.dXdt)[1]) ? p.tube.d : (integrator.dt*maximum(p.liquid.dXdt)[1])
    δv = 2p.tube.d

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
  
    L = p.tube.L
    δv = 2p.tube.d
    L_diff = δv
    
    Xpnewone = mod(p.liquid.Xp[left_index][1]+L_diff/2,L), mod(p.liquid.Xp[i][end] - L_diff/2,L)
    dXdtnewonevalue = (i != 1) ? (p.liquid.dXdt[i-1][1]*Lliquidslug[i-1] + p.liquid.dXdt[i][end]*Lliquidslug[i])/(Lliquidslug[i-1]+Lliquidslug[i]) : (p.liquid.dXdt[end][1]*Lliquidslug[end] + p.liquid.dXdt[i][end]*Lliquidslug[i])/(Lliquidslug[end]+Lliquidslug[i])
        #    println("hahaha")

        systemp = deepcopy(p)


    if i != 1
        splice!(systemp.liquid.Xp,i-1:i,[Xpnewone])
    else
        splice!(systemp.liquid.Xp,length(systemp.liquid.Xp));
        splice!(systemp.liquid.Xp,1);
        insert!(systemp.liquid.Xp,length(systemp.liquid.Xp)+1,Xpnewone)
    end

    if i != 1
        splice!(systemp.liquid.dXdt,i-1:i,[(dXdtnewonevalue,dXdtnewonevalue)])
    else
        splice!(systemp.liquid.dXdt,length(systemp.liquid.dXdt));
        splice!(systemp.liquid.dXdt,1);
        insert!(systemp.liquid.dXdt,length(systemp.liquid.dXdt)+1,(dXdtnewonevalue,dXdtnewonevalue))
    end

    splice!(systemp.vapor.δstart,i)
    splice!(systemp.vapor.δend,i)
    splice!(systemp.vapor.Lfilm_start,i)
    splice!(systemp.vapor.Lfilm_end,i)
    # splice!(systemp.vapor.Eratio,i)
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
    (i != 1) ? insert!(systemp.liquid.Xarrays,i-1,Xarraysnewone) : insert!(systemp.liquid.Xarrays,length(systemp.liquid.Xarrays)+1,Xarraysnewone);

    θarraysnewone = (i != 1) ? [p.liquid.θarrays[i-1][1:end-1]; (p.liquid.θarrays[i-1][end-1]+p.liquid.θarrays[i][2])/2 ;p.liquid.θarrays[i][2:end]] : [p.liquid.θarrays[end][1:end-1]; (p.liquid.θarrays[end][end-1]+p.liquid.θarrays[i][2]) / 2 ;p.liquid.θarrays[i][2:end]]
    splice!(systemp.liquid.θarrays,i);
    (i != 1) ? splice!(systemp.liquid.θarrays,i-1) : splice!(systemp.liquid.θarrays,length(systemp.liquid.θarrays));
    (i != 1) ? insert!(systemp.liquid.θarrays,i-1,θarraysnewone) : insert!(systemp.liquid.θarrays,length(systemp.liquid.θarrays)+1,θarraysnewone);

    return deepcopy(systemp)
end

function getmerge_flags(δv,sys)

    # only for closed loop tube
    numofliquidslug = length(sys.liquid.Xp)
    numofmergingsite = numofliquidslug

  

    merge_flags = Array{Bool,1}(undef, numofmergingsite)

    Xpvapor = getXpvapor(sys.liquid.Xp,sys.tube.L,sys.tube.closedornot)
    dXdt = sys.liquid.dXdt
    tstep = Main.tstep

    for i in 1:numofmergingsite
        left_index = i > 1 ? i-1 : numofliquidslug
     # merging bubble length threshold
        merge_flags[i] = mod(Xpvapor[i][2] + dXdt[i][1]*tstep - Xpvapor[i][1] - dXdt[left_index][2]*tstep, sys.tube.L) < δv || mod(Xpvapor[i][2] - Xpvapor[i][1],sys.tube.L) < δv 
    
    end

    return merge_flags
end
