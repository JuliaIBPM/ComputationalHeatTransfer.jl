export merging_affect!,merging_condition,nucleateboiling

function merging_affect!(integrator)
    δv = 0.01

    p = deepcopy(getcurrentsys(integrator.u,sys0));
    L = p.tube.L;

    merge_flags = getmerge_flags(δv,p)
    indexmergingsite = findall(x->x == true, merge_flags)


    for i in indexmergingsite
        p = merging(p,i)
    end



    Lvaporplug = XptoLvaporplug(p.liquid.Xp,p.tube.L,p.tube.closedornot)
    M = p.vapor.P.^(1/p.vapor.γ).* Lvaporplug

    unew=[XMδtovec(p.liquid.Xp,p.liquid.dXdt,M,p.vapor.δ); wallθtovec(p.wall.θarray); liquidθtovec(p.liquid.θarrays)];

        resize!(integrator.u,size(unew,1)::Int)
    integrator.u = deepcopy(unew)
end

function merging_condition(u,t,integrator)     # only for closed loop tube

    δv = 0.01

    sys = deepcopy(getcurrentsys(integrator.u,integrator.p));

    merge_flags = getmerge_flags(δv,sys)

    return sum(merge_flags) != 0
end

function merging(sys,i)

        # get the liquid interface velocities and lengthes for merging
    Lliquidslug = XptoLliquidslug(sys.liquid.Xp,sys.tube.L)
    Lvaporplug =    XptoLvaporplug(sys.liquid.Xp,sys.tube.L,sys.tube.closedornot)


    Xpnewone = (i != 1) ? (mod((sys.liquid.Xp[i-1][1] + Lvaporplug[i]/2),L), mod((sys.liquid.Xp[i][end] - Lvaporplug[i]/2),L)) : (mod((sys.liquid.Xp[end][1] + Lvaporplug[i]/2),L), mod((sys.liquid.Xp[i][end] - Lvaporplug[i]/2),L))
    dXdtnewonevalue = (i != 1) ? (sys.liquid.dXdt[i-1][1]*Lliquidslug[i-1] + sys.liquid.dXdt[i][end]*Lliquidslug[i])/(Lliquidslug[i-1]+Lliquidslug[i]) : (sys.liquid.dXdt[end][1]*Lliquidslug[end] + sys.liquid.dXdt[i][end]*Lliquidslug[i])/(Lliquidslug[end]+Lliquidslug[i])





        systemp = deepcopy(sys)


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

    N = (i != 1) ? length([sys.liquid.Xarrays[i-1]; sys.liquid.Xarrays[i]]) : length([sys.liquid.Xarrays[end]; sys.liquid.Xarrays[i]])


    Xarraysnewone = constructXarrays([(i != 1) ? systemp.liquid.Xp[i-1] : systemp.liquid.Xp[i]],N,0.0,sys.tube.L)[1][1]
    splice!(systemp.liquid.Xarrays,i);
    (i != 1) ? splice!(systemp.liquid.Xarrays,i-1) : splice!(systemp.liquid.Xarrays,length(systemp.liquid.Xarrays));
    (i != 1) ? insert!(systemp.liquid.Xarrays,i-1,Xarraysnewone) : insert!(systemp.liquid.Xarrays,i,Xarraysnewone);

    θarraysnewone = (i != 1) ? [sys.liquid.θarrays[i-1]; sys.liquid.θarrays[i]] : [sys.liquid.θarrays[end]; sys.liquid.θarrays[i]]
    splice!(systemp.liquid.θarrays,i);
    (i != 1) ? splice!(systemp.liquid.θarrays,i-1) : splice!(systemp.liquid.θarrays,length(systemp.liquid.θarrays));
    (i != 1) ? insert!(systemp.liquid.θarrays,i-1,θarraysnewone) : insert!(systemp.liquid.θarrays,i,θarraysnewone);


    return systemp
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
