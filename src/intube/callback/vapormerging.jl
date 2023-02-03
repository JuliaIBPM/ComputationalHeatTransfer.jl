export vaporMergingAffect!,vaporMergingCondition,vaporMerging

function vaporMergingAffect!(integrator)
    println("vapor merged!")

    p = deepcopy(getcurrentsys(integrator.u,integrator.p));
    δv = 3p.tube.d

    L = p.tube.L;

    merge_flags = getVaporMergeFlags(δv,p)
    indexmergingsite = sort(findall(x->x == true, merge_flags),rev = true)

    # println(indexmergingsite)

    # println(length(p.liquid.Xarrays))


    for i in indexmergingsite
        p = vaporMerging(p,i)
    end

    Lvaporplug = XptoLvaporplug(p.liquid.Xp,p.tube.L,p.tube.closedornot)
    # M = nondi_PtoD.(p.vapor.P) .* Lvaporplug

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

function vaporMergingCondition(u,t,integrator)     # only for closed loop tube

    sys = deepcopy(getcurrentsys(integrator.u,integrator.p));
    δv = 3sys.tube.d
    merge_flags = getVaporMergeFlags(δv,sys)

    return sum(merge_flags) != 0
    # return true
end

function vaporMerging(p,i)

    Ac = p.tube.Ac
    ρₗ = p.liquid.ρ

    # get the total mass before merging
    Mvapor_old = getMvapor(p);
    Mfilm_old = sum(getMfilm(p));
    Mliquid_old = getMliquid(p);
    
    Mold = sum(Mvapor_old + Mfilm_old + Mliquid_old)


    # get the liquid interface velocities and lengthes for merging
    Lliquidslug_old = XptoLliquidslug(p.liquid.Xp,p.tube.L)
    Lvaporplug =  XptoLvaporplug(p.liquid.Xp,p.tube.L,p.tube.closedornot)


# get compensated L of merged liquid slug for mass conservation
    left_index_after = [length(Lvaporplug)-1;1:(length(Lvaporplug)-1)]
    right_index_after = [1:length(Lvaporplug)-1;1]
    right_index = [2:length(Lvaporplug);1]


    systemp = deepcopy(p)

# delete ith liquid slug
    splice!(systemp.liquid.Xp,i)
    splice!(systemp.liquid.dXdt,i)
    splice!(systemp.liquid.Xarrays,i)
    splice!(systemp.liquid.θarrays,i)


    if i != length(p.liquid.Xp)
        splice!(systemp.vapor.δstart,right_index[i])
        splice!(systemp.vapor.δend,i)
        splice!(systemp.vapor.Lfilm_start,right_index[i])
        splice!(systemp.vapor.Lfilm_end,i)    
        splice!(systemp.vapor.P,i:i+1,[(systemp.vapor.P[i]+systemp.vapor.P[right_index[i]])/2])

        else
            # splice!(systemp.vapor.δstart,right_index[i])
            splice!(systemp.vapor.δend,i)
            # splice!(systemp.vapor.Lfilm_start,right_index[i])
            splice!(systemp.vapor.Lfilm_end,i)   
            
            splice!(systemp.vapor.δstart,length(systemp.vapor.δstart));
            splice!(systemp.vapor.δstart,1);
            insert!(systemp.vapor.δstart,1,p.vapor.δstart[end])

            splice!(systemp.vapor.Lfilm_start,length(systemp.vapor.Lfilm_start));
            splice!(systemp.vapor.Lfilm_start,1);
            insert!(systemp.vapor.Lfilm_start,1,p.vapor.Lfilm_start[end])

            splice!(systemp.vapor.P,length(systemp.vapor.P));
            splice!(systemp.vapor.P,1);
            insert!(systemp.vapor.P,1,(p.vapor.P[i]+p.vapor.P[right_index[i]])/2)
        end
    
    # get the total mass after deleting liquid slugs merging
    Mvapor_temp = getMvapor(systemp);
    Mfilm_temp = sum(getMfilm(systemp));
    Mliquid_temp = getMliquid(systemp);
        
    Mtemp = sum(Mvapor_temp + Mfilm_temp + Mliquid_temp)

    Mdiff = Mold - Mtemp


    println(sum(Mdiff))
    println(sum(Mold))
    println(sum(Mtemp))

    A = (Mdiff / Ac / (ρₗ - PtoD(systemp.vapor.P[right_index_after[i]])))
    B = (Lliquidslug_old[i] + p.vapor.Lfilm_end[i] + p.vapor.Lfilm_start[right_index[i]])

    L_lengthen =  A < B ? A : B

    Xp_new = deepcopy(systemp.liquid.Xp)

    Xp_new[right_index_after[i]] = (Xp_new[right_index_after[i]][1] - L_lengthen/2,Xp_new[right_index_after[i]][2])
    Xp_new[left_index_after[i]] = (Xp_new[left_index_after[i]][1],Xp_new[left_index_after[i]][2] + L_lengthen/2)

    systemp.liquid.Xp = Xp_new

    # println(L_lengthen)
    # println(A)
    # println(Mdiff)

    return deepcopy(systemp)
end

function getVaporMergeFlags(δv,sys)

    # only for closed loop tube
    numofliquidslug = length(sys.liquid.Xp)
    numofmergingsite = numofliquidslug
    merge_flags = Array{Bool,1}(undef, numofmergingsite)

# println(numofliquidslug)

    # Xpvapor = getXpvapor(sys.liquid.Xp,sys.tube.L,sys.tube.closedornot)
    Xp = sys.liquid.Xp
    for i in 1:numofmergingsite
     # merging bubble length threshold
        merge_flags[i] = ((Xp[i][2] - Xp[i][1]) < δv && (Xp[i][2] - Xp[i][1]) >= 0) || ((Xp[i][2] - Xp[i][1] + sys.tube.L) < δv && (Xp[i][2] - Xp[i][1]) < 0) ? true : false
    end

    return merge_flags
end
