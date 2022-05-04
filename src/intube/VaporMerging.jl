export vaporMergingAffect!,vaporMergingCondition

function vaporMergingAffect!(integrator)
    println("vapor merged!")

    p = deepcopy(getcurrentsys(integrator.u,integrator.p));
    δv = p.tube.d

    L = p.tube.L;

    merge_flags = getVaporMergeFlags(δv,p)
    indexmergingsite = sort(findall(x->x == true, merge_flags),rev = true)

    for i in indexmergingsite
        p = vaporMerging(p,i)
    end

    # println(length(p.liquid.Xarrays))

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

function vaporMergingCondition(u,t,integrator)     # only for closed loop tube

    sys = deepcopy(getcurrentsys(integrator.u,integrator.p));
    δv = sys.tube.d
    merge_flags = getVaporMergeFlags(δv,sys)

    return sum(merge_flags) != 0
    # return true
end

function vaporMerging(p,i)

        # get the liquid interface velocities and lengthes for merging
    Lliquidslug = XptoLliquidslug(p.liquid.Xp,p.tube.L)
    Lvaporplug =  XptoLvaporplug(p.liquid.Xp,p.tube.L,p.tube.closedornot)


# get compensated L of merged liquid slug for mass conservation
    left_index = i > 1 ? i-1 : length(Lvaporplug)
    right_index = i < length(Lvaporplug) ? i+1 : 1

    Mperlength_left = getMperlength(p,left_index)
    Mperlength_right = getMperlength(p,right_index)

    Mfilm = getMfilm(p)
    Mvapor = getMvapor(p)
    Mliquid = getMliquid(p)
    # Mmerged = Mfilm[i]+Mvapor[i]

    MvaporPreMerging = Mvapor[i] + Mvapor[right_index]
    MfilmPreMerging = Mfilm[i] + Mfilm[right_index]

    MvaporNew = MvaporPreMerging
    MfilmNew = MfilmPreMerging + Mliquid[i]

    LvaporNew = Lvaporplug[i] + Lliquidslug[i] + Lvaporplug[right_index]
    δareaNew = MfilmNew ./ LvaporNew ./ p.liquid.ρ
    δNewOne = getδFromδarea(p.tube.Ac,p.tube.d,δareaNew)

    PNewOne = DtoP(MvaporNew/LvaporNew/(Ac-δareaNew))

    systemp = deepcopy(p)

# delete ith liquid slug
    splice!(systemp.liquid.Xp,i)
    splice!(systemp.liquid.dXdt,i)
    splice!(systemp.liquid.Xarrays,i)
    splice!(systemp.liquid.θarrays,i)

    # splice!(systemp.vapor.δ,i)
    # splice!(systemp.vapor.P,i)

    if i != length(systemp.liquid.Xp)
        splice!(systemp.vapor.δ,i:i+1,[δNewOne])
        splice!(systemp.vapor.P,i:i+1,[PNewOne])
    else
        splice!(systemp.liquid.δ,length(systemp.liquid.δ));
        splice!(systemp.liquid.δ,1);
        insert!(systemp.liquid.δ,1,δNewOne)

        splice!(systemp.liquid.P,length(systemp.liquid.P));
        splice!(systemp.liquid.P,1);
        insert!(systemp.liquid.P,1,PNewOne)
    end

    return deepcopy(systemp)
end

function getVaporMergeFlags(δv,sys)

    # only for closed loop tube
    numofliquidslug = length(sys.liquid.Xp)
    numofmergingsite = numofliquidslug
    merge_flags = Array{Bool,1}(undef, numofmergingsite)

    # Xpvapor = getXpvapor(sys.liquid.Xp,sys.tube.L,sys.tube.closedornot)
    Xp = sys.liquid.Xp
    for i in 1:numofmergingsite
     # merging bubble length threshold
        merge_flags[i] = ((Xp[i][2] - Xp[i][1]) < δv && (Xp[i][2] - Xp[i][1]) >= 0) || ((Xp[i][2] - Xp[i][1] + sys.tube.L) < δv && (Xp[i][2] - Xp[i][1]) < 0) ? true : false
    end

    return merge_flags
end
