export fixdx_affect!,fixdx_condition
# boiling_condition,
function fixdx_condition(u,t,integrator)

        p = deepcopy(getcurrentsys(integrator.u,integrator.p));

        dξ_init = p.tube.L / p.tube.N

        δξmax = 2*dξ_init
        δξmin = 0.5*dξ_init

        sys = deepcopy(getcurrentsys(integrator.u,integrator.p));

        ReconstructFlags = getReconstructFlags(δξmin,δξmax,sys)

        return sum(ReconstructFlags) != 0
        # return true
end

function getReconstructFlags(δξmin,δξmax,p)

    # only for closed loop tube
    numofliquidslug = length(p.liquid.Xp)
    L = p.tube.L
    ReconstructFlags = Array{Bool,1}(undef, numofliquidslug)

    Xarrays = p.liquid.Xarrays

    for i in 1:numofliquidslug
     # merging bubble length threshold
        ReconstructFlags[i] = mod((Xarrays[i][2] - Xarrays[i][1]),L) < δξmin || mod((Xarrays[i][2] - Xarrays[i][1]),L) > δξmax ? true : false

    end

    return ReconstructFlags
end


function fixdx_affect!(integrator)

    p = deepcopy(getcurrentsys(integrator.u,integrator.p));
    L = p.tube.L
    δξmax = 2*p.tube.d
    δξmin = 0.5*p.tube.d

    sys = deepcopy(getcurrentsys(integrator.u,integrator.p));


    ReconstructFlags = getReconstructFlags(δξmin,δξmax,p)
    indexReconstructSite = sort(findall(x->x == true, ReconstructFlags),rev = true)

    θ_interp_liquidtowall = p.mapping.θ_interp_liquidtowall

    Lliquid = XptoLliquidslug(p.liquid.Xp,sys.tube.L)
    Nliquid =  ceil.(Int, Lliquid./p.tube.d)

    # println(p.liquid.Xarrays[indexReconstructSite[1]])
    # println(indexReconstructSite)

    for i in indexReconstructSite
        println("reconstruct dx! in", i ," at ",integrator.t)

        p.liquid.Xarrays[i] = constructoneXarray((p.liquid.Xarrays[i][1],p.liquid.Xarrays[i][end]),Nliquid[i],L)
        p.liquid.θarrays[i] = θ_interp_liquidtowall(p.liquid.Xarrays[i])
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
