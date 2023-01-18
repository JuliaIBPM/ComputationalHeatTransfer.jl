export randomXp,L_to_boiltime

function randomXp(tube::Tube,numofslugs=30,chargeratio=0.46,σ_charge=0.1)

    L = tube.L
    Lmin = tube.d

    σ_persection = σ_charge*L/sqrt(numofslugs)

    L_perslug=L/numofslugs*chargeratio
    L_persection=L/numofslugs

    Ls = abs.((rand(numofslugs) .- 0.5).*σ_persection .+ L_perslug)

    Xp1s = zeros(numofslugs);
    Xp2s = deepcopy(Xp1s);

    if minimum(Ls) > Lmin && maximum(Ls) < L_persection

        for i in eachindex(Xp1s)
            Xp1s[i] = (i-1)*L_persection
            Xp2s[i] = Xp1s[i] + Ls[i]
        end

        displacement = L*rand()

        Xp1s = mod.(Xp1s.+displacement,L)
        Xp2s = mod.(Xp2s.+displacement,L)

        else println("generation failed")
    end

    X0 = map(tuple,Xp1s,Xp2s)
    real_ratio = sum(Ls)/L

    X0,real_ratio
end

function L_to_boiltime(L_new_bubble,Rn,fluid_type,vapor::Vapor,tube::Tube)
    P = vapor.P[1]
    property = SaturationFluidProperty(fluid_type,PtoT(P));

    tube_d = tube.d
    ρₗ = property.ρₗ
    Cpₗ = property.Cpₗ 
    kₗ = property.kₗ

    Hfg = PtoHfg(P)
    Tref = PtoT(P)
    ρv = PtoD(P)
    ΔTthres = RntoΔT(Rn,Tref,fluid_type,tube_d)

    A = (2 * (ΔTthres/Tref) * Hfg*ρv/ρₗ)^0.5
    Ja = ΔTthres*Cpₗ*ρₗ/ρv/Hfg
    B = (12*kₗ/pi/ρₗ/Cpₗ)^0.5 * Ja

    t = 1e-3:1e-3:10e0
    tstar = t .* A^2 ./ B^2 

    Rplus = 2/3 .* ((tstar .+ 1).^1.5 .- (tstar).^1.5 .- 1);
    R = Rplus .* B^2 ./ A; 
    L_equivalent = (4/3) ./ tube_d .* R.^2 

    interp_Rtot = LinearInterpolation(L_equivalent, t);

    return interp_Rtot(L_new_bubble)
end
