export randomXp

function randomXp(tube;numofslugs=32,chargeratio=0.46,σ_charge=0.01)

    L = tube.L
    Lmin = tube.d

    σ_persection = σ_charge*L/sqrt(numofslugs)

    L_perslug=L/numofslugs*chargeratio
    L_persection=L/numofslugs

    # println(L_persection,σ_persection)
    # if chargeratio <= 0.5
        Ls = abs.(randn(numofslugs).*σ_persection .+ L_perslug)
    # else
    #     Ls = ((2*chargeratio-1) .+ rand(numofslugs).*(2-2*chargeratio)).*L_persection
    # end

    # println(Ls)

    Xp1s = zeros(numofslugs);
    Xp2s = deepcopy(Xp1s);

    if minimum(Ls) > Lmin && maximum(Ls) < L_persection

        for i = 1:length(Xp1s)
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
