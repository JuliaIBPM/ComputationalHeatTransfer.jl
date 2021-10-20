export randomXp

function randomXp(L;numofslugs=32,chargeratio=0.46)

    # L = sys.tube.L

    L_perslug=L/numofslugs*chargeratio
    L_persection=L/numofslugs

    Ls = (randn(numofslugs)/4 .+ 1).*L_perslug

    Xp1s = zeros(numofslugs);
    Xp2s = deepcopy(Xp1s);

    if minimum(Ls) > 0 && maximum(Ls) < L/numofslugs

        for i = 1:length(Xp1s)
            Xp1s[i] = (i-1)*L/numofslugs
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
