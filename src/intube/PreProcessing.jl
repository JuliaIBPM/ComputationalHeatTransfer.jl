export randomXp

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
