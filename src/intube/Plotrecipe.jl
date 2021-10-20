using RecipesBase

@recipe function f(val::PHPSystem;plottype="T")

    T0 = 295.0
    P0 = 220337

    # layout := (3,1)
    if plottype == "T"

            x1,y1 = stackXpTemp(val)

            legend := false
            color_palette  := :seaborn_dark
            title := "OHP temperatures"
            y1 = y1 .* T0

            return x1,y1

    elseif plottype == "ΔT"

            x2,y2 = stackXpTemp(val)

            legend := false
            color_palette  := :seaborn_dark
            title := "OHP ΔT"

            x2 = x2[1]
            y2 = y2[1]

            y2 -= nondi_PtoT(val.mapping.P_interp_liquidtowall.(x2))

            y2 = y2 .* T0

            return x2,y2


    elseif plottype == "P"

            x3,y3 = stackXpTemp(val)

            legend := false
            color_palette  := :seaborn_dark
            title := "OHP pressures"

            popfirst!(x3)
            popfirst!(y3)

            for i = 1:length(y3)
                y3[i] = val.mapping.P_interp_liquidtowall.(x3[i])
            end

            y3 = y3 .* P0

            return x3,y3
    end
end

function stackXpTemp(val::PHPSystem)
    Xpvapor = getXpvapor(val.liquid.Xp,val.tube.L,val.tube.closedornot)
    θvapor  = nondi_PtoT.(val.vapor.P)
    Xp = val.liquid.Xp

    all_θ  = []
    all_Xp = []

    push!(all_Xp,val.wall.Xarray); push!(all_θ, val.wall.θarray)

    j=1
    while j <= length(Xp)
        if Xp[j][end] >= Xp[j][1]
            push!(all_Xp,val.liquid.Xarrays[j]); push!(all_θ, val.liquid.θarrays[j])
            else
            # find the index at the end
            index = findfirst(x->x <= val.liquid.Xarrays[j][end], val.liquid.Xarrays[j])

            push!(all_Xp,val.liquid.Xarrays[j][1:index-1]); push!(all_θ, val.liquid.θarrays[j][1:index-1])
            push!(all_Xp,val.liquid.Xarrays[j][index:end]); push!(all_θ, val.liquid.θarrays[j][index:end])
        end

        j += 1
    end


    j=1
    while j <= length(Xpvapor)
        if Xpvapor[j][end] >= Xpvapor[j][1]
            push!(all_Xp,[Xpvapor[j][1],Xpvapor[j][end]]); push!(all_θ,[θvapor[j], θvapor[j]])
            else
            push!(all_Xp,[0.0,Xpvapor[j][end]]); push!(all_θ,[θvapor[j], θvapor[j]])
            push!(all_Xp,[Xpvapor[j][1],val.tube.L]); push!(all_θ,[θvapor[j], θvapor[j]])
        end

        j += 1
    end

    return all_Xp,all_θ
end
