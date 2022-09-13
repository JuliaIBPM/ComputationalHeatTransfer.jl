export constructmapping,sys_interpolation

function sys_interpolation(sys)
    X_inner = Array{Float64}(undef, 0)
    θ_inner = Array{Float64}(undef, 0)
    H_inner = Array{Float64}(undef, 0)
    X_inner_pres  = Array{Float64}(undef, 0)
    P_inner = Array{Float64}(undef, 0)

    Xp  = sys.liquid.Xp

    θ = PtoT.(sys.vapor.P)
    P = sys.vapor.P
    δstart = sys.vapor.δstart
    δend = sys.vapor.δend
    Lfilm_start = sys.vapor.Lfilm_start
    Lfilm_end = sys.vapor.Lfilm_end
    Xpvapor = getXpvapor(Xp,sys.tube.L,sys.tube.closedornot)

    H_film_start = Hfilm.(δstart,[sys])
    H_film_end = Hfilm.(δend,[sys])
    H_vapor = sys.vapor.Hᵥ
    H_liquid = sys.liquid.Hₗ

    Nvapor = length(P)

    loop_plus_index = [2:Nvapor;1]

    max_i = argmax(Xpvapor)
    max_j = argmax(Xpvapor[max_i])
  

    for i = 1:Nvapor

            if i == max_i && max_j != 2

                X_inner_loop_append,H_inner_loop_append = XHloop_append(i,H_film_start[i],H_film_end[i],sys)

                append!(X_inner,X_inner_loop_append)
                append!(θ_inner,[θ[i], θ[i], θ[i], θ[i], θ[i], θ[i], θ[i], θ[i]])
                append!(H_inner,H_inner_loop_append)

                append!(X_inner_pres,[Xpvapor[i][1], sys.tube.L, 0.0, Xpvapor[i][end]])
                append!(P_inner,[P[i], P[i], P[i], P[i]])
            else
                append!(X_inner,[Xpvapor[i][1],Xpvapor[i][1]+Lfilm_start[i],
                Xpvapor[i][1]+Lfilm_start[i],Xpvapor[i][end]-Lfilm_end[i],
                Xpvapor[i][end]-Lfilm_end[i],Xpvapor[i][end]])

                append!(θ_inner,[θ[i],θ[i],θ[i],θ[i],θ[i],θ[i]])
                append!(H_inner,[H_film_start[i], H_film_start[i],
                H_vapor,H_vapor,
                H_film_end[i],H_film_end[i]])

                append!(X_inner_pres,[Xpvapor[i][1],Xpvapor[i][end]])
                append!(P_inner,[P[i], P[i]])
    
        end
        

    if (sys.liquid.Xarrays[i][1] > sys.liquid.Xarrays[i][end])
        period_index = findmin(sys.liquid.Xarrays[i])[2]
        Xarrays_temp = deepcopy(sys.liquid.Xarrays[i])
        θarrays_temp = deepcopy(sys.liquid.θarrays[i])
        
        
        insert!(Xarrays_temp,period_index, sys.tube.L)
        insert!(Xarrays_temp,period_index+1, 0.0)
        insert!(θarrays_temp,period_index, sys.liquid.θarrays[i][period_index])
        insert!(θarrays_temp,period_index+1, sys.liquid.θarrays[i][period_index])
        append!(X_inner,Xarrays_temp)
        append!(θ_inner,θarrays_temp)
        
        append!(H_inner,H_liquid .* ones(length(sys.liquid.Xarrays[i]) + 2))
        

        append!(X_inner_pres,[sys.liquid.Xarrays[i][1], sys.tube.L, 0.0, sys.liquid.Xarrays[i][end]])
        P_inner_end = (sys.tube.L-sys.liquid.Xarrays[i][1])/mod(sys.liquid.Xarrays[i][end]-sys.liquid.Xarrays[i][1],sys.tube.L) * (P[loop_plus_index[i]]-P[i]) + P[i]
        append!(P_inner,[P[i], P_inner_end, P_inner_end, P[loop_plus_index[i]]])

        else
            append!(X_inner,sys.liquid.Xarrays[i])
            append!(θ_inner,sys.liquid.θarrays[i])
            append!(H_inner,H_liquid .* ones(length(sys.liquid.Xarrays[i])))

    end
end

        imin = argmin(X_inner)
        X_inner_final =Array{Float64}(undef, 0)
        θ_inner_final =Array{Float64}(undef, 0)
        H_inner_final =Array{Float64}(undef, 0)
        if imin != 1
            X_inner_final = [X_inner[imin:end];X_inner[1:imin-1]]
            θ_inner_final = [θ_inner[imin:end];θ_inner[1:imin-1]]
            H_inner_final = [H_inner[imin:end];H_inner[1:imin-1]]
        end

        imin = argmin(X_inner_pres)
        X_inner_pres_final =Array{Float64}(undef, 0)
        P_inner_final =Array{Float64}(undef, 0)
        if imin != 1
            X_inner_pres_final = [X_inner_pres[imin:end];X_inner_pres[1:imin-1]]
            P_inner_final = [P_inner[imin:end];P_inner[1:imin-1]]
        end
    
    ``` extend wall Xarray by adding its 0.0 point```
    extend_wall_Xarray = deepcopy(sys.wall.Xarray)
    extend_wall_θarray = deepcopy(sys.wall.θarray)

    extend_wall_Xarray = [0.0;sys.wall.Xarray;sys.tube.L]
    extend_wall_θarray = [(sys.wall.θarray[1]+sys.wall.θarray[end])/2;sys.wall.θarray;(sys.wall.θarray[1]+sys.wall.θarray[end])/2]

    Interpolations.deduplicate_knots!(X_inner_final,move_knots = true)
    Interpolations.deduplicate_knots!(extend_wall_Xarray,move_knots = true)
    Interpolations.deduplicate_knots!(X_inner_pres_final,move_knots = true)

    θ_interp_liquidtowall = LinearInterpolation(X_inner_final, θ_inner_final);

    H_interp_liquidtowall = LinearInterpolation(X_inner_final, H_inner_final);

    θ_interp_walltoliquid = LinearInterpolation(extend_wall_Xarray, extend_wall_θarray);

    P_interp_liquidtowall = LinearInterpolation(X_inner_pres_final, P_inner_final);


    return θ_interp_walltoliquid, θ_interp_liquidtowall, H_interp_liquidtowall, P_interp_liquidtowall

end

function XHloop_append(i,H_film_start,H_film_end,sys)


    Xp = sys.liquid.Xp[i]
    Xpvapor = getXpvapor(sys.liquid.Xp,sys.tube.L,sys.tube.closedornot)[i]
    Lfilm_start = sys.vapor.Lfilm_start[i]
    Lfilm_end = sys.vapor.Lfilm_end[i]

    H_vapor = sys.vapor.Hᵥ

    L = sys.tube.L

    if Xpvapor[2] > Xpvapor[1]
        println("Xp error 1")
        return "Xp error 1"
        
    elseif mod(Xpvapor[1],L) + Lfilm_start > L
        
        X_inner_loop_append = mod.([Xpvapor[1],L,
                0.0,Xpvapor[1]+Lfilm_start,
                Xpvapor[1]+Lfilm_start,Xpvapor[end]-Lfilm_end,
                Xpvapor[end]-Lfilm_end,Xpvapor[end]],L)
        X_inner_loop_append[2] = L
        
        H_inner_loop_append = [H_film_start, H_film_start,
                H_film_start, H_film_start,
                H_vapor,H_vapor,
                H_film_end,H_film_end]
        
    elseif Xpvapor[2] - Lfilm_end > 0.0
        
        X_inner_loop_append = mod.([Xpvapor[1],Xpvapor[1]+Lfilm_start,
                Xpvapor[1]+Lfilm_start, L,
                0.0, Xpvapor[end]-Lfilm_end,
                Xpvapor[end]-Lfilm_end,Xpvapor[end]],L)
        X_inner_loop_append[4] = L
        
        H_inner_loop_append = [H_film_start, H_film_start,
                H_vapor, H_vapor,
                H_vapor, H_vapor,
                H_film_end,H_film_end]
        
    elseif Xpvapor[2] > 0
        
        X_inner_loop_append = mod.([Xpvapor[1],Xpvapor[1]+Lfilm_start,
                Xpvapor[1]+Lfilm_start,  Xpvapor[end]-Lfilm_end,
                Xpvapor[end]-Lfilm_end,L,
                0.0,Xpvapor[end]],L)
        X_inner_loop_append[6] = L
        
        H_inner_loop_append = [H_film_start, H_film_start,
                H_vapor, H_vapor,
                H_film_end,H_film_end,
                H_film_end,H_film_end]
        
    else 
        println("Xp error 1")
        return "Xp error 2"
    end

return X_inner_loop_append,H_inner_loop_append
end