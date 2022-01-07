export constructmapping,sys_interpolation

function sys_interpolation(sys)
    X_inner = Array{Float64}(undef, 0)
    θ_inner = Array{Float64}(undef, 0)
    H_inner = Array{Float64}(undef, 0)
    X_inner_pres  = Array{Float64}(undef, 0)
    P_inner = Array{Float64}(undef, 0)
    mk = 0
    θ = PtoT.(sys.vapor.P)
    # θ = nondi_PtoT.(sys.vapor.P)
    P = sys.vapor.P
    # θ = real.((sys.vapor.P .+ 0im).^((sys.vapor.γ-1)/sys.vapor.γ)) # isentropic
    H_vapor = sys.vapor.k ./ sys.vapor.δ
    H_liquid = sys.liquid.Hₗ


    for i = 1:length(sys.liquid.Xarrays)

            if i != 1 && (sys.liquid.Xarrays[i][1] < sys.liquid.Xarrays[i-1][end])
                append!(X_inner,[sys.liquid.Xarrays[i-1][end], sys.tube.L, 0.0, sys.liquid.Xarrays[i][1]])
                append!(θ_inner,[θ[i], θ[i], θ[i], θ[i]])
                append!(H_inner,[H_vapor[i], H_vapor[i], H_vapor[i], H_vapor[i]])

                append!(X_inner_pres,[sys.liquid.Xarrays[i-1][end], sys.tube.L, 0.0, sys.liquid.Xarrays[i][1]])
                append!(P_inner,[P[i], P[i], P[i], P[i]])
                mk = i
            elseif i == 1 && (sys.liquid.Xarrays[i][1] < sys.liquid.Xarrays[end][end])
                append!(X_inner,[sys.liquid.Xarrays[end][end], sys.tube.L, 0.0, sys.liquid.Xarrays[i][1]])
                append!(θ_inner,[θ[i], θ[i], θ[i], θ[i]])
                append!(H_inner,[H_vapor[i], H_vapor[i], H_vapor[i], H_vapor[i]])

                append!(X_inner_pres,[sys.liquid.Xarrays[end][end], sys.tube.L, 0.0, sys.liquid.Xarrays[i][1]])
                append!(P_inner,[P[i], P[i], P[i], P[i]])
                mk = i
            elseif i != 1
                append!(X_inner,[sys.liquid.Xarrays[i-1][end],sys.liquid.Xarrays[i][1]])
                append!(θ_inner,[θ[i],θ[i]])
                append!(H_inner,[H_vapor[i], H_vapor[i]])

                append!(X_inner_pres,[sys.liquid.Xarrays[i-1][end],sys.liquid.Xarrays[i][1]])
                append!(P_inner,[P[i], P[i]])
            elseif i == 1
                append!(X_inner,[sys.liquid.Xarrays[end][end],sys.liquid.Xarrays[i][1]])
                append!(θ_inner,[θ[i],θ[i]])
                append!(H_inner,[H_vapor[i], H_vapor[i]])

                append!(X_inner_pres,[sys.liquid.Xarrays[end][end],sys.liquid.Xarrays[i][1]])
                append!(P_inner,[P[i], P[i]])
        end

# println(X_inner)
# println(sys.liquid.Xarrays)
# println(sys.liquid.Xp)

            # Harrays_temp = deepcopy(sys.liquid.Xarrays[i])

            if (sys.liquid.Xarrays[i][1] > sys.liquid.Xarrays[i][end])
                period_index = findmin(sys.liquid.Xarrays[i])[2]
                Xarrays_temp = deepcopy(sys.liquid.Xarrays[i])
                θarrays_temp = deepcopy(sys.liquid.θarrays[i])

                # println(sys.liquid.Xarrays[i][1],"\n",sys.liquid.Xarrays[i][end],"\n",length(sys.liquid.θarrays[i]))

                insert!(Xarrays_temp,period_index, sys.tube.L)
                insert!(Xarrays_temp,period_index+1, 0.0)
                insert!(θarrays_temp,period_index, sys.liquid.θarrays[i][period_index])
                insert!(θarrays_temp,period_index+1, sys.liquid.θarrays[i][period_index])
                append!(X_inner,Xarrays_temp)
                append!(θ_inner,θarrays_temp)


                append!(H_inner,H_liquid .* ones(length(sys.liquid.Xarrays[i]) + 2))

                if i != length(sys.liquid.Xarrays)
                    append!(X_inner_pres,[sys.liquid.Xarrays[i][1], sys.tube.L, 0.0, sys.liquid.Xarrays[i][end]])
                    P_inner_end = (sys.tube.L-sys.liquid.Xarrays[i][1])/mod(sys.liquid.Xarrays[i][end]-sys.liquid.Xarrays[i][1],sys.tube.L) * (P[i+1]-P[i]) + P[i]
                    append!(P_inner,[P[i], P_inner_end, P_inner_end, P[i+1]])
                elseif i == length(sys.liquid.Xarrays)
                    append!(X_inner_pres,[sys.liquid.Xarrays[i][1], sys.tube.L, 0.0, sys.liquid.Xarrays[i][end]])
                    P_inner_end = (sys.tube.L-sys.liquid.Xarrays[i][1])/mod(sys.liquid.Xarrays[i][end]-sys.liquid.Xarrays[i][1],sys.tube.L) * (P[1]-P[i]) + P[i]
                    append!(P_inner,[P[i], P_inner_end, P_inner_end, P[1]])
                end

            else
                append!(X_inner,sys.liquid.Xarrays[i])
                append!(θ_inner,sys.liquid.θarrays[i])
                append!(H_inner,H_liquid .* ones(length(sys.liquid.Xarrays[i])))
            end

    end



        # println(length(X_inner))
        # println(length(θ_inner))
    # if mk != 0
    #         if mk != 1
    #             append!(X_inner,[sys.liquid.Xarrays[i-1][end], sys.tube.L])
    #             append!(T_inner,[θ[mk],θ[mk]])
    #         elseif mk == 1
    #             append!(X_inner,[sys.liquid.Xarrays[end][end], sys.tube.L])
    #             append!(T_inner,[θ[mk],θ[mk]])
    #         end
    # end


    imin = argmin(X_inner)
    X_inner_final =Array{Float64}(undef, 0)
    θ_inner_final =Array{Float64}(undef, 0)
    H_inner_final =Array{Float64}(undef, 0)
    if imin != 1
        append!(X_inner_final,view(X_inner, imin:length(X_inner)))
        append!(X_inner_final,view(X_inner, 1:imin-1))
        append!(θ_inner_final,view(θ_inner, imin:length(θ_inner)))
        append!(θ_inner_final,view(θ_inner, 1:imin-1))
        append!(H_inner_final,view(H_inner, imin:length(H_inner)))
        append!(H_inner_final,view(H_inner, 1:imin-1))
    end

    imin = argmin(X_inner_pres)
    X_inner_pres_final =Array{Float64}(undef, 0)
    P_inner_final =Array{Float64}(undef, 0)
    if imin != 1
        append!(X_inner_pres_final,view(X_inner_pres, imin:length(X_inner_pres)))
        append!(X_inner_pres_final,view(X_inner_pres, 1:imin-1))
        append!(P_inner_final,view(P_inner, imin:length(P_inner)))
        append!(P_inner_final,view(P_inner, 1:imin-1))
    end


``` extend wall Xarray by adding its 0.0 point```
extend_wall_Xarray = deepcopy(sys.wall.Xarray)
extend_wall_θarray = deepcopy(sys.wall.θarray)

prepend!(extend_wall_Xarray,[0.0])
prepend!(extend_wall_θarray,[sys.wall.θarray[end]])

append!(extend_wall_Xarray,[sys.tube.L])
append!(extend_wall_θarray,[sys.wall.θarray[end]])

# println(length(extend_wall_Xarray))
# println(length(extend_wall_θarray))
# println(minimum(θ_inner_final))
# println(θ_inner_final[1])

# test
# X_test = zeros(length(X_inner_final)-1)
# println(length(X_test))
# for i=1:length(X_test)
#     X_test[i] = X_inner_final[i+1] - X_inner_final[i]
# end
#     println(findmin(X_test))
# println(X_inner_final[end-5:end])
    # println(length(X_test))
    # println(X_inner_final[4195:4210])
    # println(θ_inner_final[4195:4210])
# return X_inner_final
    θ_interp_liquidtowall = LinearInterpolation(X_inner_final, θ_inner_final);

    H_interp_liquidtowall = LinearInterpolation(X_inner_final, H_inner_final);

    θ_interp_walltoliquid = LinearInterpolation(extend_wall_Xarray, extend_wall_θarray);

    P_interp_liquidtowall = LinearInterpolation(X_inner_pres_final, P_inner_final);

    return θ_interp_walltoliquid, θ_interp_liquidtowall, H_interp_liquidtowall, P_interp_liquidtowall
end


function test_interpolation(sys)
    X_inner = Array{Float64}(undef, 0)
    θ_inner = Array{Float64}(undef, 0)
    H_inner = Array{Float64}(undef, 0)
    X_inner_pres  = Array{Float64}(undef, 0)
    P_inner = Array{Float64}(undef, 0)
    mk = 0
    # θ = nondi_PtoT.(sys.vapor.P)
    θ = PtoT.(sys.vapor.P)
    P = sys.vapor.P
    # θ = real.((sys.vapor.P .+ 0im).^((sys.vapor.γ-1)/sys.vapor.γ)) # isentropic
    H_vapor = sys.vapor.k ./ sys.vapor.δ
    H_liquid = sys.liquid.Hₗ



    for i = 1:length(sys.liquid.Xarrays)

            if i != 1 && (sys.liquid.Xarrays[i][1] < sys.liquid.Xarrays[i-1][end])
                append!(X_inner,[sys.liquid.Xarrays[i-1][end], sys.tube.L, 0.0, sys.liquid.Xarrays[i][1]])
                append!(θ_inner,[θ[i], θ[i], θ[i], θ[i]])
                append!(H_inner,[H_vapor[i], H_vapor[i], H_vapor[i], H_vapor[i]])

                append!(X_inner_pres,[sys.liquid.Xarrays[i-1][end], sys.tube.L, 0.0, sys.liquid.Xarrays[i][1]])
                append!(P_inner,[P[i], P[i], P[i], P[i]])
                mk = i
            elseif i == 1 && (sys.liquid.Xarrays[i][1] < sys.liquid.Xarrays[end][end])
                append!(X_inner,[sys.liquid.Xarrays[end][end], sys.tube.L, 0.0, sys.liquid.Xarrays[i][1]])
                append!(θ_inner,[θ[i], θ[i], θ[i], θ[i]])
                append!(H_inner,[H_vapor[i], H_vapor[i], H_vapor[i], H_vapor[i]])

                append!(X_inner_pres,[sys.liquid.Xarrays[end][end], sys.tube.L, 0.0, sys.liquid.Xarrays[i][1]])
                append!(P_inner,[P[i], P[i], P[i], P[i]])
                mk = i
            elseif i != 1
                append!(X_inner,[sys.liquid.Xarrays[i-1][end],sys.liquid.Xarrays[i][1]])
                append!(θ_inner,[θ[i],θ[i]])
                append!(H_inner,[H_vapor[i], H_vapor[i]])

                append!(X_inner_pres,[sys.liquid.Xarrays[i-1][end],sys.liquid.Xarrays[i][1]])
                append!(P_inner,[P[i], P[i]])
            elseif i == 1
                append!(X_inner,[sys.liquid.Xarrays[end][end],sys.liquid.Xarrays[i][1]])
                append!(θ_inner,[θ[i],θ[i]])
                append!(H_inner,[H_vapor[i], H_vapor[i]])

                append!(X_inner_pres,[sys.liquid.Xarrays[end][end],sys.liquid.Xarrays[i][1]])
                append!(P_inner,[P[i], P[i]])
        end



            # Harrays_temp = deepcopy(sys.liquid.Xarrays[i])

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
        append!(X_inner_final,view(X_inner, imin:length(X_inner)))
        append!(X_inner_final,view(X_inner, 1:imin-1))
        append!(θ_inner_final,view(θ_inner, imin:length(θ_inner)))
        append!(θ_inner_final,view(θ_inner, 1:imin-1))
        append!(H_inner_final,view(H_inner, imin:length(H_inner)))
        append!(H_inner_final,view(H_inner, 1:imin-1))
    end

    imin = argmin(X_inner_pres)
    X_inner_pres_final =Array{Float64}(undef, 0)
    P_inner_final =Array{Float64}(undef, 0)
    if imin != 1
        append!(X_inner_pres_final,view(X_inner_pres, imin:length(X_inner_pres)))
        append!(X_inner_pres_final,view(X_inner_pres, 1:imin-1))
        append!(P_inner_final,view(P_inner, imin:length(P_inner)))
        append!(P_inner_final,view(P_inner, 1:imin-1))
    end

# # println(minimum(θ_inner_final))
# # println(θ_inner_final[1])
#
#     # println(X_inner_final)
# # return X_inner_final
#     θ_interp_liquidtowall = LinearInterpolation(X_inner_final, θ_inner_final);
#
#     H_interp_liquidtowall = LinearInterpolation(X_inner_final, H_inner_final);
#
#     θ_interp_walltoliquid = LinearInterpolation(sys.wall.Xarray, sys.wall.θarray);
#
#     P_interp_liquidtowall = LinearInterpolation(X_inner_pres_final, P_inner_final);

    return X_inner_final,θ_inner_final
end
