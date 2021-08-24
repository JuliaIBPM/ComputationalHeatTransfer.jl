export constructmapping,sys_interpolation


# """
#     Create the mapping index of liquid slug corresponding to the wall
#     return Array{Tuple{Float64,Float64},1}
#
#     Xwall   :: Array{Float64,1}
#     Xarrays :: Array{Array{Float64,1},1}
# """
#
# function walltoliquidmapping(Xwall,Xarrays,closedornot,L)
#
#
# if closedornot == false
# for i = 1:length(Xarrays)
#
#     if Xarrays[i][end] < Xwall
#
#     else
#         for j = 1:length(Xarrays[i])
#             if (j == 1 && Xarrays[i][j] >= Xwall)
#                     return (i,-1)
#                     elseif Xarrays[i][j] >= Xwall && Xarrays[i][j-1] <= Xwall
#                    return (i,j)
#             end
#         end
#     end
#
# end
#
#     return (length(Xarrays)+1,-1) # for closed end tube
#
# end
#
# if closedornot == true
#
#     for i = 1:length(Xarrays)
#
# #         println((Xarrays[i][end] < Xarrays[i][1]) && ifamongone(Xwall,Xarrays[i],L))
#
#
#         # firstly deal with the case in a crossing starting point liquid slug
#         if (Xarrays[i][end] < Xarrays[i][1]) && ifamongone(Xwall,Xarrays[i],L)
#
#
#             if (Xwall <= Xarrays[i][end])
#
#                 for index = 1:length(Xarrays[i])-1
#                     j = length(Xarrays[i]) - index + 1
#
#                     if (Xarrays[i][j] <= Xwall) || ((Xarrays[i][j] >= Xwall) && (Xarrays[i][j] <= Xarrays[i][j-1]))
#                         return (i,j)
#                     end
#                 end
#
#             end
#
#             if (Xwall >= Xarrays[i][1])
#
#                 for j = 1:length(Xarrays[i])-1
#
#                     if (Xarrays[i][j] >= Xwall) || ((Xarrays[i][j] <= Xwall) && (Xarrays[i][j+1] <= Xarrays[i][j]))
#                         return (i,j)
#                     end
#                 end
#
#             end
#         end
#
#
#
#         if (Xarrays[i][end] < Xarrays[i][1]) && !ifamongone(Xwall,Xarrays[i],L)
#
#             if ((i > 1) && (Xwall >= Xarrays[i-1][end])) || ((i == 1) && (Xwall >= Xarrays[end][end]))
#                  return (i,-1)
#             end
#
#         end
#
#         # then deal with the normal liquid slug
#         if (Xarrays[i][end] >= Xarrays[i][1]) && ifamongone(Xwall,Xarrays[i])
#
#             for j = 1:length(Xarrays[i])
#
#                 if (Xarrays[i][j] >= Xwall)
#                     return (i,j)
#                 end
#             end
#         end
#
#
#         # then deal with the normal liquid slug
#         if (Xarrays[i][end] >= Xarrays[i][1]) && !ifamongone(Xwall,Xarrays[i])
#             if ((i > 1) && (Xarrays[i][1] <= Xarrays[i-1][end])) || ((i == 1) && (Xarrays[i][1] <= Xarrays[end][end]))
#                 if ((i > 1) && ((Xwall >= Xarrays[i-1][end]) || (Xwall <= Xarrays[i][1])) || ((i == 1) && ((Xwall >= Xarrays[end][end]) || (Xwall <= Xarrays[i][1]))))
#
#                             return (i,-1)
#
#                 end
#             end
#
#
#
#             if ((i > 1) && (Xarrays[i][1] >= Xwall) && (Xarrays[i-1][end] <= Xwall)) || ((i == 1) && (Xarrays[i][1] >= Xwall) && (Xarrays[end][end] <= Xwall))
#                     return (i,-1)
#                 end
#
#         end
#
#
#     end
#
#     return ("error")
#
# end
#
# end
#
# """
#     Create the mapping index of liquid slug corresponding to the wall
#     return Array{Tuple{Float64,Float64},1}
#
#     Xwallarray   :: Array{Float64,1}
#     Xliquidone   :: Array{Float64,1}
# """
#
# function liquidtowallmapping(Xliquidone,Xwallarray)
#
# for i = 2:length(Xwallarray)
#     if Xwallarray[i] >= Xliquidone && Xwallarray[i-1] <= Xliquidone
#         return (i)
#     end
# end
#     return (-1) # for closed end tube
# end
#
# """
#     A lazy way to tranfer Array{Array{Float64,1},1} to Array{Array{Int64,1},1}
#
# """
#
# function truncate(Xarrays::Array{Array{Float64,1},1})
#
#     integerXarrays = Array{Array{Int64,1},1}(undef, length(Xarrays))
#
#     for i =1:length(Xarrays)
#         integerXarrays[i] = trunc.(Int, Xarrays[i])
#     end
#     return integerXarrays
# end
#
# """
#     A function to create mapping indexes for Xarrays and Xwallarray
#
#     walltoliquid : Array{Tuple{Int64,Int64},1}
#     liquidtowall : Array{Array{Int64,1},1}
# """
#
# function constructmapping(Xarrays,Xwallarray,closedornot,L)
#     walltoliquid = Array{Tuple{Int64,Int64},1}(undef, length(Xwallarray))
#
#     for i = 1:length(Xwallarray)
#         walltoliquid[i] = walltoliquidmapping(Xwallarray[i],Xarrays,closedornot,L)
#     end
#
#     liquidtowall = truncate(Xarrays)
#
#     for i = 1:length(Xarrays)
#         for j = 1:length(Xarrays[i])
#             liquidtowall[i][j] = liquidtowallmapping(Xarrays[i][j],Xwallarray)
#         end
#     end
#
#     return walltoliquid,liquidtowall
# end

function sys_interpolation(sys)
    X_inner = Array{Float64}(undef, 0)
    θ_inner = Array{Float64}(undef, 0)
    H_inner = Array{Float64}(undef, 0)
    X_inner_pres  = Array{Float64}(undef, 0)
    P_inner = Array{Float64}(undef, 0)
    mk = 0
    θ = nondi_PtoT.(sys.vapor.P)
    P = sys.vapor.P
    # θ = real.((sys.vapor.P .+ 0im).^((sys.vapor.γ-1)/sys.vapor.γ)) # isentropic
    H_vapor = sys.vapor.Hδ ./ sys.vapor.δ
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

# println(minimum(θ_inner_final))
# println(θ_inner_final[1])

    # println(X_inner_final)
# return X_inner_final
    θ_interp_liquidtowall = LinearInterpolation(X_inner_final, θ_inner_final);

    H_interp_liquidtowall = LinearInterpolation(X_inner_final, H_inner_final);

    θ_interp_walltoliquid = LinearInterpolation(sys.wall.Xarray, sys.wall.θarray);

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
    θ = nondi_PtoT.(sys.vapor.P)
    P = sys.vapor.P
    # θ = real.((sys.vapor.P .+ 0im).^((sys.vapor.γ-1)/sys.vapor.γ)) # isentropic
    H_vapor = sys.vapor.Hδ ./ sys.vapor.δ
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
