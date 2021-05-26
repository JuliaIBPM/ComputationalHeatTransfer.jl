export constructmapping


"""
    Create the mapping index of liquid slug corresponding to the wall
    return Array{Tuple{Float64,Float64},1}

    Xwall   :: Array{Float64,1}
    Xarrays :: Array{Array{Float64,1},1}
"""

function walltoliquidmapping(Xwall,Xarrays,closedornot,L)


if closedornot == false
for i = 1:length(Xarrays)

    if Xarrays[i][end] < Xwall

    else
        for j = 1:length(Xarrays[i])
            if (j == 1 && Xarrays[i][j] >= Xwall)
                    return (i,-1)
                    elseif Xarrays[i][j] >= Xwall && Xarrays[i][j-1] <= Xwall
                   return (i,j)
            end
        end
    end

end

    return (length(Xarrays)+1,-1) # for closed end tube

end

if closedornot == true

    for i = 1:length(Xarrays)

#         println((Xarrays[i][end] < Xarrays[i][1]) && ifamongone(Xwall,Xarrays[i],L))


        # firstly deal with the case in a crossing starting point liquid slug
        if (Xarrays[i][end] < Xarrays[i][1]) && ifamongone(Xwall,Xarrays[i],L)


            if (Xwall <= Xarrays[i][end])

                for index = 1:length(Xarrays[i])-1
                    j = length(Xarrays[i]) - index + 1

                    if (Xarrays[i][j] <= Xwall) || ((Xarrays[i][j] >= Xwall) && (Xarrays[i][j] <= Xarrays[i][j-1]))
                        return (i,j)
                    end
                end

            end

            if (Xwall >= Xarrays[i][1])

                for j = 1:length(Xarrays[i])-1

                    if (Xarrays[i][j] >= Xwall) || ((Xarrays[i][j] <= Xwall) && (Xarrays[i][j+1] <= Xarrays[i][j]))
                        return (i,j)
                    end
                end

            end
        end



        if (Xarrays[i][end] < Xarrays[i][1]) && !ifamongone(Xwall,Xarrays[i],L)

            if ((i > 1) && (Xwall >= Xarrays[i-1][end])) || ((i == 1) && (Xwall >= Xarrays[end][end]))
                 return (i,-1)
            end

        end

        # then deal with the normal liquid slug
        if (Xarrays[i][end] >= Xarrays[i][1]) && ifamongone(Xwall,Xarrays[i])

            for j = 1:length(Xarrays[i])

                if (Xarrays[i][j] >= Xwall)
                    return (i,j)
                end
            end
        end


        # then deal with the normal liquid slug
        if (Xarrays[i][end] >= Xarrays[i][1]) && !ifamongone(Xwall,Xarrays[i])
            if ((i > 1) && (Xarrays[i][1] <= Xarrays[i-1][end])) || ((i == 1) && (Xarrays[i][1] <= Xarrays[end][end]))
                if ((i > 1) && ((Xwall >= Xarrays[i-1][end]) || (Xwall <= Xarrays[i][1])) || ((i == 1) && ((Xwall >= Xarrays[end][end]) || (Xwall <= Xarrays[i][1]))))

                            return (i,-1)

                end
            end



            if ((i > 1) && (Xarrays[i][1] >= Xwall) && (Xarrays[i-1][end] <= Xwall)) || ((i == 1) && (Xarrays[i][1] >= Xwall) && (Xarrays[end][end] <= Xwall))
                    return (i,-1)
                end

        end


    end

    return ("error")

end

end

"""
    Create the mapping index of liquid slug corresponding to the wall
    return Array{Tuple{Float64,Float64},1}

    Xwallarray   :: Array{Float64,1}
    Xliquidone   :: Array{Float64,1}
"""

function liquidtowallmapping(Xliquidone,Xwallarray)

for i = 2:length(Xwallarray)
    if Xwallarray[i] >= Xliquidone && Xwallarray[i-1] <= Xliquidone
        return (i)
    end
end
    return (-1) # for closed end tube
end

"""
    A lazy way to tranfer Array{Array{Float64,1},1} to Array{Array{Int64,1},1}

"""

function truncate(Xarrays::Array{Array{Float64,1},1})

    integerXarrays = Array{Array{Int64,1},1}(undef, length(Xarrays))

    for i =1:length(Xarrays)
        integerXarrays[i] = trunc.(Int, Xarrays[i])
    end
    return integerXarrays
end

"""
    A function to create mapping indexes for Xarrays and Xwallarray

    walltoliquid : Array{Tuple{Int64,Int64},1}
    liquidtowall : Array{Array{Int64,1},1}
"""

function constructmapping(Xarrays,Xwallarray,closedornot,L)
    walltoliquid = Array{Tuple{Int64,Int64},1}(undef, length(Xwallarray))

    for i = 1:length(Xwallarray)
        walltoliquid[i] = walltoliquidmapping(Xwallarray[i],Xarrays,closedornot,L)
    end

    liquidtowall = truncate(Xarrays)

    for i = 1:length(Xarrays)
        for j = 1:length(Xarrays[i])
            liquidtowall[i][j] = liquidtowallmapping(Xarrays[i][j],Xwallarray)
        end
    end

    return walltoliquid,liquidtowall
end
