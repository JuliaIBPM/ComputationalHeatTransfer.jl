export ODE_innertube

function ODE_innertube(u,p,t)

        # indexes = Int64[]
        # θarraystemp = Array[]
        #
        # for i = 1:length(u)
        #     if abs(u[i]+1e10) <= 10^(-1)
        #         push!(indexes,i)
        #     end
        # end

    index_dynamics_end = findfirst(x->abs(x+1e10) <= 10^(-1), u)

    newsys = getcurrentsys(u,p)

    dynamicsdu = dynamicsmodel(u[1:index_dynamics_end-1],newsys)



    # walldu = duwallθtovec(wallmodel(newsys))
    # walldu = duwallθtovec(duwalltemp)

    # for i = 1:length(indexes)-2
    # push!(θarraystemp, u[indexes[i+1]+1:indexes[i+2]-1])
    # end
    # push!(θarraystemp, u[indexes[end]+1:end])
    # duliquidtemp = zero.(deepcopy(θarraystemp))
    # duliquidtemp = liquidmodel(newsys)
    liquiddu = duliquidθtovec(liquidmodel(newsys))

    du = [dynamicsdu;liquiddu]

    # println(t)

    return(du)

end
