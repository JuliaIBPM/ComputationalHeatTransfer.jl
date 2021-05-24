export ODE_innertube

function ODE_innertube(u,p,t)

        indexes = Int64[]
        θliquidtemp = Array[]

        for i = 1:length(u)
            if abs(u[i]+1e10) <= 10^(-1)
                push!(indexes,i)
            end
        end


    newsys = getcurrentsys(u,p)

    dynamicsdu = dynamicsmodel(u[1:indexes[1]-1],newsys)



    duwalltemp = wallmodel(u[indexes[1]+1:indexes[2]-1],newsys)
    walldu = duwallθtovec(duwalltemp)



    for i = 1:length(indexes)-2
    push!(θliquidtemp, u[indexes[i+1]+1:indexes[i+2]-1])
    end
    push!(θliquidtemp, u[indexes[end]+1:end])
    duliquidtemp = zero.(deepcopy(θliquidtemp))
    duliquidtemp = liquidmodel(θliquidtemp,newsys)
    liquiddu = duliquidθtovec(duliquidtemp)

    du = [dynamicsdu;walldu;liquiddu]

    return(du)

end
