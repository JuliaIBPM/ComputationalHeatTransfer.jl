export ODE_innertube

function ODE_innertube(u,p,t)

    # println(dt)

    index_dynamics_end = findfirst(x->abs(x+1e10) <= 10^(-1), u)

    # println(u[1:100])

    newsys = getcurrentsys(u,p)

    dynamicsdu = dynamicsmodel(u[1:index_dynamics_end-1],newsys)

    liquiddu = duliquidθtovec(liquidmodel(newsys))

    du = [dynamicsdu;liquiddu]

    return(du)

end
