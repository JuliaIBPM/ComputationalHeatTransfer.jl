# module Tools

export getheight,
XMtovec,XMδtovec,vectoXM,vectoXMδ,
XptoLvaporplug,XptoLliquidslug,getXpvapor,XpvaportoLoverlap,
ifamongone,ifamong,settemperature!,laplacian,constructXarrays,
walltoliquidmapping,liquidtowallmapping,truncate,constructmapping,
duliquidθtovec,duwallθtovec,liquidθtovec,wallθtovec,updateXarrays,
getcurrentsys,getwallWearray,modX!,getoneXarrayindex,getsuperheat,
getmerge_flags,onePtooneT

# using ..Systems
# using LinearAlgebra

"""
    (open loop only)
    This function is a sub-function of getheight. This function is to get the actural physical height for one interface
        X     ::   the location of one interface
        L2D   ::   the length of one bend to another bend (the length in 2D)
        angle ::   the inclination angle
"""

function getoneheight(X::Float64,L2D::Float64,angle::Float64)

    Integer(mod(div(X,L2D),2.0)) == 0 ? L2D - mod(X,L2D) : mod(X,L2D)

end

"""
    (open loop only)
    This function is to get the actural physical heights for all interfaces
        Xp    ::   the locations of all interfaces
        L2D   ::   the length of one bend to another bend (the length in 2D)
        angle ::   the inclination angle
"""

function getheight(Xp::Array{Tuple{Float64,Float64},1},L2D::Float64,angle::Float64)

    height=deepcopy(Xp)

    for i =1:length(Xp)
        height[i]=(getoneheight(Xp[i][1],L2D,angle), getoneheight(Xp[i][end],L2D,angle))
    end

    return height
end

"""
    (open loop only)
    This function is to transform Xp, dXdt of the interface, and M of the vapor to form our state vector u
        Xp    ::   the locations of all interfaces
        dXdt  ::   the 1D velocity of all interfaces
        M     ::   the mass of all vapors
"""

function XMtovec(Xp::Array{Tuple{Float64,Float64},1},dXdt::Array{Tuple{Float64,Float64},1},M::Array{Float64,1})
    if (length(Xp) == length(dXdt)) && (length(Xp) + 1 == length(M))

        u=zeros(5*length(Xp)+1)

        for i = 1:length(Xp)

            # input Xp
            u[2*i-1] = Xp[i][1]
            u[2*i] = Xp[i][end]

            # input dXdt
            u[2*length(Xp) + 2*i-1] = dXdt[i][1]
            u[2*length(Xp) + 2*i] = dXdt[i][end]
        end

        for i = 1:length(M)
            # input M
            u[4*length(Xp) + i] = M[i]
        end

        return u
    end

    if (length(Xp) == length(dXdt)) && (length(Xp) == length(M))

        u=zeros(5*length(Xp))

        for i = 1:length(Xp)

            # input Xp
            u[2*i-1] = Xp[i][1]
            u[2*i] = Xp[i][end]

            # input dXdt
            u[2*length(Xp) + 2*i-1] = dXdt[i][1]
            u[2*length(Xp) + 2*i] = dXdt[i][end]
        end

        for i = 1:length(M)
            # input M
            u[4*length(Xp) + i] = M[i]
        end

        return u
    end
            println("the lengthes of X and dXdt and M do not match!")
            return "error"

end

"""
    (open loop only)
    This function is to transform Xp, dXdt of the interface, and M of the vapor to form our state vector u
        Xp    ::   the locations of all interfaces
        dXdt  ::   the 1D velocity of all interfaces
        M     ::   the mass of all vapors
        δ     ::   the thickness of film in all vapors
"""

function XMδtovec(Xp,dXdt,M,δ)

    return ([XMtovec(Xp,dXdt,M);δ])
end

"""
    (open loop only)
    This function is to transform Xp, dXdt of the interface, and M of the vapor to form our state vector u
        u    ::   the dynamic portion of state vector
"""

function vectoXM(u::Array{Float64,1})

    if mod(length(u),5) == 1
        maxindex = Integer( (length(u) - 1)/5 )

        Xp = map(tuple, zeros(maxindex), zeros(maxindex))
        dXdt = map(tuple, zeros(maxindex), zeros(maxindex))
        M = zeros(maxindex+1)

        for i = 1:maxindex

            # input Xp
            Xp[i] = (u[2*i-1],u[2*i])

            # input dXdt
            dXdt[i] = (u[2*maxindex + 2*i-1],u[2*maxindex + 2*i])
        end

        for i = 1:(maxindex+1)

            # input M
            M[i] = u[4*maxindex + i]

        end

        return Xp,dXdt,M
    end

    if mod(length(u),5) == 0
        maxindex = div(length(u),5)

        Xp = map(tuple, zeros(maxindex), zeros(maxindex))
        dXdt = map(tuple, zeros(maxindex), zeros(maxindex))
        M = zeros(maxindex)

        for i = 1:maxindex

            # input Xp
            Xp[i] = (u[2*i-1],u[2*i])

            # input dXdt
            dXdt[i] = (u[2*maxindex + 2*i-1],u[2*maxindex + 2*i])
        end

        for i = 1:maxindex

            # input M
            M[i] = u[4*maxindex + i]

        end

        return Xp,dXdt,M
    end

return "error"

end

"""
    (open loop only)
    This function is to transform Xp, dXdt of the interface, and M of the vapor to form our state vector u
        u    ::   the dynamic portion of state vector
"""

function vectoXMδ(u::Array{Float64,1})
if mod(length(u),6) == 2
    maxindex = Integer( (length(u) - 2)/6 )

    Xp = map(tuple, zeros(maxindex), zeros(maxindex))
    dXdt = map(tuple, zeros(maxindex), zeros(maxindex))
    M = zeros(maxindex+1)
    δ = zeros(maxindex+1)

    for i = 1:maxindex

        # input Xp
        Xp[i] = (u[2*i-1],u[2*i])

        # input dXdt
        dXdt[i] = (u[2*maxindex + 2*i-1],u[2*maxindex + 2*i])
    end

    for i = 1:(maxindex+1)

        # input M
        M[i] = u[4*maxindex + i]
        δ[i] = u[5*maxindex + 1 + i]
    end

    return Xp,dXdt,M,δ
end

if mod(length(u),6) == 0
    maxindex = div(length(u),6)

    Xp = map(tuple, zeros(maxindex), zeros(maxindex))
    dXdt = map(tuple, zeros(maxindex), zeros(maxindex))
    M = zeros(maxindex)
    δ = zeros(maxindex)

    for i = 1:maxindex

        # input Xp
        Xp[i] = (u[2*i-1],u[2*i])

        # input dXdt
        dXdt[i] = (u[2*maxindex + 2*i-1],u[2*maxindex + 2*i])
    end

    for i = 1:maxindex

        # input M
        M[i] = u[4*maxindex + i]
        δ[i] = u[5*maxindex + i]
    end

    return Xp,dXdt,M,δ
end

return "error"

end

"""
    (open loop only)
    This function is to transform Xp of every interface, and L of the tube to form an array of vapor length
        Xp    ::   the locations of all interfaces
        L     ::   the length of the 1D tube
"""

function XptoLvaporplug(Xp::Array{Tuple{Float64,Float64},1},L::Float64,closedornot::Bool)

if closedornot == false
    maxindex = length(Xp) + 1
    Lvaporplug = zeros(maxindex)

    Lvaporplug[1] = Xp[1][1]-0.0
    Lvaporplug[end] = L-Xp[end][end]

    if maxindex > 2
        for i = 2:maxindex-1

            Lvaporplug[i] = Xp[i][1] - Xp[i-1][end]

        end
    end

    return Lvaporplug
end

if closedornot == true
    maxindex = length(Xp)
    Lvaporplug = zeros(maxindex)

    Lvaporplug[1] = mod((Xp[1][1]-Xp[end][end]),L)
    # Lvaporplug[end] = L-Xp[end][end]

    if maxindex > 1
        for i = 2:maxindex

            Lvaporplug[i] = mod((Xp[i][1] - Xp[i-1][end]),L)

        end
    end

    return Lvaporplug
end

end

"""
    (open loop only)
    This function is to transform Xp of every interface to form an array of liquid length
        Xp    ::   the locations of all interfaces
"""

function XptoLliquidslug(Xp::Array{Tuple{Float64,Float64},1},L::Float64)

    Lliquidslug = zeros(length(Xp))


        for i = 1: length(Xp)

        Lliquidslug[i] = mod((Xp[i][end] - Xp[i][1]),L)

        end

    return Lliquidslug

end

"""
    (open loop only)
    The Xp was coupled by every liquid slug. For instance, if there is one liquid slug. Xp is a one-element tuple (Xp[1][1], Xp[1][2]).
    But sometimes we need Xp to be coupled by every vapor plug. For one liquid slug, we have two vapor plugs.
    So by adding 0 and L at the beginning and the end,
    we construct a two-element tuple ((0.0,Xp[1][1]) and ((Xp[1][2],L). Generally, for every N-element Xp, we construct an N+1 element Xpvapor
        Xp    ::   the locations of all interfaces, each element means a liquid slug.
        L     ::   the length of the 1D tube
"""

function getXpvapor(Xp,L,closedornot)

    Xpvapor=deepcopy(Xp)

    if closedornot == false
        Xpvapor[1]=(0.0,Xp[1][1])

        for i = 2:(length(Xp))
            Xpvapor[i]=(Xp[i-1][end],Xp[i][1])
        end

        push!(Xpvapor,(Xp[end][end],L))
    end

    if closedornot == true
        Xpvapor[1]=(Xp[end][end],Xp[1][1])

        for i = 2:(length(Xp))
            Xpvapor[i]=(Xp[i-1][end],Xp[i][1])
        end
    end

    return Xpvapor
end

"""
    (depreciated)
    This function aims to get the overlapping length between the vapor plug and each evaporator/condenser section and sum them up.
    It can sum up all evaporator sections or condensor sections respectively, but it cannot sum both of them at the same time.
    For example:

    XpvaportoLoverlap(Xpvapor,Xe) sums all the overlapping regions for evaporator
    XpvaportoLoverlap(Xpvapor,Xc) sums all the overlapping regions for condensor

        Xpvapor    ::   the locations of all interfaces, each element means a vapor plug.
        Xce        ::   the locations of all evaporators/condensors, each element is a evaporators/condensors section
"""

function XpvaportoLoverlap(Xpvapor,Xce)
    Loverlap = zeros(length(Xpvapor))
    for i = 1:length(Xpvapor)
        for j = 1:length(Xce)
            if ifoverlap(Xpvapor[i],Xce[j])
            Loverlap[i] += oneoverlap(Xpvapor[i],Xce[j])
            end
        end
    end

    return Loverlap
end

"""
    (depreciated)
    This is a sub-function of XpvaportoLoverlap to solve the overlapping length between a single vapor section and a single evaporator/condensor section

    oneXpvapor ::  the locations of two ends of a vapor plug
    oneXce     ::  the locations of two ends of a evaporator/condensor
"""

function oneoverlap(oneXpvapor,oneXce)
    return min(oneXpvapor[end],oneXce[end]) - max(oneXpvapor[1],oneXce[1])
end

"""
    (depreciated)
    This is a sub-function to determine if there is a overlapping region between a single vapor section and a single evaporator/condensor section

    oneXpvapor ::  the locations of two ends of a vapor plug
    oneXce     ::  the locations of two ends of a evaporator/condensor
"""

function ifoverlap(oneXpvapor,oneXce)
    return ( (oneXpvapor[end] >= oneXce[1]) || (oneXpvapor[1] <= oneXce[end]) ) && (oneXpvapor[end] > oneXce[1]) && (oneXpvapor[1] < oneXce[end])
end

"""
    This is a general sub-function of ifamong to determine if the value is in the range

    value ::  a value
    range ::  a tuple
"""

function ifamongone(value::Float64, range::Tuple{Float64,Float64})
    return (value >= range[1]) && (value <= range[end]) ? true : false
end

"""
    This is a function for a closedloop to determine if the value in in the range that crosses the end point

    value ::  a value
    range ::  an array
"""

function ifamongone(value::Float64, range::Array{Float64,1}, L::Float64)
    return ((value >= range[1]) && (value <= L)) || ((value <= range[end]) && (value >= 0.0)) ? true : false
end

"""
    This is a function to see if the value in in the range

    value ::  a value
    range ::  an array
"""

function ifamongone(value::Float64, range::Array{Float64,1})
    return (value >= range[1]) && (value <= range[end]) ? true : false
end


"""
    This is a general function to determine if the value is in any of an array of range

    value ::  a value
    range ::  an array of tuple
"""

function ifamong(value, X)

    return Bool(sum(ifamongone.(value,X)))
end

"""
    This is a temporary function to initialize wall temperature field

    θᵣ     ::  one temperature value
    xvalue ::  one x value
    sys0   ::  system struct to get Xe and Xc
"""


function settemperature!(θᵣ,xvalue,sys0)

    if ifamong(xvalue, sys0.evaporator.Xe)
        θᵣ = sys0.evaporator.θe

    elseif ifamong(xvalue, sys0.condenser.Xc)
        θᵣ = sys0.condenser.θc
    end

    return θᵣ

end

"""
    This is a function to get the laplacian of a vector field u

    For now zero-gradient boundary condition is used.

    u    ::  an array
"""


function laplacian(u,periodic=false)
    unew = deepcopy(u)

    dl = ones(length(u)-1)
    dr = dl
    d  = -2*ones(length(u))

    A = Tridiagonal(dl, d, dr)

    unew = A*u

    if periodic

        unew[1]   = u[2] - 2*u[1] + u[end]

        unew[end] = u[1] - 2*u[end] + u[end-1]

    else

    # zero gradient B.C.
    unew[1]   = unew[2]
    unew[end] = unew[end-1]
    end

    return (unew)
end

"""
    initialize X and θ field for every liquid slugs. return Array{Array{Float64, 1}, 1} and Array{Array{Float64, 1}, 1}

    X0       :: Array{Tuple{Float64,Float64},1}
    N        :: Int, the number of cells in the wall (ΔX for liquid equals ΔX for the wall)
    θinitial :: value
    L        :: tube length
"""

function constructXarrays(X0::Array{Tuple{Float64,Float64},1},N,θinitial,L)
    Xarrays=Array{Array{Float64, 1}, 1}(undef, length(X0))

    Lliquid = XptoLliquidslug(X0,L)

    Nliquid =  floor.(Int, N.*Lliquid./L)

    for i = 1:length(Xarrays)
        if X0[i][1] < X0[i][2]
            Xarrays[i] = range(X0[i][1], X0[i][2], length=Nliquid[i])
        else
            Xarrays[i] = range(X0[i][1], X0[i][2]+L, length=Nliquid[i]) .- L
            Xarrays[i] = mod.(Xarrays[i], L)
        end
    end

    θarrays = deepcopy(Xarrays)
    for i = 1:length(θarrays)
        θarrays[i][:] .= θinitial
    end

    return(Xarrays,θarrays)
end

"""
    initialize X and θ field for wall, return Array{Float64, 1} and Array{Float64, 1}

    X0       :: Array{Tuple{Float64,Float64},1}
    N        :: Int, the number of cells in the wall (ΔX for liquid equals ΔX for the wall)
    θinitial :: value
    L        :: tube length
"""

function constructXarrays(L::Float64,N,θinitial)
    Xwallarray = Array{Float64, 1}(undef, N)
    Xwallarray = range(0, L, length=N)

    θwallarray = deepcopy(Xwallarray)
    θwallarray = range(θinitial, θinitial, length=N)

    return(Xwallarray,θwallarray)
end

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


"""
    A bunch of functions to transfer θ to state vector rate du
"""

function duliquidθtovec(duθarrays)
    return vcat(map(duwallθtovec, duθarrays)...)
end

function duwallθtovec(duθwall)
    return [0.0; duθwall]
end

"""
    A bunch of functions to transfer θ to state vector u
"""

function liquidθtovec(θarrays)
    return vcat(map(wallθtovec, θarrays)...)
end

function wallθtovec(θwall)
    return [-1e10; θwall]
end

"""
    When having a new Xp because of dynamics, Xarrays need to be updated, too.
"""

function updateXarrays(Xp,θarrays,L)

    Xarrays = deepcopy(θarrays)

    for i = 1:length(Xarrays)
        if Xp[i][1] < Xp[i][2]
            Xarrays[i] = range(Xp[i][1], Xp[i][2], length=length(Xarrays[i]))
        else
            Xarrays[i] = range(Xp[i][1], Xp[i][2]+L, length=length(Xarrays[i])) .- L
            Xarrays[i] = mod.(Xarrays[i], L)
        end
    end

    # for i = 1:length(Xp)
    #     Xarrays[i] = LinRange(Xp[i][1],Xp[i][2],length(θarrays[i]))
    # end

    return Xarrays
end

"""
    give a new u and an old system, return a new system sysnew
"""

function getcurrentsys(u,sys0)

        indexes = Int64[]
        θliquidrec = Array[]

        for i = 1:length(u)
            if abs(u[i]+1e10) <= 10^(-1)
                push!(indexes,i)
            end
        end


    Xp,dXdt,M,δ = vectoXMδ(u[1:indexes[1]-1])
    modX!(Xp,sys0.tube.L)
    θwallrec = u[indexes[1]+1:indexes[2]-1]

    for i = 1:length(indexes)-2
    push!(θliquidrec, u[indexes[i+1]+1:indexes[i+2]-1])
    end
    push!(θliquidrec, u[indexes[end]+1:end])

    sysnew = deepcopy(sys0)

    sysnew.liquid.Xp = Xp
    sysnew.liquid.dXdt = dXdt
    sysnew.liquid.θarrays = θliquidrec
    sysnew.liquid.Xarrays = updateXarrays(Xp,sysnew.liquid.θarrays,sysnew.tube.L)


    Lvaporplug = XptoLvaporplug(Xp,sys0.tube.L,sys0.tube.closedornot)
    γ = sys0.vapor.γ
    P = real.((M./Lvaporplug .+ 0im).^(γ))
    sysnew.vapor.P = P
    sysnew.vapor.δ = δ

    sysnew.wall.θarray = θwallrec

    walltoliquid, liquidtowall = constructmapping(sysnew.liquid.Xarrays ,sysnew.wall.Xarray, sysnew.tube.closedornot, sysnew.tube.L)
    sysnew.mapping = Mapping(walltoliquid,liquidtowall)

    return sysnew
end

function getwallWearray(Xarray,p::PHPSystem)

    Wearray = zero(deepcopy(Xarray))


    for i = 1:length(Xarray)
        for j = 1:length(p.evaporator.Xe)
            if ifamongone(Xarray[i],p.evaporator.Xe[j])
                Wearray[i] = p.evaporator.We[j]
            end
        end
    end

    return Wearray
end

function modX!(Xp,L)
    for i = 1:length(Xp)
         Xp[i] = mod.(Xp[i],L)
    end

    return Xp
end

function getoneXarrayindex(X,Xarray)
    for i = 1:length(Xarray)
        if (X >= Xarray[i] && X <= Xarray[i+1])
            return i
        end
    end

    return length(Xarray)
end

function getsuperheat(Xstation,sys)

    wallindex = getoneXarrayindex(Xstation,sys.wall.Xarray)
    liquidindex = sys.mapping.walltoliquid[wallindex]

    Δθ = sys.wall.θarray[wallindex] - sys.liquid.θarrays[liquidindex[1]][liquidindex[2]]

    return Δθ
end

function getmerge_flags(δv,sys)

    # only for closed loop tube
    numofliquidslug = length(sys.liquid.Xp)
    numofmergingsite = numofliquidslug
    merge_flags = Array{Bool,1}(undef, numofmergingsite)

    Xpvapor = getXpvapor(sys.liquid.Xp,sys.tube.L,sys.tube.closedornot)

    for i in 1:numofmergingsite
     # merging bubble length threshold
        merge_flags[i] = ((Xpvapor[i][2] - Xpvapor[i][1]) < δv && (Xpvapor[i][2] - Xpvapor[i][1]) >= 0) || ((Xpvapor[i][2] - Xpvapor[i][1] + sys.tube.L) < δv && (Xpvapor[i][2] - Xpvapor[i][1]) < 0) ? true : false

    end

    return merge_flags
end

function onePtooneT(P,constant)
    1/(1-(log(P)/constant))
end

# end
