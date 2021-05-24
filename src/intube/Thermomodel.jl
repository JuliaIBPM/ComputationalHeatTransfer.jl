# module Thermomodel

export dMdtdynamicsmodel,wallmodel,liquidmodel,dynamicsmodel
# zhang2002model!,dMdtzhang2002model,dynamicsmodel

# using ..Systems,..Tools

function dynamicsmodel(u::Array{Float64,1},p::PHPSystem)

if p.tube.closedornot == false

    du = zero(deepcopy(u))

    (Xp,dXdt0,M,δ)=vectoXMδ(u)


    numofliquidslug =  Integer( (length(u) - 2)/6 )
    sys = deepcopy(p)

    γ = sys.liquid.γ
    ω0 = sys.liquid.ω0
    ℘ = sys.liquid.℘
    Lvaporplug = XptoLvaporplug(Xp,sys.tube.L,sys.tube.closedornot)
    height = getheight(Xp,sys.tube.L2D,sys.tube.angle)
    Xpvapor = getXpvapor(Xp,sys.tube.L,sys.tube.closedornot)



    P = real.((M./Lvaporplug .+ 0im).^(γ))
    θ = real.((P .+ 0im).^((γ-1)/γ))


    for i = 1:numofliquidslug
        du[2*i-1] = u[2*numofliquidslug+2*i-1]
        du[2*i] = du[2*i-1]
        du[2*numofliquidslug + 2*i-1] = -32*u[2*numofliquidslug + 2*i-1] - (ω0[1]^2)*(0.5*(height[i][end]-height[i][1])) + ℘[1]*(P[i]-P[i+1])
        # use ℘[1] and ω0[1] for now, in the future change them to one value rather than a array.
        du[2*numofliquidslug + 2*i] = du[2*numofliquidslug + 2*i-1]

    end

# not sure which one is better
    du[4*numofliquidslug+1:5*numofliquidslug+1] .= dMdtdynamicsmodel(Xpvapor,θ,sys)
    du[5*numofliquidslug+2:end] .= [0.0]

    return du
end

if p.tube.closedornot == true

        du = zero(deepcopy(u))

        (Xp,dXdt0,M,δ)=vectoXMδ(u)


        numofliquidslug =  Integer( (length(u))/6 )
        sys = deepcopy(p)

        γ = sys.liquid.γ
        ω0 = sys.liquid.ω0
        ℘ = sys.liquid.℘
        Lvaporplug = XptoLvaporplug(Xp,sys.tube.L,sys.tube.closedornot)
        height = getheight(Xp,sys.tube.L2D,sys.tube.angle)
        Xpvapor = getXpvapor(Xp,sys.tube.L,sys.tube.closedornot)



        P = real.((M./Lvaporplug .+ 0im).^(γ))
        θ = real.((P .+ 0im).^((γ-1)/γ))

        for i = 1:numofliquidslug
            du[2*i-1] = u[2*numofliquidslug+2*i-1]
            du[2*i] = du[2*i-1]

            if i != numofliquidslug
                du[2*numofliquidslug + 2*i-1] = -32*u[2*numofliquidslug + 2*i-1] - (ω0[1]^2)*(0.5*(height[i][end]-height[i][1])) + ℘[1]*(P[i]-P[i+1])
            else
                du[2*numofliquidslug + 2*i-1] = -32*u[2*numofliquidslug + 2*i-1] - (ω0[1]^2)*(0.5*(height[i][end]-height[i][1])) + ℘[1]*(P[i]-P[1])
            end
            # use ℘[1] and ω0[1] for now, in the future change them to one value rather than a array.
            du[2*numofliquidslug + 2*i] = du[2*numofliquidslug + 2*i-1]

        end

    # not sure which one is better
        du[4*numofliquidslug+1:5*numofliquidslug] .= dMdtdynamicsmodel(Xpvapor,θ,sys)
        du[5*numofliquidslug+1:end] .= [0.0]

        return du
end


end

function dMdtdynamicsmodel(Xpvapor::Array{Tuple{Float64,Float64},1},θ::Array{Float64,1},sys::PHPSystem)

    dMdt=zeros(length(Xpvapor))

    He = sys.evaporator.He

    Hc = sys.condenser.Hc

        dx = sys.wall.Xarray[2]-sys.wall.Xarray[1]


    for i = 1:length(Xpvapor)
        indexes = findall( x -> (mod(x[1],length(Xpvapor)) == mod(i,length(Xpvapor)) && (x[end] == -1)), sys.mapping.walltoliquid)

            for j in length(indexes)

                # print("index",indexes,"\n")
                # print("θarray",sys.wall.θarray,"\n")

                dMdt[i] += He*dx*(sys.wall.θarray[indexes[j]] - θ[i])

            end


    end



    return dMdt

end


function wallmodel(θarray::Array{Float64,1},p::PHPSystem)
    sys = deepcopy(p)

    du = zero(deepcopy(θarray))

    γ = sys.vapor.γ
    Hₗ = sys.liquid.Hₗ
    He = sys.evaporator.He
    dx = sys.wall.Xarray[2]-sys.wall.Xarray[1]

    Xarray = sys.wall.Xarray
    Wearray = getwallWearray(Xarray,sys)

    Hwc = sys.condenser.Hwc
    θc  = sys.condenser.θc
    Xc  = sys.condenser.Xc
    hevisidec=ifamong.(Xarray,[Xc])

    H = zero(deepcopy(θarray))
    θarray_temp_flow = zero(deepcopy(θarray))

        for i = 1:length(θarray)

            index = sys.mapping.walltoliquid[i]

            if index[2] == -1

                if sys.tube.closedornot == false
                    P = sys.vapor.P[index[1]]
                end


                if sys.tube.closedornot == true
                    if index[1] > length(sys.vapor.P)
                        P = sys.vapor.P[index[1]-length(sys.vapor.P)]
                    else
                        P = sys.vapor.P[index[1]]
                    end
                end

                θarray_temp_flow[i] = real.((P .+ 0im).^((γ-1)/γ))

                H[i] = He
            else
                θliquidarrays = sys.liquid.θarrays
                θarray_temp_flow[i] = θliquidarrays[index[1]][index[2]]

                H[i] = Hₗ
            end
        end


        du = sys.wall.α .* laplacian(θarray,sys.tube.closedornot) ./ dx ./ dx + H .* (θarray_temp_flow - θarray) .* dx + Wearray .* dx + hevisidec .* Hwc .* (θc .- θarray) .* dx
        return du
end

function liquidmodel(θarrays,p::PHPSystem)
    sys = deepcopy(p)

    du = zero.(deepcopy(θarrays))

    γ = sys.vapor.γ
    Hₗ = sys.liquid.Hₗ



    θarray_temp_wall = zero.(deepcopy(θarrays))
    for i = 1:length(θarrays)

        dx = sys.wall.Xarray[2]-sys.wall.Xarray[1]

        indexes = sys.mapping.liquidtowall[i]

        for j = 1:length(indexes)
            θarray_temp_wall[i][j] = sys.wall.θarray[indexes[j]]
        end

        du[i] = sys.wall.α .* laplacian(θarrays[i]) ./ dx ./ dx + Hₗ .* (θarray_temp_wall[i] - θarrays[i]) .* dx
    end


    return du
end


"""
    Depreciated
    get the array of evaporator's heat flux along the wall if the 1D tube wall model is used.
"""

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


# """
#     This function is required by "DifferentialEquation.jl" Package.
#         du :: an empty state vector to be derived at the end of this function.
#         u  :: the state vector input from the last time step.
#         p  :: some important parameters which do not belong to u
#         t  :: the current time
#
#     This approach solves the same set of dimensionless governing equations as
#     "https://doi.org/10.1016/S0017-9310(01)00348-9"(Zhang et al. 2002)
#     Instead of solving them piece by piece. Our current approach solve them as a whole system.
# """
#
# function zhang2002model!(du::Array{Float64,1},u::Array{Float64,1},p::PHPSystem,t::Float64)
#
#
#     (Xp,dXdt0,M)=vectoXM(u)
#
#
#     numofliquidslug =  Integer( (length(u) - 1)/5 )
#     sys0 = p
#
#     γ = sys0.liquid.γ
#     ω0 = sys0.liquid.ω0
#     ℘ = sys0.liquid.℘
#     Lvaporplug = XptoLvaporplug(Xp,sys0.tube.L,sys0.tube.closedornot)
#
#     height = getheight(Xp,sys0.tube.L2D,sys0.tube.angle)
#     Xpvapor = getXpvapor(Xp,sys0.tube.L,sys0.tube.closedornot)
#
#
#
#     # get P from M and γ
#     # P = (M./Lvaporplug).^(γ)
#     # P = zeros(size(M))
#     # for i in length(M)
#     #        P[i] = M[i] > 0 ? (M[i]./Lvaporplug[i]).^(γ) : -(-M[i]./Lvaporplug[i]).^(γ)
#     #    end
#     P = real.((M./Lvaporplug .+ 0im).^(γ))
#
#     # get θ from P and γ
#     # θ = zeros(size(P))
#     # for i in length(P)
#     #        θ[i] = P[i] > 0 ? P[i].^((γ-1)/γ) : -(-P[i]).^((γ-1)/γ)
#     #    end
#     θ = real.((P .+ 0im).^((γ-1)/γ))
#     # θ = P.^((γ-1)/γ)
#
#
#     for i = 1:numofliquidslug
#         du[2*i-1] = u[2*numofliquidslug+2*i-1]
#         du[2*i] = du[2*i-1]
#
#         du[2*numofliquidslug + 2*i-1] = -32*u[2*numofliquidslug + 2*i-1] - (ω0[i]^2)*(0.5*(height[i][end]-height[i][1])) + ℘[i]*(P[i]-P[i+1])
#         du[2*numofliquidslug + 2*i] = du[2*numofliquidslug + 2*i-1]
#
#     end
#
#
# # not sure which one is better
#         du[4*numofliquidslug+1:5*numofliquidslug+1] .= dMdtzhang2002model(Xpvapor,θ,sys0)
#         # du[4*numofliquidslug+1:5*numofliquidslug+1] .= dMdtconstantH(Xpvapor,θ,sys0)
#
#     return du
#
# end
#
# """
#     This function is a sub-function of zhang2002model. This function is the phase change part,
#         Xpvapor ::  storing the location of two ends of each vapor
#         θ       ::  the temperature in each vapor.
#         sys0    ::  the struct that stores every needed initial conditions and boundary conditions
#
#     In this case. He and Hc are preset constants.
# """
#
# function dMdtzhang2002model(Xpvapor::Array{Tuple{Float64,Float64},1},θ::Array{Float64,1},sys0::PHPSystem)
#
#     dMdt=zeros(length(Xpvapor))
#
#
#     Xe = sys0.evaporator.Xe
#     He = sys0.evaporator.He
#     θe = sys0.evaporator.θe
#
#     Xc = sys0.condenser.Xc
#     Hc = sys0.condenser.Hc
#     θc = sys0.condenser.θc
#
#     Levapoverlap=XpvaportoLoverlap(Xpvapor,Xe)
#     Lcondoverlap=XpvaportoLoverlap(Xpvapor,Xc)
#
#
#     # May not be right for multi liquid flow
#     for i = 1:length(Xpvapor)
#         if Lcondoverlap[i] < 1e-8
#             dMdt[i] = He * Levapoverlap[i] * (θe-θ[i])
#         else
#             dMdt[i] = -Hc * Lcondoverlap[i] * (θ[i]-θc)
#         end
#     end
#
#     return dMdt
#
# end
#
# function dMdtconstantH(Xpvapor::Array{Tuple{Float64,Float64},1},θ::Array{Float64,1},sys0::PHPSystem)
#
#     dMdt=zeros(length(Xpvapor))
#
#
#     Xe = sys0.evaporator.Xe
#     He = sys0.evaporator.He
#     θe = sys0.evaporator.θe
#
#     Xc = sys0.condenser.Xc
#     Hc = sys0.condenser.Hc
#     θc = sys0.condenser.θc
#
#     Levapoverlap=XpvaportoLoverlap(Xpvapor,Xe)
#     Lcondoverlap=XpvaportoLoverlap(Xpvapor,Xc)
#
#     for i = 1:length(Xpvapor)
#         dMdt[i] = He * Levapoverlap[i] * (θe-θ[i]) -Hc * Lcondoverlap[i] * (θ[i]-θc)
#     end
#     return dMdt
#
# end


#
#
# end
