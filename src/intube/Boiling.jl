# module Boiling
#
# using ..Systems,..Tools

export boiling_affect!,nucleateboiling
# boiling_condition,
function boiling_condition(u,t,integrator)
    # t_to_nondi_t = 0.2831486159429433
    # t_interval = 0.1 * t_to_nondi_t

    t_interval = 0.1

    ϵ = 1e-5


    # println(t)
    # println((abs(mod(t,t_interval)-t_interval) < ϵ))
    # println(mod(t,t_interval) < ϵ)
    return (abs(mod(t,t_interval)-t_interval) < ϵ) || mod(t,t_interval) < ϵ
    # return mod(t,0.01*t_to_nondi_t)
    # mod(t,t_interval) - ϵ
end

function boiling_affect!(integrator)

    Δθthreshold = integrator.p.wall.ΔTthres
    # t_to_nondi_t = 0.2831486159429433
    # """
    # modified here!!!!!!!!!!!!!!!!!!!
    # """

# println("Boiled at " ,integrator.t/t_to_nondi_t)
    # *0.1 """modified here!!!!!!!!!!!!!!!!!!!""" modified here!!!!!!!!!!!!!!!!!!!```


    p = deepcopy(getcurrentsys(integrator.u,integrator.p))

    for i = 1:length(p.wall.Xstations)

        if ifamong(p.wall.Xstations[i], p.liquid.Xp, p.tube.L) && suitable_for_boiling(p,i)

                    # println(ifamong(p.wall.Xstations[i], p.liquid.Xp) && suitable_for_boiling(p,i))

            # println(p.wall.Xstations[i])
            # println(length(p.liquid.Xp))
            Δθ = getsuperheat(p.wall.Xstations[i],p)
            # println(Δθ)
            if Δθ > Δθthreshold
                # println("Boiled! at ",p.wall.Xstations[i], " on ", integrator.t/t_to_nondi_t)
                println("Boiled! on ",p.wall.Xstations[i], " at ", integrator.t)#
# ```modified here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!```
#                 θinsert = p.mapping.θ_interp_walltoliquid.(Xstations[i])
#                 Pinsert = nondi_TtoP.(θinsert)

                Pinsert = p.mapping.P_interp_liquidtowall(Xstations[i])
                # θinsert = nondi_PtoT.(Pinsert)
                θinsert = PtoT.(Pinsert)


                p = nucleateboiling(p,(p.wall.Xstations[i]-p.tube.d,p.wall.Xstations[i]+p.tube.d),Pinsert) # P need to be given from energy equation
            end
        end


        # println(p.vapor.P)

    end

    Lvaporplug = XptoLvaporplug(p.liquid.Xp,p.tube.L,p.tube.closedornot)
    # M = p.vapor.P.^(1/p.vapor.γ).* Lvaporplug
    # M = nondi_PtoD.(p.vapor.P) .* Lvaporplug
    Ac = p.tube.Ac
    M = PtoD.(p.vapor.P) .* Lvaporplug .* Ac


    unew=[XMδtovec(p.liquid.Xp,p.liquid.dXdt,M,p.vapor.δ);liquidθtovec(p.liquid.θarrays)];

#     set_u!(integrator,  unew)
    resize!(integrator.u,length(unew))
    integrator.u = deepcopy(unew)


end

function nucleateboiling(sys,Xvapornew,Pinsert)
    ρ = deepcopy(sys.liquid.ρ)
    d = deepcopy(sys.tube.d)
    Xp = deepcopy(sys.liquid.Xp)
    dXdt = deepcopy(sys.liquid.dXdt)
    δ = deepcopy(sys.vapor.δ)
    P = deepcopy(sys.vapor.P)
    L = sys.tube.L
    Ac = sys.tube.Ac
    # γ = sys.vapor.γ
    Xarrays = sys.liquid.Xarrays
    θarrays = sys.liquid.θarrays
    closedornot = sys.tube.closedornot


    Lvaporplug =    XptoLvaporplug(Xp,sys.tube.L,sys.tube.closedornot)

    M = PtoD.(P) .* Lvaporplug .* Ac
    # M = nondi_PtoD.(P) .* Lvaporplug
    # M = P.^(1/γ).* Lvaporplug

    index = getinsertindex(Xp,(Xvapornew[2]+Xvapornew[1])/2,sys.tube.L,sys.tube.closedornot)




    # Xvapornew = (1.50,1.60)
    # Pinsert = 1.5
    # ρinsert = nondi_PtoD(Pinsert)
    ρinsert = PtoD(Pinsert)

    Linsert = Xvapornew[end] - Xvapornew[1]
    Minsert = ρinsert .* Linsert



    """let's do constant film thickness for now!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"""
    δnew = insert!(δ,index+1,δ[index])
    # δnew = getnewδ(δ,index,Xvapornew,ρ,d,Minsert,closedornot) # mass conservation



    Xpnew = getnewXp(Xp,index,Xvapornew,closedornot)
    # Mnew = getnewM(M,index,Minsert,closedornot)
    #
    # Lvaporplugnew = XptoLvaporplug(Xpnew,L,closedornot)
    # Pnew = nondi_DtoP.(Mnew./Lvaporplugnew)

    Pnew = insert!(P,index+1,Pinsert)

    Xarraysnew = getnewXarrays(index,Xp,Xpnew,Xarrays,L,closedornot)
    θarraysnew = getnewθarrays(index,Xp,Xpnew,Xarrays,θarrays,L,closedornot)

# only for open loop!
    dXdtnew = deepcopy(dXdt) # momentum conservation
    insert!(dXdtnew,index+1,dXdtnew[index])






    sysnew = deepcopy(sys)

    sysnew.liquid.Xp = Xpnew
    sysnew.liquid.dXdt = dXdtnew
    sysnew.liquid.Xarrays = Xarraysnew
    sysnew.liquid.θarrays = θarraysnew
    sysnew.vapor.P = Pnew
    sysnew.vapor.δ = δnew

    # walltoliquid,liquidtowall = constructmapping(sysnew.liquid.Xarrays ,sysnew.wall.Xarray, sysnew.tube.closedornot, sysnew.tube.L)
    # print(Xarraysnew[5],"\n",θarraysnew[5],"\n")
    # sysnew.mapping = Mapping(walltoliquid,liquidtowall)


    θ_interp_walltoliquid, θ_interp_liquidtowall, H_interp_liquidtowall, P_interp_liquidtowall = sys_interpolation(sysnew)
    sysnew.mapping = Mapping(θ_interp_walltoliquid, θ_interp_liquidtowall, H_interp_liquidtowall, P_interp_liquidtowall)

return sysnew
end

#
# function getnewδ(δ,index,Xvapornew,ρ,d,Minsert,closedornot)
#     Linsert = Xvapornew[end] - Xvapornew[1]
#
#     crossAfilms = getcrossAδ.([d],δ)
#     insertcrossA = (closedornot && index == length(crossAfilms)) ?  0.5(crossAfilms[index] + crossAfilms[1]) - Minsert/ρ/Linsert : 0.5(crossAfilms[index] + crossAfilms[index+1]) - Minsert/ρ/Linsert
#     insert!(crossAfilms,index+1,insertcrossA)
#
#     δnew = crossAtoδ.([d],crossAfilms)
# end
#
# function getnewMδ(δ,index,Xvapornew,ρ,d,Minsert,closedornot)
#
# end
#
#
# function crossAtoδ(d,crossA)
#     C = crossA # nondimensional area
#     δ = 1 - sqrt(1-C)
# end
#
#
# function getcrossAδ(d,δ)
#     crossAouter = 1
#     crossAinner = (1-δ)^2
#     crossAδ     = crossAouter - crossAinner
# end



 function getnewθarrays(index,Xp,Xpnew,Xarrays,θarrays,L,closedornot)
    θarraysnew = deepcopy(θarrays)
    Xarraysnew = deepcopy(Xarrays)
    arrayindex = getarrayindex(Xpnew[index][2],Xarrays[index])

    θarraysnewleft = θarrays[index][1:arrayindex]
    # append!(θarraysnewleft, θarrays[index][arrayindex])

    θarraysnewright= θarrays[index][arrayindex+1:end]
    insert!(θarraysnewright, 1, θarrays[index][arrayindex])
    # insert!(θarraysnewright, 1, θarrays[index][arrayindex])


    splice!(θarraysnew, index)
    insert!(θarraysnew, index,θarraysnewleft)
    insert!(θarraysnew, index+1,θarraysnewright)
end


function getnewXarrays(index,Xp,Xpnew,Xarrays,L,closedornot)
    Xarraysnew = deepcopy(Xarrays)
    arrayindex = getarrayindex(Xpnew[index][2],Xarrays[index])
    # println(arrayindex)
    # println(Xpnew[index][2])
    # println(Xarrays[index])

    # Xarraysnewleft = LinRange(Xpnew[index][1],Xpnew[index][2],arrayindex+1)
    # Xarraysnewright= LinRange(Xpnew[index+1][1],Xpnew[index+1][2],length(Xarrays[index])-arrayindex+2)

    Xarraysnewleft = constructoneXarray(Xpnew[index],arrayindex,L)
    Xarraysnewright = constructoneXarray(Xpnew[index+1],length(Xarrays[index])-arrayindex+1,L)

    splice!(Xarraysnew, index)
    insert!(Xarraysnew, index,Xarraysnewleft)
    insert!(Xarraysnew, index+1,Xarraysnewright)
end


function getinsertindex(Xp,Xvapornew_center,L,closedornot)

if !closedornot
    for index = 1:length(Xp)
        if Xp[index][1] <= Xvapornew_center && Xp[index][2] >= Xvapornew_center
            return index
        end
    end
end

if closedornot
    for index = 1:length(Xp)
        if Xp[index][1] <= Xvapornew_center && Xp[index][2] >= Xvapornew_center && Xp[index][2] >= Xp[index][1]
            return index
        end

        if (mod(Xvapornew_center - Xp[index][1],L) < mod(Xp[index][2] - Xp[index][1],L)) && Xp[index][2] <= Xp[index][1]
            return index
        end
    end
end
# println((Xvapornew[2]+Xvapornew[1])/2)
# println(Xp)
        return NaN
end

function getarrayindex(X,Xarray)

for arrayindex = 1:length(Xarray)
    if X >= Xarray[arrayindex] && X <= Xarray[arrayindex+1]
        return arrayindex
end
    end
        return NaN
end

function getnewXp(Xp,index,Xvapornew,closedornot)

    Xpnew = deepcopy(Xp)

    Linsert = Xvapornew[end] - Xvapornew[1]

    insertXp1=(Xp[index][1]-Linsert/2,Xvapornew[1])
    insertXp2=(Xvapornew[2],Xp[index][2]+Linsert/2)

    splice!(Xpnew, index)
    insert!(Xpnew, index,insertXp1)
    insert!(Xpnew, index+1,insertXp2)

    return Xpnew
end

# function getnewδ(δ,Pinsert,index,Xvapornew,closedornot)
#
#     Linsert = Xvapornew[end] - Xvapornew[1]
#     δnew = deepcopy(δ)
#
#     δavg = (δ[index]+δ[index+1])/2
#     insert!(δnew, index+1,δavg)
#
#     return δnew
# end


function getnewM(M,index,Minsert,closedornot)

    Mnew = deepcopy(M)

    insert!(Mnew,index+1,Minsert)

    return Mnew
end


function getsuperheat(Xstation,sys)


    # P = sys.mapping.P_interp_liquidtowall(2.0)
    # # println(Xstation)

    # Δθ = sys.mapping.θ_interp_walltoliquid(Xstation) - nondi_PtoT(sys.mapping.P_interp_liquidtowall(Xstation))
    Δθ = sys.mapping.θ_interp_walltoliquid(Xstation) - PtoT(sys.mapping.P_interp_liquidtowall(Xstation))
    return Δθ
end


```
    get the index of element of Xarray closest to X.
    closed loop considered
```
function getoneXarrayindex(X,Xarray)
    for i = 1:length(Xarray)
        # considering the closed loop B.C.
        if (Xarray[i] >= Xarray[i+1]) && ((Xarray[i] <= X) || (Xarray[i+1] >= X))
            return i
        end
        if (X >= Xarray[i] && X <= Xarray[i+1])
            return i
        end
    end

    return length(Xarray)
end

function suitable_for_boiling(p,i)
    suitable_flag =  false
    index_max = length(p.liquid.Xp)


        index = getinsertindex(p.liquid.Xp,p.wall.Xstations[i],p.tube.L,p.tube.closedornot)
        Lvaporplug = XptoLvaporplug(p.liquid.Xp,p.tube.L,p.tube.closedornot)

        L_vapor_left =  Lvaporplug[index]
        L_vapor_right = (index == index_max) ? Lvaporplug[1] : Lvaporplug[index + 1]

        suitable_flag = (2*p.tube.d < L_vapor_left) && (2*p.tube.d < L_vapor_right) ? true : false



        L_liquid_left =  mod(p.wall.Xstations[i] - p.liquid.Xp[index][1],p.tube.L)
        L_liquid_right = mod(p.liquid.Xp[index][2] - p.wall.Xstations[i],p.tube.L)

        if (5*p.tube.d > L_liquid_left) || (5*p.tube.d > L_liquid_right)
            suitable_flag = false
        end

    return suitable_flag
end


#
# end
