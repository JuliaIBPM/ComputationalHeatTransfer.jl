# module Boiling
#
# using ..Systems,..Tools

export boiling_condition,boiling_affect!,nucleateboiling

function boiling_condition(u,t,integrator)
    return (abs(mod(t,0.01)-0.01) < 1e-6) || mod(t,0.01) < 1e-6
end

function boiling_affect!(integrator)
#     println("hahhaha")
    Δθthreshold = 0.01

    p = deepcopy(getcurrentsys(integrator.u,integrator.p))

    for i = 1:length(p.wall.Xstations)

        if ifamong(p.wall.Xstations[i], p.liquid.Xp)
            Δθ = getsuperheat(p.wall.Xstations[i],p)
            if Δθ > Δθthreshold

                # get the insert pressure.
                wallindex = getoneXarrayindex(p.wall.Xstations[i],p.wall.Xarray)
                liquidindex = p.mapping.walltoliquid[wallindex]
                θinsert = p.liquid.θarrays[liquidindex[1]][liquidindex[2]]
                Pinsert = θinsert.^(γ/(γ-1))

                p = nucleateboiling(p,(p.wall.Xstations[i]-p.tube.d/2,p.wall.Xstations[i]+p.tube.d/2),Pinsert) # P need to be given from energy equation
            end
        end




    end

    Lvaporplug = XptoLvaporplug(p.liquid.Xp,p.tube.L,p.tube.closedornot)
    M = p.vapor.P.^(1/p.vapor.γ).* Lvaporplug

    unew=[XMδtovec(p.liquid.Xp,p.liquid.dXdt,M,p.vapor.δ); wallθtovec(p.wall.θarray); liquidθtovec(p.liquid.θarrays)];

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
    γ = sys.vapor.γ
    Xarrays = sys.liquid.Xarrays
    θarrays = sys.liquid.θarrays
    closedornot = sys.tube.closedornot

    Lvaporplug =    XptoLvaporplug(Xp,sys.tube.L,sys.tube.closedornot)
    M = P.^(1/γ).* Lvaporplug

    index = getinsertindex(Xp,Xvapornew,closedornot)





    # Xvapornew = (1.50,1.60)
    # Pinsert = 1.5
    ρinsert = Pinsert.^(1/γ)
    Linsert = Xvapornew[end] - Xvapornew[1]
    Minsert = Pinsert.^(1/γ).* Linsert

    # index = getinsertindex(Xp,Xvapornew)
    δnew = getnewδ(δ,index,Xvapornew,ρ,d,Minsert,closedornot) # mass conservation
    Xpnew = getnewXp(Xp,index,Xvapornew,closedornot)
    Mnew = getnewM(M,index,Minsert,closedornot)

    Lvaporplugnew = XptoLvaporplug(Xpnew,L,closedornot)
    Pnew = (Mnew./Lvaporplugnew).^γ

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

    walltoliquid,liquidtowall = constructmapping(sysnew.liquid.Xarrays ,sysnew.wall.Xarray, sysnew.tube.closedornot, sysnew.tube.L)
    # print(typeof(walltoliquid),"\n",typeof(walltoliquid),"\n")
    sysnew.mapping = Mapping(walltoliquid,liquidtowall)

return sysnew
end


function getnewδ(δ,index,Xvapornew,ρ,d,Minsert,closedornot)
    Linsert = Xvapornew[end] - Xvapornew[1]

    crossAfilms = getcrossAδ.([d],δ)
    insertcrossA = (closedornot && index == length(crossAfilms)) ?  0.5(crossAfilms[index] + crossAfilms[1]) - Minsert/ρ/Linsert : 0.5(crossAfilms[index] + crossAfilms[index+1]) - Minsert/ρ/Linsert
    insert!(crossAfilms,index+1,insertcrossA)

    δnew = crossAtoδ.([d],crossAfilms)
end

function getnewMδ(δ,index,Xvapornew,ρ,d,Minsert,closedornot)

end


function crossAtoδ(d,crossA)
    C = crossA # nondimensional area
    δ = 1 - sqrt(1-C)
end


function getcrossAδ(d,δ)
    crossAouter = 1
    crossAinner = (1-δ)^2
    crossAδ     = crossAouter - crossAinner
end



 function getnewθarrays(index,Xp,Xpnew,Xarrays,θarrays,L,closedornot)
    θarraysnew = deepcopy(θarrays)
    Xarraysnew = deepcopy(Xarrays)
    arrayindex = getarrayindex(Xpnew[index][2],Xarrays[index])

    θarraysnewleft = θarrays[index][1:arrayindex]
    θarraysnewright= θarrays[index][arrayindex+1:end]

    splice!(θarraysnew, index)
    insert!(θarraysnew, index,θarraysnewleft)
    insert!(θarraysnew, index+1,θarraysnewright)
end


function getnewXarrays(index,Xp,Xpnew,Xarrays,L,closedornot)
    Xarraysnew = deepcopy(Xarrays)
    arrayindex = getarrayindex(Xpnew[index][2],Xarrays[index])

    Xarraysnewleft = LinRange(Xpnew[index][1],Xpnew[index][2],arrayindex)
    Xarraysnewright= LinRange(Xpnew[index+1][1],Xpnew[index+1][2],length(Xarrays[index])-arrayindex)

    splice!(Xarraysnew, index)
    insert!(Xarraysnew, index,Xarraysnewleft)
    insert!(Xarraysnew, index+1,Xarraysnewright)
end


function getinsertindex(Xp,Xvapornew,closedornot)

for index = 1:length(Xp)
    if Xp[index][1] <= Xvapornew[1] && Xp[index][2] >= Xvapornew[2]
        return index
end
    end
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

    wallindex = getoneXarrayindex(Xstation,sys.wall.Xarray)
    liquidindex = sys.mapping.walltoliquid[wallindex]

    Δθ = sys.wall.θarray[wallindex] - sys.liquid.θarrays[liquidindex[1]][liquidindex[2]]

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


#
# end
