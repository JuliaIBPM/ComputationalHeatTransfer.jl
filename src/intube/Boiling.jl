export boiling_affect!,nucleateboiling,boiling_condition
# boiling_condition,
function boiling_condition(u,t,integrator)
    t_interval = 0.1

    ϵ = 1e-5

    return (abs(mod(t,t_interval)-t_interval) < ϵ) || mod(t,t_interval) < ϵ
end

function boiling_affect!(integrator)

    Δθthreshold = integrator.p.wall.ΔTthres

    p = deepcopy(getcurrentsys(integrator.u,integrator.p))
  
    Δθ_array = getsuperheat.(p.wall.Xstations,[p])
    superheat_flag = Δθ_array .> Δθthreshold

    b_count = 0;
    for i = 1:length(p.wall.Xstations)

        if ifamong(p.wall.Xstations[i], p.liquid.Xp, p.tube.L) && suitable_for_boiling(p,i) && superheat_flag[i]

            Δθ = getsuperheat(p.wall.Xstations[i],p)
          
                push!(Main.boil_hist,[i,integrator.t]);
                b_count += 1;

                Pinsert = p.mapping.P_interp_liquidtowall(p.wall.Xstations[i])

                p = nucleateboiling(p,(p.wall.Xstations[i]-2p.tube.d,p.wall.Xstations[i]+2p.tube.d),Pinsert) # P need to be given from energy equation
        end
    end

    Lvaporplug = XptoLvaporplug(p.liquid.Xp,p.tube.L,p.tube.closedornot)
    Ac = p.tube.Ac
    d = p.tube.d

    δstart = p.vapor.δstart
    δend = p.vapor.δend

    Lfilm_start = p.vapor.Lfilm_start
    Lfilm_end = p.vapor.Lfilm_end

    δarea_start = Ac .* (1 .- ((d .- 2*δstart) ./ d) .^ 2);
    δarea_end = Ac .* (1 .- ((d .- 2*δend) ./ d) .^ 2);

    volume_vapor = Lvaporplug .* Ac - Lfilm_start .* δarea_start - Lfilm_end .* δarea_end
    M = PtoD.(p.vapor.P) .* volume_vapor


    unew=[XMδLtovec(p.liquid.Xp,p.liquid.dXdt,M,p.vapor.δstart,δend,Lfilm_start,Lfilm_end);liquidθtovec(p.liquid.θarrays)];

    resize!(integrator.u,length(unew))
    integrator.u = deepcopy(unew)
end

function nucleateboiling(sys,Xvapornew,Pinsert)
    ρ = deepcopy(sys.liquid.ρ)
    Ac = sys.tube.Ac
    d = deepcopy(sys.tube.d)
    Xp = deepcopy(sys.liquid.Xp)
    dXdt = deepcopy(sys.liquid.dXdt)
    δstart = deepcopy(sys.vapor.δstart)
    δend = deepcopy(sys.vapor.δend)
    # Astart = getδarea(Ac,d,δstart)
    # Aend = getδarea(Ac,d,δend)
    δfilm_deposit = deepcopy(sys.vapor.δfilm_deposit)
    P = deepcopy(sys.vapor.P)
    L = sys.tube.L
    
    Xarrays = sys.liquid.Xarrays
    θarrays = sys.liquid.θarrays
    closedornot = sys.tube.closedornot
    Lfilm_start = deepcopy(sys.vapor.Lfilm_start)
    Lfilm_end = deepcopy(sys.vapor.Lfilm_end)


    Lvaporplug = XptoLvaporplug(Xp,sys.tube.L,sys.tube.closedornot)

    index = getinsertindex(Xp,(Xvapornew[2]+Xvapornew[1])/2,sys.tube.L,sys.tube.closedornot)

    # ρinsert = PtoD(Pinsert)
    Mvapor = getMvapor(sys)

    Linsert = mod(Xvapornew[end] - Xvapornew[1],L)

    """let's do constant film thickness for now!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"""
    δstart_new = insert!(δstart,index+1,δfilm_deposit)
    δend_new = insert!(δend,index+1,δfilm_deposit)

    Nvapor = length(P)
    loop_plus_index = [2:Nvapor;1]
    loop_plus_index_new = [3:Nvapor+1;1:2]

    Lfilm_start_new = insert!(Lfilm_start,index+1,Linsert/4)
    Lfilm_end_new = insert!(Lfilm_end,index+1,Linsert/4)

    max_index = findmax([Lfilm_start_new[index],Lvaporplug[index] - Lfilm_end_new[index] - Lfilm_start_new[index], Lfilm_end_new[index]])[2]
    if max_index == 1 && Lfilm_start_new[index] > Linsert/2
        splice!(Lfilm_start_new,index,Lfilm_start_new[index]-Linsert/2)
    elseif max_index == 3 && Lfilm_end_new[index] > Linsert/2
        splice!(Lfilm_end_new,index,Lfilm_end_new[index]-Linsert/2)
    elseif max_index == 2 && (Lvaporplug[index] - Lfilm_end_new[index] - Lfilm_start_new[index]) < Linsert/2
        println("pure vapor too short")
    end

    max_index = findmax([Lfilm_start_new[loop_plus_index_new[index]],Lvaporplug[loop_plus_index[index]] - Lfilm_end_new[loop_plus_index_new[index]] - Lfilm_start_new[loop_plus_index_new[index]], Lfilm_end_new[loop_plus_index_new[index]]])[2]
    if max_index == 1 && Lfilm_start_new[loop_plus_index_new[index]] > Linsert/2
        splice!(Lfilm_start_new,loop_plus_index_new[index],Lfilm_start_new[loop_plus_index_new[index]]-Linsert/2)
    elseif max_index == 3 && Lfilm_end_new[loop_plus_index_new[index]] > Linsert/2
        splice!(Lfilm_end_new,loop_plus_index_new[index],Lfilm_end_new[loop_plus_index_new[index]]-Linsert/2)
    elseif max_index == 2 && (Lvaporplug[loop_plus_index[index]] - Lfilm_end_new[loop_plus_index_new[index]] - Lfilm_start_new[loop_plus_index_new[index]]) < Linsert/2
        println("pure vapor too short")
    end

    Lliquid_adjust = 0
    Xpnew = getnewXp(Xp,index,Xvapornew,Lliquid_adjust,L,closedornot)

    Lvaporplug_new = XptoLvaporplug(Xpnew,sys.tube.L,sys.tube.closedornot)
    Astart_new = getδarea(Ac,d,δstart_new)
    Aend_new = getδarea(Ac,d,δend_new)
    Volumevapor_new = getVolumevapor(Ac,Astart_new,Aend_new,Lvaporplug_new,Lfilm_start_new,Lfilm_end_new)

    Pnew_left = DtoP(Mvapor[index]/Volumevapor_new[index])
    Pnew_right = DtoP(Mvapor[loop_plus_index[index]]/Volumevapor_new[loop_plus_index_new[index]])

    # println(P[index-2:index+5])

    # avoid collapse right after boiling
    Pnew = insert!(P,index+1,maximum([Pinsert,Pnew_left,Pnew_right]))
    splice!(Pnew,index,Pnew_left)
    splice!(Pnew,loop_plus_index_new[index],Pnew_right)

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
    sysnew.vapor.δstart = δstart_new
    sysnew.vapor.δend = δend_new
    sysnew.vapor.Lfilm_start = Lfilm_start_new
    sysnew.vapor.Lfilm_end = Lfilm_end_new

    θ_interp_walltoliquid, θ_interp_liquidtowall, H_interp_liquidtowall, P_interp_liquidtowall = sys_interpolation(sysnew)
    sysnew.mapping = Mapping(θ_interp_walltoliquid, θ_interp_liquidtowall, H_interp_liquidtowall, P_interp_liquidtowall)

    # println(Pnew[index-2:index+5])

return sysnew
end


 function getnewθarrays(index,Xp,Xpnew,Xarrays,θarrays,L,closedornot)
    θarraysnew = deepcopy(θarrays)
    Xarraysnew = deepcopy(Xarrays)
    arrayindex = getarrayindex(Xpnew[index][2],Xarrays[index])

    θarraysnewleft = θarrays[index][1:arrayindex]

    θarraysnewright= θarrays[index][arrayindex+1:end]
    insert!(θarraysnewright, 1, θarrays[index][arrayindex])

    splice!(θarraysnew, index)
    insert!(θarraysnew, index,θarraysnewleft)
    insert!(θarraysnew, index+1,θarraysnewright)
end


function getnewXarrays(index,Xp,Xpnew,Xarrays,L,closedornot)
    Xarraysnew = deepcopy(Xarrays)
    arrayindex = getarrayindex(Xpnew[index][2],Xarrays[index])

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

function getnewXp(Xp,index,Xvapornew,Lliquid_adjust,L,closedornot)

    Xpnew = deepcopy(Xp)

    Linsert = mod(Xvapornew[end] - Xvapornew[1],L)

    insertXp1=mod.((Xp[index][1]-Linsert/2+Lliquid_adjust/2,Xvapornew[1]),L)
    insertXp2=mod.((Xvapornew[2],Xp[index][2]+Linsert/2-Lliquid_adjust/2),L)

    splice!(Xpnew, index)
    insert!(Xpnew, index,insertXp1)
    insert!(Xpnew, index+1,insertXp2)

    return Xpnew
end

# simplified non mass conservation
function getnewXp(Xp,index,Xvapornew,L)

    Xpnew = deepcopy(Xp)

    insertXp1=mod.((Xp[index][1],Xvapornew[1]),L)
    insertXp2=mod.((Xvapornew[2],Xp[index][2]),L)

    splice!(Xpnew, index)
    insert!(Xpnew, index,insertXp1)
    insert!(Xpnew, index+1,insertXp2)

    return Xpnew
end


function getnewM(M,index,Minsert,closedornot)

    Mnew = deepcopy(M)

    insert!(Mnew,index+1,Minsert)

    return Mnew
end


function getsuperheat(Xstation,sys)
    Δθ = sys.mapping.θ_interp_walltoliquid(Xstation) - PtoT(sys.mapping.P_interp_liquidtowall(Xstation))
    return Δθ
end


```
    get the index of element of Xarray closest to X.
    closed loop considered
```
function getoneXarrayindex(X,Xarray)
    for i = 1:length(Xarray)
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
    suitable_flag =  true
    index_max = length(p.liquid.Xp)


        index = getinsertindex(p.liquid.Xp,p.wall.Xstations[i],p.tube.L,p.tube.closedornot)
        Lvaporplug = XptoLvaporplug(p.liquid.Xp,p.tube.L,p.tube.closedornot)

        L_vapor_left =  Lvaporplug[index]
        L_vapor_right = (index == index_max) ? Lvaporplug[1] : Lvaporplug[index + 1]

        suitable_flag = (5*p.tube.d < L_vapor_left) && (5*p.tube.d < L_vapor_right) ? true : false


        L_liquid_left =  mod(p.wall.Xstations[i] - p.liquid.Xp[index][1],p.tube.L)
        L_liquid_right = mod(p.liquid.Xp[index][2] - p.wall.Xstations[i],p.tube.L)

                # println(L_liquid_left)
                # println(L_liquid_right)

        if (10*p.tube.d > L_liquid_left) || (10*p.tube.d > L_liquid_right)
            suitable_flag = false
        end

    return suitable_flag
end


#
# end
