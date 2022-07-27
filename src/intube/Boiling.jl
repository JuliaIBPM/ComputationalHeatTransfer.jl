# module Boiling
#
# using ..Systems,..Tools
# using CSV

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
    # println(length(p.liquid.Xp))

    Δθ_array = getsuperheat.(p.wall.Xstations,[p])
    superheat_flag = Δθ_array .> Δθthreshold

    b_count = 0;
    for i = 1:length(p.wall.Xstations)

        if ifamong(p.wall.Xstations[i], p.liquid.Xp, p.tube.L) && suitable_for_boiling(p,i) && superheat_flag[i]

            Δθ = getsuperheat(p.wall.Xstations[i],p)
            # println(i)
            # if Δθ > Δθthreshold
  
                push!(boil_hist,[i,integrator.t]);
                b_count += 1;
# ```modified here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!```
#                 θinsert = p.mapping.θ_interp_walltoliquid.(Xstations[i])
#                 Pinsert = nondi_TtoP.(θinsert)

                Pinsert = p.mapping.P_interp_liquidtowall(Xstations[i])
                # θinsert = nondi_PtoT.(Pinsert)
                # θinsert = PtoT.(Pinsert)


                p = nucleateboiling(p,(p.wall.Xstations[i]-2p.tube.d,p.wall.Xstations[i]+2p.tube.d),Pinsert) # P need to be given from energy equation
            # end
        end


        # println(p.vapor.P)

    end

    # print("boil number=",b_count)
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

    # M = PtoD.(p.vapor.P) .* Lvaporplug .* Ac .* ((d .- 2 .* δ) ./ d) .^2


    unew=[XMδLtovec(p.liquid.Xp,p.liquid.dXdt,M,p.vapor.δstart,δend,Lfilm_start,Lfilm_end);liquidθtovec(p.liquid.θarrays)];

#     set_u!(integrator,  unew)
    resize!(integrator.u,length(unew))
    integrator.u = deepcopy(unew)


    # println(length(p.liquid.Xp))
end

function nucleateboiling(sys,Xvapornew,Pinsert)
    ρ = deepcopy(sys.liquid.ρ)
    d = deepcopy(sys.tube.d)
    Xp = deepcopy(sys.liquid.Xp)
    dXdt = deepcopy(sys.liquid.dXdt)
    # δ = deepcopy(sys.vapor.δ)
    δstart = deepcopy(sys.vapor.δstart)
    δend = deepcopy(sys.vapor.δend)
    δfilm_deposit = deepcopy(sys.vapor.δfilm_deposit)
    P = deepcopy(sys.vapor.P)
    L = sys.tube.L
    Ac = sys.tube.Ac
    # γ = sys.vapor.γ
    Xarrays = sys.liquid.Xarrays
    θarrays = sys.liquid.θarrays
    closedornot = sys.tube.closedornot
    Lfilm_start = deepcopy(sys.vapor.Lfilm_start)
    Lfilm_end = deepcopy(sys.vapor.Lfilm_end)

    # δarea_start = Ac .* (1 .- ((d .- 2*δstart) ./ d) .^ 2);
    # δarea_end = Ac .* (1 .- ((d .- 2*δend) ./ d) .^ 2);


    Lvaporplug =    XptoLvaporplug(Xp,sys.tube.L,sys.tube.closedornot)

    # volume_vapor = Lvaporplug .* Ac - Lfilm_start .* δarea_start - Lfilm_end .* δarea_end
    # M = PtoD.(p.vapor.P) .* volume_vapor
    # M = PtoD.(P) .* Lvaporplug .* Ac  .* ((d .- 2 .* δ) ./ d) .^2
    # M = nondi_PtoD.(P) .* Lvaporplug
    # M = P.^(1/γ).* Lvaporplug

    index = getinsertindex(Xp,(Xvapornew[2]+Xvapornew[1])/2,sys.tube.L,sys.tube.closedornot)




    # Xvapornew = (1.50,1.60)
    # Pinsert = 1.5
    # ρinsert = nondi_PtoD(Pinsert)
    ρinsert = PtoD(Pinsert)

    Linsert = mod(Xvapornew[end] - Xvapornew[1],L)
    # Mvaporinsert = ρinsert .* Linsert .* Ac  .* ((d .- 2 .* δ[index]) ./ d) .^2



    """let's do constant film thickness for now!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"""
    δstart_new = insert!(δstart,index+1,δfilm_deposit)
    δend_new = insert!(δend,index+1,δfilm_deposit)

    # L0threshold = 4e-4

    Nvapor = length(P)
    loop_plus_index = [2:Nvapor;1]
    loop_plus_index_new = [3:Nvapor+1;1:2]

    # println(length(P))

    # if  Lvaporplug[index] - (Lfilm_start[index] + Lfilm_end[index]) < Linsert/2
    #     (Lfilm_start[index] > Lfilm_end[index]) ? Lfilm_start[index] -= Linsert/2 : Lfilm_end[index] -= Linsert/2
    # end
    # if Lvaporplug[loop_plus_index[index]] - (Lfilm_start[loop_plus_index[index]] + Lfilm_end[loop_plus_index[index]]) > L0threshold
    # (   Lfilm_start[loop_plus_index[index]] > Lfilm_end[loop_plus_index[index]] )  ? Lfilm_start[loop_plus_index[index]] -= Linsert/2 : Lfilm_end[loop_plus_index[index]] -= Linsert/2
    # end
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
    #    println("hahaha")
    max_index = findmax([Lfilm_start_new[loop_plus_index_new[index]],Lvaporplug[loop_plus_index[index]] - Lfilm_end_new[loop_plus_index_new[index]] - Lfilm_start_new[loop_plus_index_new[index]], Lfilm_end_new[loop_plus_index_new[index]]])[2]
    if max_index == 1 && Lfilm_start_new[loop_plus_index_new[index]] > Linsert/2
        splice!(Lfilm_start_new,loop_plus_index_new[index],Lfilm_start_new[loop_plus_index_new[index]]-Linsert/2)
    elseif max_index == 3 && Lfilm_end_new[loop_plus_index_new[index]] > Linsert/2
        # println("hahaha")
        # println(index)
        # println(Lfilm_end_new)
        # # println(Lfilm_end_new[loop_plus_index_new[index]])
        splice!(Lfilm_end_new,loop_plus_index_new[index],Lfilm_end_new[loop_plus_index_new[index]]-Linsert/2)
        # println(Lfilm_end_new[loop_plus_index_new[index]])
    elseif max_index == 2 && (Lvaporplug[loop_plus_index[index]] - Lfilm_end_new[loop_plus_index_new[index]] - Lfilm_start_new[loop_plus_index_new[index]]) < Linsert/2
        println("pure vapor too short")
    end
    # δnew = insert!(δ,index+1,δ[index])
    # δnew = getnewδ(δ,index,Xvapornew,ρ,d,Minsert,closedornot) # mass conservation
    # δarea = getδarea(Ac,d,δ[index])
    # Mfilminsert = ρ.*δarea.*Linsert

    # Minsert = Mfilminsert + Mvaporinsert

    # left_index = index
    # right_index = index < length(sys.vapor.δstart) ? index+1 : 1

    # Mperlength_left = getMperlength(sys,left_index)
    # Mperlength_right = getMperlength(sys,right_index)

    # Lliquid_adjust = (Minsert - Mperlength_left*Linsert/2 - Mperlength_right*Linsert/2) / (ρ*Ac - Mperlength_left/2 - Mperlength_right/2)
    # # Xpnew =
    # Mnew = getnewM(M,index,Minsert,closedornot)
    #
    # Lvaporplugnew = XptoLvaporplug(Xpnew,L,closedornot)
    # Pnew = nondi_DtoP.(Mnew./Lvaporplugnew)

    Lliquid_adjust = 0
    Xpnew = getnewXp(Xp,index,Xvapornew,Lliquid_adjust,L,closedornot)
    # simplified Xpnew
    # Xpnew = getnewXp(Xp,index,Xvapornew,L)
    Pnew = insert!(P,index+1,Pinsert)


    # # modify
    # Pnew[index] = Pnew[index] + 10000
    # Pnew[loop_plus_index[index]] =  Pnew[loop_plus_index[index]] + 20000
    # Pnew[loop_plusplus_index[index]] = Pnew[loop_plusplus_index[index]] + 10000

    Xarraysnew = getnewXarrays(index,Xp,Xpnew,Xarrays,L,closedornot)
    θarraysnew = getnewθarrays(index,Xp,Xpnew,Xarrays,θarrays,L,closedornot)



# only for open loop!
    dXdtnew = deepcopy(dXdt) # momentum conservation
    insert!(dXdtnew,index+1,dXdtnew[index])

    # dXdt_boil_adjust_plus = dXdtnew[right_index] .+ (0.5,0.5)
    # dXdtnew[right_index] = dXdt_boil_adjust_plus
    # dXdt_boil_adjust_minus = dXdtnew[index] .- (0.5,0.5)
    # dXdtnew[index] = dXdt_boil_adjust_minus

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

    
    # sysnew.vapor.δ = δnew

    # walltoliquid,liquidtowall = constructmapping(sysnew.liquid.Xarrays ,sysnew.wall.Xarray, sysnew.tube.closedornot, sysnew.tube.L)
    # print(Xarraysnew[5],"\n",θarraysnew[5],"\n")
    # sysnew.mapping = Mapping(walltoliquid,liquidtowall)


    θ_interp_walltoliquid, θ_interp_liquidtowall, H_interp_liquidtowall, P_interp_liquidtowall = sys_interpolation(sysnew)
    sysnew.mapping = Mapping(θ_interp_walltoliquid, θ_interp_liquidtowall, H_interp_liquidtowall, P_interp_liquidtowall)

return sysnew
end


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
    # suitable_flag =  false
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
