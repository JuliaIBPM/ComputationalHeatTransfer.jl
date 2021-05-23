# # module Postprocessing
#
# export soltoResult,soltoMatrxResult
#
# using ..Systems,..Tools
#
# """
#     reserved for future use
# """
#
# function soltoResult(sol,sys0) #only good for one calculation per time point, not good for onec calculation for all time
#
#     γ = sys0.liquid.γ
#     numofliquidslug =  Integer( (size(sol)[1]-1)/5  )
#
#     MatrxXp=sol[1:2*numofliquidslug,:]
#     MatrxdXdt=sol[2*numofliquidslug+1:4*numofliquidslug,:]
#     M=sol[4*numofliquidslug+1:5*numofliquidslug+1,:]
#
#     Xp = Array{Tuple{Float64, Float64}}((undef), Integer(size(MatrxXp)[1]/2), size(MatrxXp)[2])
#     dXdt = Array{Tuple{Float64, Float64}}((undef), size(Xp))
#     Lvaporplug = zeros(size(M))
#
#     # transfer from 2D array to array of tuple for Xp and dXdt
#     for i in 1:size(Xp)[1]
#         for j in 1:size(Xp)[2]
#             Xp[i,j]=((MatrxXp[2*i-1,j]),(MatrxXp[2*i,j]))
#             dXdt[i,j]=((MatrxdXdt[2*i-1,j]),(MatrxdXdt[2*i,j]))
#         end
#     end
#
#     for j in 1:size(M)[2]
#         Lvaporplug[:,j]=XptoLvaporplug(Xp[:,j], sys0.tube.L, sys0.tube.closedornot)
#     end
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
#
#     # θ = P.^((γ-1)/γ)
#
#
#     # convert a matrix to array of array
#     Xp=mapslices(x->[x], MatrxXp, dims=2)[:]
#     dXdt=mapslices(x->[x], MatrxdXdt, dims=2)[:]
#     M=mapslices(x->[x], M, dims=2)[:]
#     P=mapslices(x->[x], P, dims=2)[:]
#     θ=mapslices(x->[x], θ, dims=2)[:]
#
#     # if the input sol is merely a state vector u Array{Float64,1} rather than ODESolution{}
#     # the result then does not contain t (we set the t to be -1.0)
#     if typeof(sol) == Array{Float64,1}
#         return PHPResult(-1.0, vec(hcat(Xp...)),vec(hcat(dXdt...)),vec(hcat(P...)),vec(hcat(θ...)),vec(hcat(M...)))
#     elseif typeof(sol) == Array{Float64,2}
#         return PHPResult(-1.0,Xp,dXdt,P,θ,M)
#     else
#         return PHPResult(sol.t,Xp,dXdt,P,θ,M)
#     end
#     # result=typeof(sol) == Array{Float64,1} ? PHPResult(-1.0, vec(hcat(Xp...)),vec(hcat(dXdt...)),vec(hcat(P...)),vec(hcat(θ...)),vec(hcat(M...))) : PHPResult(sol.t,Xp,dXdt,P,θ,M)
# end
#
#
# function soltoMatrxResult(sol,sys0) #only good for one calculation per time point, not good for onec calculation for all time
#
#     indexes = Int64[]
#     θliquidrec = Array[]
#
#     for i = 1:length(sol[1])
#         if abs(sol[1][i]+1e10) <= 10^(-1)
#             push!(indexes,i)
#         end
#     end
#
#
#     γ = sys0.liquid.γ
#     # numofliquidslug =  Integer( (size(sol)[1]-1)/5  )
#     # numofliquidslug =  Integer( (indexes[1]-2)/5  )
#     numofliquidslug =  Integer( (indexes[1]-3)/6  )
#
#     MatrxXp=sol[1:2*numofliquidslug,:]
#     MatrxdXdt=sol[2*numofliquidslug+1:4*numofliquidslug,:]
#     MatrxM=sol[4*numofliquidslug+1:5*numofliquidslug+1,:]
#     Matrxδ=sol[5*numofliquidslug+2:6*numofliquidslug+2,:]
#
#     return MatrxXp, MatrxdXdt, MatrxM, Matrxδ
# end
#
# # end
