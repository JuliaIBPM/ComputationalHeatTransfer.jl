temperature(T,σ,x,sys::ILMSystem,t) = T .+ get_background_temp(sys.phys_params)
@snapshotoutput temperature