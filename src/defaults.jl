function default_convection_velocity(v,t,base_cache,phys_params)
    fill!(v,0)
end

const DEFAULT_CONVECTION_VELOCITY_FUNCTION = default_convection_velocity
const DEFAULT_BACKGROUND_TEMPERATURE = 0.0

function get_convection_velocity_function(phys_params::Dict)
    return get(phys_params,"convection velocity model",DEFAULT_CONVECTION_VELOCITY_FUNCTION)
end

function get_heating_model(forcing::Dict)
	return get(forcing,"heating models",nothing)
end

function get_heating_model(forcing::Nothing)
	return get_heating_model(Dict())
end

function get_background_temp(phys_params::Dict)
    return get(phys_params,"background temperature",DEFAULT_BACKGROUND_TEMPERATURE)
end