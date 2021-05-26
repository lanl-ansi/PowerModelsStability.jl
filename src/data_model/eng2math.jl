""
function _eng2math_inverter_bus!(data_math::Dict{String,<:Any}, data_eng::Dict{String,<:Any})
    if haskey(data_math, "bus")
        for (_,bus) in data_math["bus"]
            if !haskey(bus, "inverter_bus")
                bus["inverter_bus"] = false
            end
        end
    end
end


const _pms_eng2math_passthrough = Dict{String,Vector{String}}(
    "bus" => String["mp", "mq", "tau", "inverter_bus"]
)


transform_data_model(
    data_eng::Dict{String,<:Any};
    eng2math_extensions::Vector{<:Function}=Function[],
    kwargs...)::Dict{String,Any} = _PMD.transform_data_model(
        data_eng;
        eng2math_extensions=[_eng2math_inverter_bus!, eng2math_extensions...],
        eng2math_passthrough=_pms_eng2math_passthrough,
        kwargs...)
