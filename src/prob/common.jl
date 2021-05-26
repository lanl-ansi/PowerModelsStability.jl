instantiate_mc_model(
    data::Dict{String,<:Any},
    model_type::Type,
    build_method::Function;
    eng2math_extensions::Vector{<:Function}=Function[],
    kwargs...) = _PMD.instantiate_mc_model(
        data,
        model_type,
        build_method;
        eng2math_extensions=[_eng2math_inverter_bus!, eng2math_extensions...],
        eng2math_passthrough=_pms_eng2math_passthrough,
        kwargs...)

solve_mc_opf(
    data::Dict{String,<:Any},
    model_type::Type,
    solver;
    eng2math_extensions::Vector{<:Function}=Function[],
    kwargs...) = _PMD.solve_mc_opf(
        data,
        model_type,
        solver;
        eng2math_extensions=[_eng2math_inverter_bus!, eng2math_extensions...],
        eng2math_passthrough=_pms_eng2math_passthrough,
        kwargs...)
