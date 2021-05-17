"voltage angle bounds"
function apply_voltage_angle_bounds!(eng::Dict{String,<:Any}, vad::Real)
    if haskey(eng, "line")
        for (_,line) in eng["line"]
            line["vad_lb"] = fill(-vad, length(line["f_connections"]))
            line["vad_ub"] = fill( vad, length(line["f_connections"]))
        end
    end
end

"inverter bus voltage balanced assumptions"
function apply_balanced_inverters!(pm, nw::Int, inv_bus_list)
    # force the inverter bus voltage magnitude the same and angle 120 degrees apart
    for inv_bus in inv_bus_list
        idx = pm.data["bus_lookup"][inv_bus]
        vi = _PMD.var(pm, nw, :vi, idx)
        vr = _PMD.var(pm, nw, :vr, idx)
        # use a rectangular formulation
        JuMP.@NLconstraint(pm.model, vi[1]^2 + vr[1]^2 == vi[2]^2 + vr[2]^2)
        JuMP.@NLconstraint(pm.model, vi[1]^2 + vr[1]^2 == vi[3]^2 + vr[3]^2)
        JuMP.@NLconstraint(pm.model, vi[2]*vr[1] - vi[1]*vr[2] == sqrt(3)*(vi[1]*vi[2] + vr[1]* vr[2]))
        JuMP.@NLconstraint(pm.model, vi[1]*vr[3] - vi[3]*vr[1] == sqrt(3)*(vi[1]*vi[3] + vr[1]* vr[3]))
    end
end

"Run the opf solution process and obtain the MATHEMATICAL model"
function run_mc_model(data::Dict{String,<:Any}, model_type::Type, solver, build_mc::Function; ref_extensions::Vector{<:Function}=Vector{Function}([]), make_si=!get(data, "per_unit", false), multinetwork::Bool=false, kwargs...)
    dataOut = Dict()
    if get(data, "data_model", _PMD.MATHEMATICAL) == _PMD.ENGINEERING
        data_math = _PMD.transform_data_model(data; build_multinetwork=multinetwork)
        for i in keys(data_math["bus"])
            try
                busName = data_math["bus"][i]["name"]
                for j in keys(data["bus"][busName])
                    if j in ["mp","mq","tau","inverter_bus"]
                        data_math["bus"][i][j] = data["bus"][busName][j]
                    end
                end
                # if it is not an inverter bus
                if !("inverter_bus" in keys(data_math["bus"][i]))
                    data_math["bus"][i]["inverter_bus"] = false
                end
            catch
                data_math["bus"][i]["inverter_bus"] = false
            end
        end
        result = _PM.run_model(data_math, model_type, solver, build_mc; ref_extensions=[_PMD.ref_add_arcs_transformer!, _PMD.ref_add_connections!, ref_extensions...], multiconductor=true, multinetwork=multinetwork, kwargs...)
        dataOut = deepcopy(data_math)
    elseif get(data, "data_model", _PMD.MATHEMATICAL) == _PMD.MATHEMATICAL
        result = _PM.run_model(data, model_type, solver, build_mc; ref_extensions=[_PMD.ref_add_arcs_transformer!, _PMD.ref_add_connections!, ref_extensions...], multiconductor=true, multinetwork=multinetwork, kwargs...)
        dataOut = deepcopy(data)
    end

    return result, dataOut
end

function run_mc_opf(data::Union{Dict{String,<:Any},String}, model_type::Type, solver; kwargs...)
    return run_mc_model(data, model_type, solver, _PMD.build_mc_opf; kwargs...)
end

"Obtain the MATHEMATICAL model without running the opf solution process"
function transform_data_model(mpData::Dict{String,<:Any})::Dict{String,<:Any}
    apply_voltage_angle_bounds!(mpData,1);
    mpData_math = _PMD.transform_data_model(mpData)
    data_fields = ["mp","mq","inverter_bus","tau"]
    # add necessary data field to the math model as well
    for i in keys(mpData_math["bus"])
        try
            busName = mpData_math["bus"][i]["name"]
            for j in keys(mpData["bus"][busName])
                if j in data_fields
                    mpData_math["bus"][i][j] = mpData["bus"][busName][j]
                end
            end
            # if it is not an inverter bus
            if !("inverter_bus" in keys(mpData_math["bus"][i]))
                mpData_math["bus"][i]["inverter_bus"] = false
            end
        catch
            mpData_math["bus"][i]["inverter_bus"] = false
        end
    end
    return mpData_math
end

"Make all buses that are type 2 or 3 inverter-based"
function convert_gen(mpData_math::Dict{String,<:Any},inverter_data_math::Dict{<:Any,<:Any})::Dict{String,<:Any}
    for i in keys(mpData_math["bus"])
        if mpData_math["bus"][i]["bus_type"] in [2,3]
            mpData_math["bus"][i]["inverter_bus"] = true
            for key in ["mp", "mq", "tau"]
                mpData_math["bus"][i][key] = inverter_data_math[i][key]
            end
        end
    end
    return mpData_math;
end
