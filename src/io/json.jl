"Parses json from iostream or string"
function parse_json(io::Union{IO,String})::Dict{String,Any}
    inverter_data = JSON.parsefile(io)

    PMS._fix_inverter_matrices!(inverter_data)

    return inverter_data
end


""
function _fix_inverter_matrices!(inverter_data::Dict{String,<:Any})
    if haskey(inverter_data, "inverters")
        for inverter in inverter_data["inverters"]
            inverter["r"] = Matrix{Float64}(hcat(inverter["r"]...))
            inverter["x"] = Matrix{Float64}(hcat(inverter["x"]...))
        end
    end
end
