""
function parse_file(dss_file::String, inverter_file::String; pop_solar::Bool=false, kwargs...)
    pmd_data = _PMD.parse_file(dss_file; kwargs...)

    inverter_data = PMS.parse_json(inverter_file)

    PMS.add_inverters!(pmd_data, inverter_data; pop_solar=pop_solar)

    for (k, v) in inverter_data
        if k != "inverters"
            if haskey(pmd_data, k)
                Memento.warn(_LOGGER, "'$k' from inverter data already exists in network data, skipping")
            else
                pmd_data[k] = v
            end
        end
    end

    return pmd_data
end
