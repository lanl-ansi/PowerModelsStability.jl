"add the inverters from the read-in json dictionary"
function add_inverters!(pmd_data::Dict{String,<:Any}, inverter_data::Dict{String,<:Any}; pop_solar::Bool=false)
    for (invInd, inverter) in enumerate(get(inverter_data, "inverters", []))
        bus_gen_terms = get(pmd_data, "is_kron_reduced", false) ? collect(1:3) : collect(1:4)

        # add bus/bus_lookup
        _PMD.add_bus!(
            pmd_data,
            "inverter_$(invInd)";
            terminals=bus_gen_terms,
            grounded=get(pmd_data, "is_kron_reduced", false) ? Int[] : [4],
            rg=[0.0],
            xg=[0.0],
            mp=inverter["mp"],
            mq=inverter["mq"],
            tau=inverter["tau"],
            inverter_bus=true)

        # add gen
        _PMD.add_generator!(
            pmd_data,
            "invGen_$(invInd)",
            "inverter_$(invInd)",
            Int[1,2,3];
            pg=Vector{Float64}(inverter["pg"]),
            qg=Vector{Float64}(inverter["qg"]),
            pg_ub=Vector{Float64}(inverter["pg_ub"]),
            qg_ub=Vector{Float64}(inverter["qg_ub"]),
            pg_lb=Vector{Float64}(inverter["pg_lb"]),
            qg_lb=Vector{Float64}(inverter["qg_lb"]),
            vg=Vector{Float64}(inverter["vg"]),
            phases=3
        )
        # TODO bug in PMD?
        pmd_data["generator"]["invGen_$(invInd)"]["connections"] = bus_gen_terms

        # add connecting branch
        _PMD.add_line!(
            pmd_data,
            "invLine_$(invInd)",
            inverter["busID"],
            "inverter_$(invInd)",
            Int[1,2,3],
            Int[1,2,3],
            rs=inverter["r"],
            xs=inverter["x"]
        )

    end

    if pop_solar
        # pop out the solar since it is replaced by the inverters
        delete!(pmd_data,"solar")
    end
end
