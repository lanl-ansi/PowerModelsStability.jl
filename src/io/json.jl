# read json files for the inverters
"Parses json from iostream or string"
function parse_json(io::Union{IO,String})::Dict{String,Any}
    inverter_data = JSON.parsefile(io)
    for i in 1:length(inverter_data["inverters"])
        rtemp = inverter_data["inverters"][i]["r"]
        xtemp = inverter_data["inverters"][i]["x"]
        inverter_data["inverters"][i]["r"] = reshape(hcat(rtemp...), (length(rtemp[1]), length(rtemp)))
        inverter_data["inverters"][i]["x"] = reshape(hcat(xtemp...), (length(xtemp[1]), length(xtemp)))
    end

    return inverter_data
end

"add the inverters from the read-in json dictionary"
function add_inverters!(mpData, inverter_data, pop_solar = false)
    invInd = 0
    for item in get(inverter_data, "inverters", [])
        invInd += 1
        # add bus/bus_lookup
        _PMD.add_bus!(mpData, "inverter_$(invInd)", terminals = [1,2,3,4], grounded = [4], rg = [0.0], xg = [0.0],
            mp = item["mp"], mq = item["mq"], tau = item["tau"], inverter_bus = true)

        # add gen
        _PMD.add_generator!(mpData, "invGen_$(invInd)", "inverter_$(invInd)", Vector{Int}([1,2,3]),
            pg = Array{Float64,1}(item["pg"]), qg = Array{Float64,1}(item["qg"]),
            pg_ub = Array{Float64,1}(item["pg_ub"]), qg_ub = Array{Float64,1}(item["qg_ub"]),
            pg_lb = Array{Float64,1}(item["pg_lb"]), qg_lb = Array{Float64,1}(item["qg_lb"]),
            vg = Array{Float64,1}(item["vg"]), phases = 3)
        mpData["generator"]["invGen_$(invInd)"]["connections"] = Vector{Int}([1,2,3,4])

        # add connecting branch
        _PMD.add_line!(mpData, "invLine_$(invInd)", item["busID"], "inverter_$(invInd)",
            Vector{Int}([1,2,3]), Vector{Int}([1,2,3]), rs = item["r"], xs = item["x"])

    end
    if pop_solar
        # pop out the solar since it is replaced by the inverters
        delete!(mpData,"solar")
    end
end
