# read json files for the inverters
"Parses json from iostream or string"
function parse_json(io::Union{IO,String})
    inverter_data = JSON.parsefile(io);
    for i in 1:length(inverter_data["inverters"])
        rtemp = inverter_data["inverters"][i]["r"];
        xtemp = inverter_data["inverters"][i]["x"];
        inverter_data["inverters"][i]["r"] = reshape(hcat(rtemp...), (length(rtemp[1]), length(rtemp)));
        inverter_data["inverters"][i]["x"] = reshape(hcat(xtemp...), (length(xtemp[1]), length(xtemp)));
    end

    return inverter_data;
end

function add_inverters(inverter_data, mpData)
    # add the inverters from the read-in json dictionary
    invInd = 0;
    for item in inverter_data["inverters"]
        invInd += 1;
        # add bus/bus_lookup
        add_bus!(mpData, "inverter_$(invInd)", terminals = [1,2,3,4], grounded = [4], rg = [0.0], xg = [0.0],
            mp = item["mp"], mq = item["mq"], tau = item["tau"], bus_type = 4);

        # add gen
        add_generator!(mpData, "invGen_$(invInd)", "inverter_$(invInd)", Vector{Int}([1,2,3]),
            pg = Array{Float64,1}(item["pg"]), qg = Array{Float64,1}(item["qg"]),
            pg_ub = Array{Float64,1}(item["pg_ub"]), qg_ub = Array{Float64,1}(item["qg_ub"]),
            pg_lb = Array{Float64,1}(item["pg_lb"]), qg_lb = Array{Float64,1}(item["qg_lb"]),
            vg = Array{Float64,1}(item["vg"]), phases = 3);
        mpData["generator"]["invGen_$(invInd)"]["connections"] = Vector{Int}([1,2,3,4]);

        # add connecting branch
        add_line!(mpData, "invLine_$(invInd)", item["busID"], "inverter_$(invInd)",
            Vector{Int}([1,2,3]), Vector{Int}([1,2,3]), rs = item["r"], xs = item["x"]);

    end
    return mpData;
end

function transform_data_model(mpData)
    mpData_math = PowerModelsDistribution.transform_data_model(mpData);
    data_fields = ["mp","mq","omega0","tau"];
    # add necessary data field to the math model as well
    for i in keys(mpData_math["bus"])
        if mpData_math["bus"][i]["name"] in keys(mpData["bus"])
            for j in keys(mpData["bus"][mpData_math["bus"][i]["name"]])
                # input all keys that
                if j in data_fields
                    mpData_math["bus"][i][j] = mpData["bus"][mpData_math["bus"][i]["name"]][j];
                end
            end
            # change the inverter type
            if "bus_type" in keys(mpData["bus"][mpData_math["bus"][i]["name"]])
                if mpData["bus"][mpData_math["bus"][i]["name"]]["bus_type"] == 4
                    mpData_math["bus"][i]["bus_type"] = 4;
                end
            end
        end
    end
    return mpData_math;
end
