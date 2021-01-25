# read json files for the inverters
"Parses json from iostream or string"
function parse_json(io::Union{IO,String})
    inverter_data = JSON.parsefile(io);

    return inverter_data;
end

function add_inverters(inverter_data, mpData)
    # add the inverters from the read-in json dictionary
    invInd = 0;
    for item in inverter_data["inverters"]
        invInd += 1;
        # add bus/bus_lookup
        add_bus!(mpData, "inverter_$(invInd)", terminals = [1,2,3,4], grounded = [4], rg = [0.0], xg = [0.0],
            mp = item["mp"], mq = item["mq"], tau = item["tau"], omega0 = item["omega0"], bus_type = 4);

        # add gen
        add_generator!(mpData, "invGen_$(invInd)", "inverter_$(invInd)", Vector{Int}([1,2,3]),
            pg = Array{Float64,1}(item["pg"]), qg = Array{Float64,1}(item["qg"]),
            pg_ub = Array{Float64,1}(item["pg_ub"]), qg_ub = Array{Float64,1}(item["qg_ub"]),
            pg_lb = Array{Float64,1}(item["pg_lb"]), qg_lb = Array{Float64,1}(item["qg_lb"]),
            vg = Array{Float64,1}(item["vg"]), phases = 3);
        mpData["generator"]["invGen_$(invInd)"]["connections"] = Vector{Int}([1,2,3,4]);

        # add connecting branch
        r = zeros(3,3);
        x = zeros(3,3);
        for i in 1:3
            for j in 1:3
                r[i,j] = item["r"][i][j];
                x[i,j] = item["x"][i][j];
            end
        end

        add_line!(mpData, "invLine_$(invInd)", item["busID"], "inverter_$(invInd)",
            Vector{Int}([1,2,3]), Vector{Int}([1,2,3]), rs = r, xs = x);

    end
    return mpData;
end
