"Obtain the branches that are between two buses and buses with inverters"
function preproc(mpData::Dict{String,<:Any})::Tuple
    brList = []
    invList = []
    invLine = Dict()
    invConnected = Dict()
    for k in keys(mpData["branch"])
        t_bus = mpData["branch"][k]["t_bus"]
        f_bus = mpData["branch"][k]["f_bus"]
        if !(mpData["bus"]["$(t_bus)"]["inverter_bus"]) && !(mpData["bus"]["$(f_bus)"]["inverter_bus"])
            push!(brList,k)
        elseif mpData["bus"]["$(t_bus)"]["inverter_bus"]
            push!(invList,"$(f_bus)")
            if "$(f_bus)" in keys(invConnected)
                # there has already been an inverter connected to f_bus
                push!(invConnected["$(f_bus)"],"$(t_bus)")
            else
                invConnected["$(f_bus)"] = ["$(t_bus)"]
            end
            invLine["$(t_bus)"] = k
        else
            push!(invList,"$(t_bus)")
            if "$(t_bus)" in keys(invConnected)
                # there has already been an inverter connected to f_bus
                push!(invConnected["$(t_bus)"],"$(f_bus)")
            else
                invConnected["$(t_bus)"] = ["$(f_bus)"]
            end
            invLine["$(f_bus)"] = k
        end
    end

    # obtain the list of buses with loads, assuming balanced loads
    loadList = Dict()
    vnomList = Dict()
    loadConnections = Dict()
    for i in keys(mpData["load"])
        # println(i," ",mpData["load"][i]["load_bus"])
        if !("$(mpData["load"][i]["load_bus"])" in keys(loadList))
            loadList["$(mpData["load"][i]["load_bus"])"] = [[0.0,0.0,0.0],[0.0,0.0,0.0]]
            loadList["$(mpData["load"][i]["load_bus"])"][1][mpData["load"][i]["connections"]] += mpData["load"][i]["pd"]
            loadList["$(mpData["load"][i]["load_bus"])"][2][mpData["load"][i]["connections"]] += mpData["load"][i]["qd"]
            vnomList["$(mpData["load"][i]["load_bus"])"] = [0.0,0.0,0.0]
            vnomList["$(mpData["load"][i]["load_bus"])"][mpData["load"][i]["connections"]] += (mpData["load"][i]["pd"] .!= 0.0)*mpData["load"][i]["vnom_kv"]
            loadConnections["$(mpData["load"][i]["load_bus"])"] = deepcopy(mpData["load"][i]["connections"])
        else
            loadList["$(mpData["load"][i]["load_bus"])"][1][mpData["load"][i]["connections"]] += mpData["load"][i]["pd"]
            vnomList["$(mpData["load"][i]["load_bus"])"][mpData["load"][i]["connections"]] += (mpData["load"][i]["pd"] .!= 0.0)*mpData["load"][i]["vnom_kv"]
            for ic in mpData["load"][i]["connections"]
                if !(ic in loadConnections["$(mpData["load"][i]["load_bus"])"])
                    append!(loadConnections["$(mpData["load"][i]["load_bus"])"],ic)
                end
            end
        end
    end

    # get the bus list
    busList = []
    for i in keys(mpData["bus"])
        if !(mpData["bus"][i]["inverter_bus"])
            push!(busList, i)
        end
    end

    return busList, brList, invList, invConnected, invLine, loadList, vnomList, loadConnections
end

"Obtain the load parameters from the model data"
function procLoad(mpData, loadList, vnomList, omega0, loadConnections)
    θList = [0,2π/3,-2π/3]
    load_R = Dict()
    load_X = Dict()
    load_L = Dict()
    # Get the load resistance and reactance
    for i in keys(loadList)
        # i is the bus which the load connects to
        load_R[i] = zeros(3,3)
        load_X[i] = zeros(3,3)
        for j in 1:3
            if j in loadConnections[i]
                if vnomList[i][j] != 0
                    loadRX = vnomList[i][j]^2/(loadList[i][1][j]-im*loadList[i][2][j]);
                    load_R[i][j,j] = loadRX.re
                    load_X[i][j,j] = loadRX.im
                end
            end
        end
        load_L[i] = load_X[i]./(omega0)
    end
    return load_L,load_R,load_X
end
