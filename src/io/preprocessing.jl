function preproc(mpData)
    # obtain the branches that are between two buses and buses with inverters
    brList = [];
    invList = [];
    invLine = Dict();
    invConnected = Dict();
    for k in keys(mpData["line"])
        t_bus = mpData["line"][k]["t_bus"];
        f_bus = mpData["line"][k]["f_bus"];
        if !("bus_type" in keys(mpData["bus"]["$(t_bus)"])) && !("bus_type" in keys(mpData["bus"]["$(f_bus)"]))
            push!(brList,k);
        elseif ("bus_type" in keys(mpData["bus"]["$(t_bus)"]))
            if mpData["bus"]["$(t_bus)"]["bus_type"] == 4
                push!(invList,"$(f_bus)");
                if "$(f_bus)" in keys(invConnected)
                    # there has already been an inverter connected to f_bus
                    push!(invConnected["$(f_bus)"],"$(t_bus)");
                else
                    invConnected["$(f_bus)"] = ["$(t_bus)"];
                end
                invLine["$(t_bus)"] = k;
            end
        else
            push!(invList,"$(t_bus)");
            if "$(f_bus)" in keys(invConnected)
                # there has already been an inverter connected to f_bus
                push!(invConnected["$(t_bus)"],"$(f_bus)");
            else
                invConnected["$(t_bus)"] = ["$(f_bus)"];
            end
            invLine["$(f_bus)"] = k;
        end
    end

    # obtain the list of buses with loads, assuming balanced loads
    loadList = Dict();
    vnomList = Dict();
    for i in keys(mpData["load"])
        println(mpData["load"][i]["load_bus"]);
        if !("$(mpData["load"][i]["load_bus"])" in keys(loadList))
            loadList["$(mpData["load"][i]["load_bus"])"] = [[0.0,0.0,0.0],[0.0,0.0,0.0]];
            loadList["$(mpData["load"][i]["load_bus"])"][1][mpData["load"][i]["connections"]] += mpData["load"][i]["pd"];
            loadList["$(mpData["load"][i]["load_bus"])"][2][mpData["load"][i]["connections"]] += mpData["load"][i]["qd"];
            vnomList["$(mpData["load"][i]["load_bus"])"] = [0.0,0.0,0.0];
            vnomList["$(mpData["load"][i]["load_bus"])"][mpData["load"][i]["connections"]] += (mpData["load"][i]["pd"] .!= 0.0)*mpData["load"][i]["vnom_kv"];
        else
            loadList["$(mpData["load"][i]["load_bus"])"][1][mpData["load"][i]["connections"]] += mpData["load"][i]["pd"];
            vnomList["$(mpData["load"][i]["load_bus"])"][mpData["load"][i]["connections"]] += (mpData["load"][i]["pd"] .!= 0.0)*mpData["load"][i]["vnom_kv"];
        end
    end

    # get the bus list
    busList = [];
    for i in keys(mpData["bus"])
        if mpData["bus"][i]["bus_type"] != 4
            push!(busList, i);
        end
    end

    return busList, brList, invList, invConnected, invLine, loadList, vnomList;
end

function procLoad(mpData, loadList, vnomList, ω0)
    θList = [0,2π/3,-2π/3];
    load_R = Dict();
    load_X = Dict();
    load_L = Dict();
    # Get the load resistance and reactance
    for i in keys(loadList)
        load_R[i] = zeros(3,3);
        load_X[i] = zeros(3,3);
        for j in 1:3
            load_R[i][j,j] = loadList[i][1][j]/vnomList[i][j]^2;
            load_X[i][j,j] = loadList[i][2][j]/((vnomList[i][j]*cos(θList[j]))^2 - (vnomList[i][j]*sin(θList[j]))^2);
        end
        load_L[i] = load_X[i]./(ω0);
    end
    return load_L,load_R,load_X;
end
