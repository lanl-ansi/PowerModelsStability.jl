# read json files for the inverters
"Parses json from iostream or string"
function parse_json(io::Union{IO,String})
    inverter_data = JSON.parsefile(io);

    return inverter_data;
end

function add_inverters(inverter_data)
    # add the inverters from the read-in json dictionary
    for item in inverter_data
        # add bus/bus_lookup

        # add gen

        # add connecting branch

        # add map
    end
end
