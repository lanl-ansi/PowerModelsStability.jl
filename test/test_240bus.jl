@testset "test the Iowa 240 bus system" begin
    # load the data
    filePath = "./data/iowa240/Master_hse_der.dss"
    inverter_data = parse_json("./data/iowa240/case240_inverters.json")

    # solve the opf problem
    mpData = PMD.parse_file(filePath)
    add_inverters!(mpData, inverter_data, true)
    inverter_bus_list = [iKey for iKey in keys(mpData["bus"]) if occursin("inverter",iKey)]

    # obtain the opf solution and the opf model
    apply_voltage_angle_bounds!(mpData,1);
    opfSol,mpData_math = run_mc_opf(mpData, PMD.ACRPowerModel, ipopt_solver; solution_processors=[PMD.sol_data_model!])
    pm = PMD.instantiate_mc_model(mpData_math, PMD.ACRPowerModel, PMD.build_mc_opf;
        ref_extensions=[PMD.ref_add_arcs_transformer!, PMD.ref_add_connections!],solution_processors=[PMD.sol_data_model!])

    # change the virtual bus to be an inverter bus for test
    omega0 = inverter_data["omega0"]
    rN = inverter_data["rN"]
    invData = invData_proc(mpData_math, omega0)

    Atot = obtainGlobal_multi(mpData_math,opfSol,rN,omega0,invData)
    eigValList, eigVectorList = eigen(Atot)
    statusTemp = true
    vioList = []
    for eigInd in 1:length(eigValList)
        eig = eigValList[eigInd]
        if eig.re > 0
            statusTemp = false
            push!(vioList,eigVectorList[:,eigInd])
        end
    end

    @test statusTemp
    @test isempty(vioList)

    if !statusTemp
        Amg = obtainGlobal_var(mpData_math, pm, rN, omega0)
        constraint_stability(pm, 0, vioList, Amg)
    end
end
opfSol = _PM.optimize_model!(pm, optimizer = ipopt_solver, solution_processors = [PMD.sol_data_model!])
