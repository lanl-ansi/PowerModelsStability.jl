@testset "test the IEEE 13 bus system" begin
    # load the data
    filePath = "./data/ieee13/IEEE13Nodeckt_David_1.dss"
    inverter_data = parse_json("./data/case13_inverters.json")
    inverter_data_math = Dict()
    inverter_data_math["12"] = Dict()
    inverter_data_math["12"]["tau"] = inverter_data["inverters"][1]["tau"]
    inverter_data_math["12"]["mp"] = inverter_data["inverters"][1]["mp"]
    inverter_data_math["12"]["mq"] = inverter_data["inverters"][1]["mq"]

    # solve the opf problem
    mpData = PMD.parse_file(filePath)
    opfSol,mpData_math = run_mc_opf(mpData, PMD.ACRPowerModel, ipopt_solver; solution_processors=[PMD.sol_data_model!])
    mpData_math = convert_gen(mpData_math,inverter_data_math)
    opfSol,mpData_math = run_mc_opf(mpData_math, PMD.ACRPowerModel, ipopt_solver; solution_processors=[PMD.sol_data_model!])

    # obtain the opf solution and the opf model
    # apply_voltage_angle_bounds!(mpData,1);        !!!! This is not feasible
    pm = PMD.instantiate_mc_model(mpData_math, PMD.ACRPowerModel, PMD.build_mc_opf;
        ref_extensions=[PMD.ref_add_arcs_transformer!, PMD.ref_add_connections!])
    opfSol_pm = PowerModels.optimize_model!(pm,optimizer = ipopt_solver;solution_processors=[PMD.sol_data_model!])

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
        Amg = obtainGlobal_var(mpData,pm,omega0,mP,mQ,tau,rN)
        constraint_stability(pm, 0, vioList, Amg)
    end
end
