@testset "test the two bus system" begin
    
    println(">>>>>>> Testing 2-bus system <<<<<<<<")

    # load the data
    filePath = "../test/data/case2_diag.dss"
    inverter_file = "../test/data/case2_inverters.json"

    # solve the opf problem
    mpData = PMS.parse_file(filePath, inverter_file)
    mpData["settings"]["sbase_default"] = 1e5

    # obtain the opf solution and the opf model
    mpData_math = PMS.transform_data_model(mpData)
    opfSol = PMS.solve_mc_opf(mpData_math, PMD.ACPUPowerModel, ipopt_solver)
    pm = PMS.instantiate_mc_model(mpData_math, PMD.ACPUPowerModel, PMD.build_mc_opf)

    # change the virtual bus to be an inverter bus for test
    omega0 = mpData["omega0"]
    rN = mpData["rN"]

    Atot = PMS.get_global_stability_matrix(mpData_math, opfSol, omega0, rN)
    eigValList = LA.eigvals(Atot)
    eigVectorList = LA.eigvecs(Atot)
    statusTemp = true
    vioList = []
    for eigInd in 1:length(eigValList)
        eig = eigValList[eigInd]
        if eig.re > 0
            statusTemp = false
            push!(vioList, eigVectorList[eigInd,:])
        end
    end

    @test statusTemp
    @test isempty(vioList)

    if !statusTemp
        Amg = PMS.obtainGlobal_var(mpData_math,pm,rN,omega0)
        PMS.constraint_stability(pm, 0, vioList, Amg)
    end
end
