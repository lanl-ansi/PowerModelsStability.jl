@testset "test the two bus system" begin
    # load the data
    filePath = "../test/data/case2_diag.dss"
    inverter_file = "../test/data/case2_inverters.json"

    # solve the opf problem
    mpData = parse_file(filePath, inverter_file)

    # obtain the opf solution and the opf model
    mpData_math = transform_data_model(mpData)
    opfSol = solve_mc_opf(mpData_math, PMD.ACPUPowerModel, ipopt_solver)
    pm = instantiate_mc_model(mpData_math, PMD.ACPUPowerModel, PMD.build_mc_opf)

    # change the virtual bus to be an inverter bus for test
    omega0 = mpData["omega0"]
    rN = mpData["rN"]

    Atot = obtainGlobal_multi(mpData_math,opfSol,omega0,rN)
    eigValList = eigvals(Atot)
    eigVectorList = eigvecs(Atot)
    statusTemp = true
    vioList = []
    for eigInd in 1:length(eigValList)
        eig = eigValList[eigInd]
        if eig.re > 0
            statusTemp = false
            push!(vioList,eigVectorList[eigInd,:])
        end
    end

    @test statusTemp
    @test isempty(vioList)

    if !statusTemp
        Amg = obtainGlobal_var(mpData,pm,rN,omega0)
        constraint_stability(pm, 0, vioList, Amg)
    end
end
