@info "running two-bus system tests"

@testset "test two-bus system" begin
    mpData = PMD.parse_file("../test/data/case2_diag.dss"; data_model=PMD.MATHEMATICAL)
    opfSol = PMD.run_mc_opf(mpData, PMD.ACPPowerModel, ipopt_solver)

    omega0 = 2 * pi * 60
    mP = Dict()
    mQ = Dict()
    mP["3"] = 0.3
    mQ["3"] = 0.3
    tau = Dict()
    tau["3"] = 1e-4
    rN = 1000

    @testset "A" begin
        inverterList = ["3"]

        # obtain the A matrix for every single inverter
        A, errorList = obtainA_inverter(mpData, opfSol, inverterList, omega0, mP, mQ, tau)
        statusList = Dict()
        for i in inverterList
            # perform the eigenvalue post-analysis
            eigValList = eigvals(A[i])
            statusTemp = true
            for eig in eigValList
                if eig.re > 0
                    statusTemp = false
                end
            end
            statusList[i] = statusTemp
        end
    end

    @testset "Atot" begin
        Atot = obtainGlobal(mpData, opfSol, omega0, mP, mQ, tau, rN)
        eigValList = eigvals(Atot)
        statusTemp = true
        for eig in eigValList
            if eig.re > 0
                statusTemp = false
            end
        end
    end
end
