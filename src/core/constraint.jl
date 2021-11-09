# require to load inverter.jl first
"Given a matrix A with NL expression and an eigenvector v, calculate the multiplication v^\\top A v"
function matrixMulti(mp, A, v)
    n = length(v)
    multiExpression = JuMP.@NLexpression(mp, sum(sum(v[i] * A[i,j] * v[j] for j in 1:n) for i in 1:n))

    return multiExpression
end

"Obtain the small-signal stability matrix used to generate eigen-value constraints"
function obtainGlobal_var(mpData, pm, rN, omega0)
    busList, brList, invList, invConnected, invLine, loadList, vnomList, loadConnections = PMS.preproc(mpData)
    load_L, load_R, load_X = PMS.procLoad(mpData, loadList, vnomList, omega0, loadConnections)

    inverters = []
    invBusDict = Dict()

    for i in invList
        currentInv = invConnected[i]
        # obtain a list of inverter buses
        append!(inverters, currentInv)
        for invItem in currentInv
            invBusDict[invItem] = i
        end
    end

    # obtain A matrix
    Asub = PMS.obtainA_inverter_global_var(mpData, pm, rN, omega0, busList, invList, invLine, invConnected, inverters, invBusDict)
    Bsub = PMS.obtainB_inverter_global(mpData, rN, omega0, busList, brList, invList, invLine, invConnected, inverters, invBusDict)
    Csub = PMS.obtainC_inverter_global(mpData, rN, omega0, busList, brList, invList, invLine, invConnected, load_L, loadConnections, inverters, invBusDict)
    Dsub = PMS.obtainD_inverter_global(mpData, rN, omega0, busList, brList, invList, invLine, inverters, invBusDict)
    Esub = PMS.obtainE_inverter_global(mpData, rN, omega0, brList)
    Fsub = PMS.obtainF_inverter_global(mpData, rN, omega0, busList, brList, loadList, loadConnections)
    Gsub = PMS.obtainG_inverter_global(mpData, rN, omega0, busList, brList, invList, loadList, load_L, loadConnections, invLine, inverters, invBusDict)
    Hsub = PMS.obtainH_inverter_global(mpData, rN, omega0, busList, brList, load_L, loadConnections)
    Isub = PMS.obtainI_inverter_global(mpData, rN, omega0, busList, brList, loadList, load_L, load_R, load_X, loadConnections)

    Atot = PMS.combineSub(busList, brList, inverters, invBusDict, Asub, Bsub, Csub, Dsub, Esub, Fsub, Gsub, Hsub, Isub, 2)
    
    return Atot
end

"Create constraints to append to the power models"
function constraint_stability(pm::_PMD.AbstractUnbalancedPowerModel, nw::Int, eigenVectorList::Array{Any,1}, Amg::Array{Any,2})
    # create constraint generate by the eigenvalue analysis
    for ev in eigenVectorList
        evMulti_real = PMS.matrixMulti(pm.model, Amg, real(ev))
        evMulti_imag = PMS.matrixMulti(pm.model, Amg, imag(ev))
        JuMP.@NLconstraint(pm.model, evMulti_real + evMulti_imag <= 0)
    end
end
