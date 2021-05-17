# require to load inverter.jl first
"Given a matrix A with NL expression and an eigenvector v, calculate the multiplication v^\\top A v"
function matrixMulti(mp, A, v)
    nonzeroList = []
    for vi in 1:length(v)
        if v[vi] != 0
            push!(nonzeroList,vi)
        end
    end
    n = length(nonzeroList)
    multiExpression = JuMP.@NLexpression(mp, sum(sum(v[nonzeroList[i]] * A[nonzeroList[i],nonzeroList[j]] * v[nonzeroList[j]] for j in 1:n
            if (A[nonzeroList[i],nonzeroList[j]] != 0)&(v[nonzeroList[i]] != 0)&(v[nonzeroList[j]] != 0)) for i in 1:n))
    return multiExpression
end

"Generate the nonlinear cuts"
function add_stability_cut!(pm, A, v)
    vReal = real(v)
    nonzeroListR = []
    for vi in 1:length(vReal)
        if vReal[vi] != 0
            push!(nonzeroListR,vi)
        end
    end
    nReal = length(nonzeroListR)

    vImag = imag(v)
    nonzeroListI = []
    for vi in 1:length(vImag)
        if vImag[vi] != 0
            push!(nonzeroListI,vi)
        end
    end
    nImag = length(nonzeroListI)

    # JuMP.@NLconstraint(pm.model, sum(sum(vReal[nonzeroListR[i]] * A[nonzeroListR[i],nonzeroListR[j]] * vReal[nonzeroListR[j]] for j in 1:nReal
    #         if (A[nonzeroListR[i],nonzeroListR[j]] != 0)&(vReal[nonzeroListR[i]] != 0)&(vReal[nonzeroListR[j]] != 0)) for i in 1:nReal) +
    #         sum(sum(vImag[nonzeroListI[i]] * A[nonzeroListI[i],nonzeroListI[j]] * vImag[nonzeroListI[j]] for j in 1:nImag
    #         if (A[nonzeroListI[i],nonzeroListI[j]] != 0)&(vImag[nonzeroListI[i]] != 0)&(vImag[nonzeroListI[j]] != 0)) for i in 1:nImag)
    #         <= 0)
    # lhsVar = JuMP.@NLexpression(pm.model, 0)
    lhsConst = 0
    varListReal = []
    varListImag = []
    for i in 1:nReal
        for j in 1:nReal
            if (vReal[nonzeroListR[i]] != 0)&(vReal[nonzeroListR[j]] != 0)&(A[nonzeroListR[i],nonzeroListR[j]] != 0)
                if typeof(A[nonzeroListR[i],nonzeroListR[j]]) == JuMP.NonlinearExpression
                    if abs(vReal[nonzeroListR[i]]*vReal[nonzeroListR[j]]) > 1e-6
                        # lhsVar += vReal[nonzeroListR[i]] * A[nonzeroListR[i],nonzeroListR[j]] * vReal[nonzeroListR[j]]
                        push!(varListReal,(i,j))
                    end
                else
                    if abs(vReal[nonzeroListR[i]] * A[nonzeroListR[i],nonzeroListR[j]] * vReal[nonzeroListR[j]]) > 1e-6
                        lhsConst += vReal[nonzeroListR[i]] * A[nonzeroListR[i],nonzeroListR[j]] * vReal[nonzeroListR[j]]
                    end
                end
            end
        end
    end
    for i in 1:nImag
        for j in 1:nImag
            if (vImag[nonzeroListI[i]] != 0)&(vImag[nonzeroListI[j]] != 0)&(A[nonzeroListI[i],nonzeroListI[j]] != 0)
                if typeof(A[nonzeroListR[i],nonzeroListR[j]]) == JuMP.NonlinearExpression
                    if abs(vImag[nonzeroListI[i]]*vImag[nonzeroListI[j]]) > 1e-6
                        # lhsVar += vImag[nonzeroListI[i]] * A[nonzeroListI[i],nonzeroListI[j]] * vImag[nonzeroListI[j]]
                        push!(varListReal,(i,j))
                    end
                else
                    if abs(vImag[nonzeroListI[i]] * A[nonzeroListI[i],nonzeroListI[j]] * vImag[nonzeroListI[j]]) > 1e-6
                        lhsConst += vImag[nonzeroListI[i]] * A[nonzeroListI[i],nonzeroListI[j]] * vImag[nonzeroListI[j]]
                    end
                end
            end
        end
    end
    JuMP.@NLconstraint(pm.model, sum(vReal[nonzeroListR[item[1]]] * A[nonzeroListR[item[1]],nonzeroListR[item[2]]] * vReal[nonzeroListR[item[2]]] for item in varListReal) +
        sum(vImag[nonzeroListI[item[1]]] * A[nonzeroListI[item[1]],nonzeroListI[item[2]]] * vImag[nonzeroListI[item[2]]] for item in varListImag) +
        lhsConst <= 0)
end

"Obtain the small-signal stability matrix used to generate eigen-value constraints"
function obtainGlobal_var(mpData, pm, rN, omega0, invData)
    busList, brList, invList, invConnected, invLine, loadList, vnomList, loadConnections, load_L, load_R, load_X, inverters, invBusDict = invData

    # obtain A matrix
    Asub = obtainA_inverter_global_var(mpData, pm, rN, omega0, busList, invList, invLine, invConnected, inverters, invBusDict)
    Bsub = obtainB_inverter_global(mpData, rN, omega0, busList, brList, invList, invLine, invConnected, inverters, invBusDict)
    Csub = obtainC_inverter_global(mpData, rN, omega0, busList, brList, invList, invLine, invConnected, load_L, loadConnections, inverters, invBusDict)
    Dsub = obtainD_inverter_global(mpData, rN, omega0, busList, brList, invList, invLine, inverters, invBusDict)
    Esub = obtainE_inverter_global(mpData, rN, omega0, brList)
    Fsub = obtainF_inverter_global(mpData, rN, omega0, busList, brList, loadList, loadConnections)
    Gsub = obtainG_inverter_global(mpData, rN, omega0, busList, brList, invList, loadList, load_L, loadConnections, invLine, inverters, invBusDict)
    Hsub = obtainH_inverter_global(mpData, rN, omega0, busList, brList, load_L, loadConnections)
    Isub = obtainI_inverter_global(mpData, rN, omega0, busList, brList, loadList, load_L, load_R, load_X, loadConnections)

    Atot = combineSub(busList, brList, inverters, invBusDict, Asub, Bsub, Csub, Dsub, Esub, Fsub, Gsub, Hsub, Isub, 2)
    return Atot
end

"Rounding very small values in the eigenvector to be 0"
function round0eig(ev)
    evNew = deepcopy(ev)
    for i in 1:length(ev)
        item = ev[i];
        if abs(real(item)) <= 1e-5
            realItem = 0
        else
            realItem = real(item)
        end
        if abs(imag(item)) <= 1e-5
            imagItem = 0
        else
            imagItem = imag(item)
        end
        evNew[i] = realItem + im*imagItem
    end
    return evNew;
end

"Create constraints to append to the power models"
function constraint_stability(pm::PowerModels.AbstractPowerModel, nw::Int, eigenVectorList::Array{Any,1}, Amg::Array{Any,2})
    # create constraint generate by the eigenvalue analysis
    for ev in eigenVectorList
        evRound = round0eig(ev)
        # evMulti_real = matrixMulti(pm.model, Amg, real(evRound))
        # evMulti_imag = matrixMulti(pm.model, Amg, imag(evRound))
        add_stability_cut!(pm, Amg, evRound)
        # JuMP.@NLconstraint(pm.model, evMulti_real + evMulti_imag <= 0)
    end
end
