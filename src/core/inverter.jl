import LinearAlgebra: eigvals, eigvecs, dot, inv, norm

# inverter ODE construction, this script allows multiple inverter buses connect to the same bus

"given a vector with NL expression and a constant matrix C, calculate the multiplication Cv"
function vectorMulti(mp, C, v)
    
    n = length(v)
    multiExpression = [JuMP.@NLexpression(mp, sum(C[i,j] * v[j] for j in 1:n)) for i in 1:n]

    return multiExpression
end

"invert the matrix"
function inverseMat(originMat, connection)
    
    invM = zeros(3, 3)

    if norm(originMat) > 0
        try
            invOri = inv(originMat)
            invM[connection,connection] = invOri
        catch
            for i in 1:length(connection)
                for j in 1:length(connection)
                    if abs(1 / originMat[i,j]) < Inf
                        invM[connection[i],connection[j]] = 1 / originMat[i,j]
                    end
                end
            end
        end
    end
    
    return invM
end

"Obtain the submatrix A, with the opf solution given"
function obtainA_inverter_global(mpData, opfSol, rN, omega0, busList, invList, invLine, invConnected, inverters, invBusDict)
    
    A = Dict()

    for invItem in inverters
        # initialize the matrix to be zero matrices, if it is a regular bus
        A[invItem,invItem] = zeros(9, 9)
    end

    # for every inverter
    for i in inverters
        # there should be only one branch from the generator
        iInd = parse(Int64, i)
        link = mpData["branch"][invLine[i]]
        L = link["br_x"] ./ omega0
        delta0 = opfSol["solution"]["bus"][i]["va"][1]
        vd0 = opfSol["solution"]["bus"][i]["vm"] .* cos.(opfSol["solution"]["bus"][i]["va"])
        vq0 = opfSol["solution"]["bus"][i]["vm"] .* sin.(opfSol["solution"]["bus"][i]["va"])
        P0 = 0
        Q0 = 0
        for g in keys(mpData["gen"])
            if "$(mpData["gen"][g]["gen_bus"])" == i
                P0 += sum(opfSol["solution"]["gen"][g]["pg"])
                Q0 += sum(opfSol["solution"]["gen"][g]["qg"])
                # Question: is this consistent with the calculation from line flow?
            end
        end
        id0 = zeros(3)
        iq0 = zeros(3)
        for j in 1:3
            id0[j], iq0[j] = inv([vd0[j] vq0[j];vq0[j] -vd0[j]]) * [P0 / 3;Q0 / 3]
        end
        dvd_delta = -vq0
        dvq_delta = vd0
        dvd_v = cos.(opfSol["solution"]["bus"][i]["va"])
        dvq_v = sin.(opfSol["solution"]["bus"][i]["va"])
        dP_delta = dot(dvd_delta, id0) + dot(dvq_delta, iq0)
        dQ_delta = -dot(dvd_delta, iq0) + dot(dvq_delta, id0)
        dP_v = dot(dvd_v, id0) + dot(dvq_v, iq0)
        dQ_v = -dot(dvd_v, iq0) + dot(dvq_v, id0)
        dP_id = vd0
        dP_iq = vq0
        dQ_id = vq0
        dQ_iq = -vd0

        A[i,i][1,2] = 1
        A[i,i][2,:] = -1 / mpData["bus"][i]["tau"] * [mpData["bus"][i]["mp"] * dP_delta
                            1
                            mpData["bus"][i]["mp"] * dP_v
                            mpData["bus"][i]["mp"] * dP_id[1]; mpData["bus"][i]["mp"] * dP_id[2]; mpData["bus"][i]["mp"] * dP_id[3]
                            mpData["bus"][i]["mp"] * dP_iq[1]; mpData["bus"][i]["mp"] * dP_iq[2]; mpData["bus"][i]["mp"] * dP_iq[3]
                            ]
        A[i,i][3,:] = -1 / mpData["bus"][i]["tau"] * [mpData["bus"][i]["mq"] * dQ_delta
                            0
                            1 + mpData["bus"][i]["mq"] * dQ_v
                            mpData["bus"][i]["mq"] * dQ_id[1]; mpData["bus"][i]["mq"] * dQ_id[2]; mpData["bus"][i]["mp"] * dQ_id[3]
                            mpData["bus"][i]["mq"] * dQ_iq[1]; mpData["bus"][i]["mq"] * dQ_iq[2]; mpData["bus"][i]["mp"] * dQ_iq[3]
                            ]
        A[i,i][4:6,1] = inv(L) * dvd_delta
        A[i,i][4:6,3] = inv(L) * dvd_v
        A[i,i][4:6,4:6] = -inv(L) * link["br_r"] - rN * inv(L)
        A[i,i][4:6,7:9] = inv(L) * link["br_x"]
        A[i,i][7:9,1] = inv(L) * dvq_delta
        A[i,i][7:9,3] = inv(L) * dvq_v
        A[i,i][7:9,4:6] = -inv(L) * link["br_x"]
        A[i,i][7:9,7:9] = -inv(L) * link["br_r"] - rN * inv(L)
    end

    # for every pair of inverters in the same bus
    for iBus in invList
        if length(invConnected[iBus]) > 1
            for i in 1:length(invConnected[iBus])
                for j in (i + 1):length(invConnected[iBus])
                    A[i,j] = zeros(9, 9)
                    link = mpData["branch"][invLine[i]]
                    L = link["br_x"] ./ omega0
                    A[i,j][4:6,4:6] = -inv(L) * link["br_r"] - rN * inv(L)
                    A[i,j][7:9,7:9] = -inv(L) * link["br_r"] - rN * inv(L)

                    link = mpData["branch"][invLine[j]]
                    L = link["br_x"] ./ omega0
                    A[j,i] = zeros(9, 9)
                    A[j,i][4:6,4:6] = -inv(L) * link["br_r"] - rN * inv(L)
                    A[j,i][7:9,7:9] = -inv(L) * link["br_r"] - rN * inv(L)
                end
            end
        end
    end

    return A
end

"Obtain the submatrix A, with the opf problem variables in"
function obtainA_inverter_global_var(mpData, pm, rN, omega0, busList, invList, invLine, invConnected, inverters, invBusDict)
    
    A = Dict()
    for invItem in inverters
        A[invItem,invItem] = Array{Any,2}(zeros(9, 9))
    end

    for iBus in inverters
        # there should be only one branch from the generator
        iInd = parse(Int64, iBus)
        link = mpData["branch"][invLine[iBus]]
        L = link["br_x"] ./ omega0
        i = mpData["bus"][iBus]["index"]
        delta0 = _PMD.var(pm, 0, :va, i)[1]
        vmList = _PMD.var(pm, 0, :vm, i)
        vaList = _PMD.var(pm, 0, :va, i)
        vd0 = [JuMP.@NLexpression(pm.model,vmList.data[j] * cos(vaList.data[j])) for j in 1:length(vmList.data)]
        vq0 = [JuMP.@NLexpression(pm.model,vmList.data[j] * sin(vaList.data[j])) for j in 1:length(vmList.data)]
        P0 = JuMP.@NLexpression(pm.model,0)
        Q0 = JuMP.@NLexpression(pm.model,0)
        for g in keys(mpData["gen"])
            if "$(mpData["gen"][g]["gen_bus"])" == i
                P0 += sum(_PMD.var(pm, 0, :pg, g))
                Q0 += sum(_PMD.var(pm, 0, :qg, g))
            end
        end
        id0 = Array{Any,1}([0,0,0])
        iq0 = Array{Any,1}([0,0,0])
        for j in 1:3
            id0[j] = JuMP.@NLexpression(pm.model,(vd0[j] * P0 / 3 + vq0[j] * Q0 / 3) / (vmList.data[j]^2))
            iq0[j] = JuMP.@NLexpression(pm.model,(vq0[j] * P0 / 3 - vd0[j] * Q0 / 3) / (vmList.data[j]^2))
        end
        dvd_delta = [JuMP.@NLexpression(pm.model,-vq0[j]) for j in 1:length(vq0)]
        dvq_delta = [JuMP.@NLexpression(pm.model,vd0[j]) for j in 1:length(vd0)]
        dvd_v = [JuMP.@NLexpression(pm.model,cos(vaList.data[j])) for j in 1:length(vaList.data)]
        dvq_v = [JuMP.@NLexpression(pm.model,sin(vaList.data[j])) for j in 1:length(vaList.data)]
        dP_delta = JuMP.@NLexpression(pm.model,sum(dvd_delta[j] * id0[j] for j in 1:length(dvd_delta)) + sum(dvq_delta[j] * iq0[j] for j in 1:length(dvq_delta)))
        dQ_delta = JuMP.@NLexpression(pm.model,-sum(dvd_delta[j] * iq0[j] for j in 1:length(dvd_delta)) + sum(dvq_delta[j] * id0[j] for j in 1:length(dvq_delta)))
        dP_v = JuMP.@NLexpression(pm.model,sum(dvd_v[j] * id0[j] for j in 1:length(dvd_v)) + sum(dvq_v[j] * iq0[j] for j in 1:length(dvq_v)))
        dQ_v = JuMP.@NLexpression(pm.model,-sum(dvd_v[j] * iq0[j] for j in 1:length(dvd_v)) + sum(dvq_v[j] * id0[j] for j in 1:length(dvq_v)))
        dP_id = vd0
        dP_iq = vq0
        dQ_id = vq0
        dQ_iq = [JuMP.@NLexpression(pm.model,-vd0[i]) for i in 1:length(vd0)]

        A[iBus,iBus][1,2] = 1
        A[iBus,iBus][2,:] = [JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * mpData["bus"][iBus]["mp"] * dP_delta)
                            -1 / mpData["bus"][iBus]["tau"]
                            JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * mpData["bus"][iBus]["mp"] * dP_v)
                            JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * mpData["bus"][iBus]["mp"] * dP_id[1])
                            JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * mpData["bus"][iBus]["mp"] * dP_id[2])
                            JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * mpData["bus"][iBus]["mp"] * dP_id[3])
                            JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * mpData["bus"][iBus]["mp"] * dP_iq[1])
                            JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * mpData["bus"][iBus]["mp"] * dP_iq[2])
                            JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * mpData["bus"][iBus]["mp"] * dP_iq[3])
                            ]
        A[iBus,iBus][3,:] = [JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * mpData["bus"][iBus]["mq"] * dQ_delta)
                            0
                            JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * (1 + mpData["bus"][iBus]["mq"] * dQ_v))
                            JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * (mpData["bus"][iBus]["mq"] * dQ_id[1]))
                            JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * (mpData["bus"][iBus]["mq"] * dQ_id[2]))
                            JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * (mpData["bus"][iBus]["mp"] * dQ_id[3]))
                            JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * (mpData["bus"][iBus]["mq"] * dQ_iq[1]))
                            JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * (mpData["bus"][iBus]["mq"] * dQ_iq[2]))
                            JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * (mpData["bus"][iBus]["mp"] * dQ_iq[3]))
                            ]
        A[iBus,iBus][4:6,1] = PMS.vectorMulti(pm.model, inv(L), dvd_delta)
        A[iBus,iBus][4:6,3] = PMS.vectorMulti(pm.model, inv(L), dvd_v)
        A[iBus,iBus][4:6,4:6] = -inv(L) * link["br_r"] - rN * inv(L)
        A[iBus,iBus][4:6,7:9] = inv(L) * link["br_x"]
        A[iBus,iBus][7:9,1] = PMS.vectorMulti(pm.model, inv(L), dvq_delta)
        A[iBus,iBus][7:9,3] = PMS.vectorMulti(pm.model, inv(L), dvq_v)
        A[iBus,iBus][7:9,4:6] = -inv(L) * link["br_x"]
        A[iBus,iBus][7:9,7:9] = -inv(L) * link["br_r"] - rN * inv(L)
    end

    # for every pair of inverters in the same bus
    for iBus in invList
        if length(invConnected[iBus]) > 1
            for i in 1:length(invConnected[iBus])
                for j in (i + 1):length(invConnected[iBus])
                    A[i,j] = zeros(9, 9)
                    link = mpData["branch"][invLine[i]]
                    L = link["br_x"] ./ omega0
                    A[i,j][4:6,4:6] = -inv(L) * link["br_r"] - rN * inv(L)
                    A[i,j][7:9,7:9] = -inv(L) * link["br_r"] - rN * inv(L)

                    link = mpData["branch"][invLine[j]]
                    L = link["br_x"] ./ omega0
                    A[j,i] = zeros(9, 9)
                    A[j,i][4:6,4:6] = -inv(L) * link["br_r"] - rN * inv(L)
                    A[j,i][7:9,7:9] = -inv(L) * link["br_r"] - rN * inv(L)
                end
            end
        end
            end

    return A
end

"Obtain the submatrix B, with the opf solution given"
function obtainB_inverter_global(mpData, rN, omega0, busList, brList, invList, invLine, invConnected, inverters, invBusDict)

    B = Dict()

    for i in inverters
        iBus = invBusDict[i]
        for k in brList
            # initialize the matrix to be zero matrices
            B[i,k] = zeros(9, 6)
            iInd = parse(Int64, iBus)
            L = mpData["branch"][invLine[i]]["br_x"] ./ omega0
            if mpData["branch"][k]["f_bus"] == iInd
                B[i,k][4:6,1:3] = rN * PMS.inverseMat(L, mpData["branch"][invLine[i]]["f_connections"])
                B[i,k][7:9,4:6] = rN * PMS.inverseMat(L, mpData["branch"][invLine[i]]["f_connections"])
            elseif mpData["branch"][k]["t_bus"] == iInd
                B[i,k][4:6,1:3] = -rN * PMS.inverseMat(L, mpData["branch"][invLine[i]]["t_connections"])
                B[i,k][7:9,4:6] = -rN * PMS.inverseMat(L, mpData["branch"][invLine[i]]["t_connections"])
            end
        end
    end

    return B
end

"Obtain the submatrix C, with the opf solution given"
function obtainC_inverter_global(mpData, rN, omega0, busList, brList, invList, invLine, invConnected, load_L, loadConnections, inverters, invBusDict)
    
    C = Dict()

    for i in inverters
        # initialize the matrix to be zero matrices
        C[i,invBusDict[i]] = zeros(9, 6)
        iInd = parse(Int64, i)
        if invBusDict[i] in keys(load_L)
            L = load_L[invBusDict[i]]
            C[i,invBusDict[i]][4:6,1:3] = rN * PMS.inverseMat(L, loadConnections[invBusDict[i]])
            C[i,invBusDict[i]][7:9,4:6] = rN * PMS.inverseMat(L, loadConnections[invBusDict[i]])
        end
    end

    return C
end

    "Obtain the submatrix D, with the opf solution given"
function obtainD_inverter_global(mpData, rN, omega0, busList, brList, invList, invLine, inverters, invBusDict)
    
    D = Dict()

    for k in brList
        L = mpData["branch"][k]["br_x"] ./ omega0
        for i in inverters
            iBus = invBusDict[i]
            # it is a regular bus
            D[k,i] = zeros(6, 9)
            if iBus == "$(mpData["branch"][k]["f_bus"])"
                D[k,i][1:3,4:6] = rN * PMS.inverseMat(L, mpData["branch"][k]["f_connections"])
                D[k,i][4:6,7:9] = rN * PMS.inverseMat(L, mpData["branch"][k]["f_connections"])
            end
            if iBus == "$(mpData["branch"][k]["t_bus"])"
                D[k,i][1:3,4:6] = -rN * PMS.inverseMat(L, mpData["branch"][k]["t_connections"])
                D[k,i][4:6,7:9] = -rN * PMS.inverseMat(L, mpData["branch"][k]["t_connections"])
            end
        end
    end

    return D
end

"Obtain the submatrix E, with the opf solution given"
function obtainE_inverter_global(mpData, rN, omega0, brList)
    E = Dict()

    for k1 in brList
        L = mpData["branch"][k1]["br_x"] ./ omega0
        for k2 in brList
            E[k1,k2] = zeros(6, 6)
            if k1 == k2
                E[k1,k2][1:3,1:3] = -2 * rN * PMS.inverseMat(L, mpData["branch"][k1]["f_connections"]) -
                    PMS.inverseMat(L * mpData["branch"][k1]["br_r"], mpData["branch"][k1]["f_connections"])
                E[k1,k2][1:3,4:6] = PMS.inverseMat(L * mpData["branch"][k1]["br_x"], mpData["branch"][k1]["f_connections"])
                E[k1,k2][4:6,1:3] = -PMS.inverseMat(L * mpData["branch"][k1]["br_x"], mpData["branch"][k1]["f_connections"])
                E[k1,k2][4:6,4:6] = -2 * rN * PMS.inverseMat(L, mpData["branch"][k1]["f_connections"]) -
                    PMS.inverseMat(L * mpData["branch"][k1]["br_r"], mpData["branch"][k1]["f_connections"])

            elseif (mpData["branch"][k1]["f_bus"] == mpData["branch"][k2]["f_bus"])
                E[k1,k2][1:3,1:3] = -rN * PMS.inverseMat(L, mpData["branch"][k1]["f_connections"])
                E[k1,k2][4:6,4:6] = -rN * PMS.inverseMat(L, mpData["branch"][k1]["f_connections"])

            elseif (mpData["branch"][k1]["t_bus"] == mpData["branch"][k2]["t_bus"])
                E[k1,k2][1:3,1:3] = -rN * PMS.inverseMat(L, mpData["branch"][k1]["t_connections"])
                E[k1,k2][4:6,4:6] = -rN * PMS.inverseMat(L, mpData["branch"][k1]["t_connections"])

            elseif mpData["branch"][k1]["f_bus"] == mpData["branch"][k2]["t_bus"]
                if mpData["branch"][k1]["t_bus"] == mpData["branch"][k2]["f_bus"]
                    E[k1,k2][1:3,1:3] = 2 * rN * PMS.inverseMat(L, mpData["branch"][k1]["f_connections"])
                    E[k1,k2][4:6,4:6] = 2 * rN * PMS.inverseMat(L, mpData["branch"][k1]["f_connections"])
                else
                    E[k1,k2][1:3,1:3] = rN * PMS.inverseMat(L, mpData["branch"][k1]["f_connections"])
                    E[k1,k2][4:6,4:6] = rN * PMS.inverseMat(L, mpData["branch"][k1]["f_connections"])
                end
            end
        end
    end

    return E
end

"Obtain the submatrix F, with the opf solution given"
function obtainF_inverter_global(mpData, rN, omega0, busList, brList, loadList, loadConnections)
    
    F = Dict()

    for k in brList
        L = mpData["branch"][k]["br_x"] ./ omega0
        for i in busList
            F[k,i] = zeros(6, 6)
            # obtain line L
            if i in keys(loadList)
                if mpData["branch"][k]["f_bus"] == parse(Int64, i)
                    F[k,i][1:3,1:3] = -rN * PMS.inverseMat(L, mpData["branch"][k]["f_connections"])
                    F[k,i][4:6,4:6] = -rN * PMS.inverseMat(L, mpData["branch"][k]["f_connections"])
                elseif mpData["branch"][k]["t_bus"] == parse(Int64, i)
                    F[k,i][1:3,1:3] = rN * PMS.inverseMat(L, mpData["branch"][k]["t_connections"])
                    F[k,i][4:6,4:6] = rN * PMS.inverseMat(L, mpData["branch"][k]["t_connections"])
                end
            end
        end
    end

    return F
end

"Obtain the submatrix G, with the opf solution given"
function obtainG_inverter_global(mpData, rN, omega0, busList, brList, invList, loadList, load_L, loadConnections, invLine, inverters, invBusDict)
    G = Dict()

    for i in inverters
        iBus = invBusDict[i]
        G[iBus,i] = zeros(6, 9)
        if iBus in keys(loadList)
            L = load_L[iBus]
            G[iBus,i][1:3,4:6] = rN * PMS.inverseMat(L[loadConnections[iBus],loadConnections[iBus]], loadConnections[iBus])
            G[iBus,i][4:6,7:9] = rN * PMS.inverseMat(L[loadConnections[iBus],loadConnections[iBus]], loadConnections[iBus])
        end
    end

    return G
end

"Obtain the submatrix H, with the opf solution given"
function obtainH_inverter_global(mpData, rN, omega0, busList, brList, load_L, loadConnections)
    H = Dict()

    for i in busList
        for k in brList
            H[i,k] = zeros(6, 6)
            # obtain line L
            if i in keys(load_L)
                L = load_L[i]
                try
                    # if L != 0
                    if mpData["branch"][k]["f_bus"] == parse(Int64, i)
                        H[i,k][1:3,1:3] = -rN * PMS.inverseMat(L[loadConnections[i],loadConnections[i]], loadConnections[i])
                        H[i,k][4:6,4:6] = -rN * PMS.inverseMat(L[loadConnections[i],loadConnections[i]], loadConnections[i])
                    elseif mpData["branch"][k]["t_bus"] == parse(Int64, i)
                        H[i,k][1:3,1:3] = rN * PMS.inverseMat(L[loadConnections[i],loadConnections[i]], loadConnections[i])
                        H[i,k][4:6,4:6] = rN * PMS.inverseMat(L[loadConnections[i],loadConnections[i]], loadConnections[i])
                    end
                catch
                    Memento.warn(_LOGGER, "L[$(i)] is not invertable! for branch $(k)")
                end
            end
        end
    end

    return H
end

"Obtain the submatrix I, with the opf solution given"
function obtainI_inverter_global(mpData, rN, omega0, busList, brList, loadList, load_L, load_R, load_X, loadConnections)
    I = Dict()

    for i in busList
        I[i] = zeros(6, 6)
        if i in keys(loadList)
            L = load_L[i]
            try
                I[i][1:3,1:3] = -rN * PMS.inverseMat(L[loadConnections[i],loadConnections[i]], loadConnections[i]) -
                    PMS.inverseMat(L[loadConnections[i],loadConnections[i]], loadConnections[i]) * load_R[i]
                I[i][1:3,4:6] = PMS.inverseMat(L[loadConnections[i],loadConnections[i]], loadConnections[i]) * load_X[i]
                I[i][4:6,1:3] = -PMS.inverseMat(L[loadConnections[i],loadConnections[i]], loadConnections[i]) * load_X[i]
                I[i][1:3,4:6] = -rN * PMS.inverseMat(L[loadConnections[i],loadConnections[i]], loadConnections[i]) -
                    PMS.inverseMat(L[loadConnections[i],loadConnections[i]], loadConnections[i]) * load_R[i]
            catch
                Memento.warn(_LOGGER, "L[$(i)] is not invertable!")
            end
        end
    end

    return I
end

"Combine all submatrices to one stability control matrix"
function combineSub(busList, brList, inverters, invBusDict, A, B, C, D, E, F, G, H, I, type)
    
    n = length(busList)
    ni = length(inverters)
    m = length(brList)

    if type == 1
        Atot = zeros(9 * ni + 6 * m + 6 * n, 9 * ni + 6 * m + 6 * n)
        else
        Atot = Array{Any,2}(zeros(9 * ni + 6 * m + 6 * n, 9 * ni + 6 * m + 6 * n))
    end

    # Fill in the total matrix
    for i in 1:ni
        # fill in A
        for j in 1:ni
            if (inverters[i], inverters[j]) in keys(A)
                Atot[9 * (i - 1) + 1:9 * i,9 * (j - 1) + 1:9 * j] = A[inverters[i],inverters[j]]
            end
        end
        for j in 1:m
            # fill in B
            Atot[9 * (i - 1) + 1:9 * i,9 * ni + 6 * (j - 1) + 1:9 * ni + 6 * j] = B[inverters[i],brList[j]]
            # fill in D
            Atot[9 * ni + 6 * (j - 1) + 1:9 * ni + 6 * j,9 * (i - 1) + 1:9 * i] = D[brList[j],inverters[i]]
        end

        for j in 1:n
            if busList[j] == invBusDict[inverters[i]]
                # fill in C
                Atot[9 * (i - 1) + 1:9 * i,9 * ni + 6 * m + 6 * (j - 1) + 1:9 * ni + 6 * m + 6 * j] = C[inverters[i],busList[j]]
                # fill in G
                Atot[9 * ni + 6 * m + 6 * (j - 1) + 1:9 * ni + 6 * m + 6 * j,9 * (i - 1) + 1:9 * i] = G[busList[j],inverters[i]]
            end
        end
    end

    for i in 1:n
        # fill in I
        Atot[9 * ni + 6 * m + 6 * (i - 1) + 1:9 * ni + 6 * m + 6 * i,9 * ni + 6 * m + 6 * (i - 1) + 1:9 * ni + 6 * m + 6 * i] = I[busList[i]]
        for j in 1:m
            # fill in F
            Atot[9 * ni + 6 * (j - 1) + 1:9 * ni + 6 * j,9 * ni + 6 * m + 6 * (i - 1) + 1:9 * ni + 6 * m + 6 * i] = F[brList[j],busList[i]]
            # fill in H
            Atot[9 * ni + 6 * m + 6 * (i - 1) + 1:9 * ni + 6 * m + 6 * i,9 * ni + 6 * (j - 1) + 1:9 * ni + 6 * j] = H[busList[i],brList[j]]
        end
    end

    for j1 in 1:m
        for j2 in 1:m
        # fill in E
            Atot[9 * ni + 6 * (j1 - 1) + 1:9 * ni + 6 * j1,9 * ni + 6 * (j2 - 1) + 1:9 * ni + 6 * j2] = E[brList[j1],brList[j2]]
        end
    end

    return Atot
end

"obtain the global matrix of the small-signal stability control"
function obtainGlobal_multi(mpData, opfSol, rN, omega0)
    
    # preprocessing
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
    Asub = PMS.obtainA_inverter_global(mpData, opfSol, rN, omega0, busList, invList, invLine, invConnected, inverters, invBusDict)
    Bsub = PMS.obtainB_inverter_global(mpData, rN, omega0, busList, brList, invList, invLine, invConnected, inverters, invBusDict)
    Csub = PMS.obtainC_inverter_global(mpData, rN, omega0, busList, brList, invList, invLine, invConnected, load_L, loadConnections, inverters, invBusDict)
    Dsub = PMS.obtainD_inverter_global(mpData, rN, omega0, busList, brList, invList, invLine, inverters, invBusDict)
    Esub = PMS.obtainE_inverter_global(mpData, rN, omega0, brList)
    Fsub = PMS.obtainF_inverter_global(mpData, rN, omega0, busList, brList, loadList, loadConnections)
    Gsub = PMS.obtainG_inverter_global(mpData, rN, omega0, busList, brList, invList, loadList, load_L, loadConnections, invLine, inverters, invBusDict)
    Hsub = PMS.obtainH_inverter_global(mpData, rN, omega0, busList, brList, load_L, loadConnections)
    Isub = PMS.obtainI_inverter_global(mpData, rN, omega0, busList, brList, loadList, load_L, load_R, load_X, loadConnections)

    Atot = PMS.combineSub(busList, brList, inverters, invBusDict, Asub, Bsub, Csub, Dsub, Esub, Fsub, Gsub, Hsub, Isub, 1)

    return Atot
end
