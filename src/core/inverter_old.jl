# inverter ODE construction
function vectorMulti(mp, C, v)
    # given a vector with NL expression and a constant matrix C, calculate the multiplication Cv
    n = length(v)
    multiExpression = [JuMP.JuMP.@NLexpression(mp, sum(C[i,j] * v[j] for j in 1:n)) for i in 1:n]
    return multiExpression
end

function obtainA_inverter_global(mpData, opfSol, rN, omega0, busList, invList, invLine, invConnected)
    A = Dict()

    for i in busList
        # initialize the matrix to be zero matrices, if it is a regular bus
        A[i] = zeros(9, 9)
    end

    for busConnected in invList
        # there should be only one branch from the generator
        i = invConnected[busConnected][1]
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

        A[busConnected][1,2] = 1
        A[busConnected][2,:] = -1 / mpData["bus"][i]["tau"] * [mpData["bus"][i]["mp"] * dP_delta
                            1
                            mpData["bus"][i]["mp"] * dP_v
                            mpData["bus"][i]["mp"] * dP_id[1]; mpData["bus"][i]["mp"] * dP_id[2]; mpData["bus"][i]["mp"] * dP_id[3]
                            mpData["bus"][i]["mp"] * dP_iq[1]; mpData["bus"][i]["mp"] * dP_iq[2]; mpData["bus"][i]["mp"] * dP_iq[3]
                            ]
        A[busConnected][3,:] = -1 / mpData["bus"][i]["tau"] * [mpData["bus"][i]["mq"] * dQ_delta
                            0
                            1 + mpData["bus"][i]["mq"] * dQ_v
                            mpData["bus"][i]["mq"] * dQ_id[1]; mpData["bus"][i]["mq"] * dQ_id[2]; mpData["bus"][i]["mp"] * dQ_id[3]
                            mpData["bus"][i]["mq"] * dQ_iq[1]; mpData["bus"][i]["mq"] * dQ_iq[2]; mpData["bus"][i]["mp"] * dQ_iq[3]
                            ]
        A[busConnected][4:6,1] = inv(L) * dvd_delta
        A[busConnected][4:6,3] = inv(L) * dvd_v
        A[busConnected][4:6,4:6] = -inv(L) * link["br_r"] - rN * inv(L)
        A[busConnected][4:6,7:9] = inv(L) * link["br_x"]
        A[busConnected][7:9,1] = inv(L) * dvq_delta
        A[busConnected][7:9,3] = inv(L) * dvq_v
        A[busConnected][7:9,4:6] = -inv(L) * link["br_x"]
        A[busConnected][7:9,7:9] = -inv(L) * link["br_r"] - rN * inv(L)
    end
    return A
end

function obtainA_inverter_global_var(mpData, pm, omega0, mP, mQ, tau, rN, busList, invList, invLine, invConnected)
    A = Dict()

    for i in busList
        # initialize the matrix to be zero matrices, if it is a regular bus
        A[i] = Array{Any,2}(zeros(9, 9))
    end

    for busConnected in invList
        # there should be only one branch from the generator
        iBus = invConnected[busConnected][0]
        iInd = parse(Int64, iBus)
        link = mpData["branch"][iBus]
        L = link["br_x"] ./ omega0
        i = mpData["bus"][iBus]["index"]
        delta0 = _PM.var(pm, 0, :va, i)[1]
        vmList = _PM.var(pm, 0, :vm, i)
        vaList = _PM.var(pm, 0, :va, i)
        vd0 = [JuMP.@NLexpression(pm.model,vmList.data[j] * cos(vaList.data[j])) for j in 1:length(vmList.data)]
        vq0 = [JuMP.@NLexpression(pm.model,vmList.data[j] * sin(vaList.data[j])) for j in 1:length(vmList.data)]
        P0 = JuMP.@NLexpression(pm.model,0)
        Q0 = JuMP.@NLexpression(pm.model,0)
        for g in keys(mpData["gen"])
            if "$(mpData["gen"][g]["gen_bus"])" == i
                P0 += sum(_PM.var(pm, 0, :pg, g))
                Q0 += sum(_PM.var(pm, 0, :qg, g))
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

        A[busConnected][1,2] = 1
        A[busConnected][2,:] = [JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * mpData["bus"][iBus]["mp"] * dP_delta)
                            -1 / mpData["bus"][iBus]["tau"]
                            JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * mpData["bus"][iBus]["mp"] * dP_v)
                            JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * mpData["bus"][iBus]["mp"] * dP_id[1])
                            JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * mpData["bus"][iBus]["mp"] * dP_id[2])
                            JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * mpData["bus"][iBus]["mp"] * dP_id[3])
                            JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * mpData["bus"][iBus]["mp"] * dP_iq[1])
                            JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * mpData["bus"][iBus]["mp"] * dP_iq[2])
                            JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * mpData["bus"][iBus]["mp"] * dP_iq[3])
                            ]
        A[busConnected][3,:] = [JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * mpData["bus"][iBus]["mq"] * dQ_delta)
                            0
                            JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * (1 + mpData["bus"][iBus]["mq"] * dQ_v))
                            JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * (mpData["bus"][iBus]["mq"] * dQ_id[1]))
                            JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * (mpData["bus"][iBus]["mq"] * dQ_id[2]))
                            JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * (mpData["bus"][iBus]["mp"] * dQ_id[3]))
                            JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * (mpData["bus"][iBus]["mq"] * dQ_iq[1]))
                            JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * (mpData["bus"][iBus]["mq"] * dQ_iq[2]))
                            JuMP.@NLexpression(pm.model,-1 / mpData["bus"][iBus]["tau"] * (mpData["bus"][iBus]["mp"] * dQ_iq[3]))
                            ]
        A[busConnected][4:6,1] = vectorMulti(pm.model, inv(L), dvd_delta)
        A[busConnected][4:6,3] = vectorMulti(pm.model, inv(L), dvd_v)
        A[busConnected][4:6,4:6] = -inv(L) * link["br_r"] - rN * inv(L)
        A[busConnected][4:6,7:9] = inv(L) * link["br_x"]
        A[busConnected][7:9,1] = vectorMulti(pm.model, inv(L), dvq_delta)
        A[busConnected][7:9,3] = vectorMulti(pm.model, inv(L), dvq_v)
        A[busConnected][7:9,4:6] = -inv(L) * link["br_x"]
        A[busConnected][7:9,7:9] = -inv(L) * link["br_r"] - rN * inv(L)
    end
    return A
end

function obtainB_inverter_global(mpData, rN, omega0, busList, brList, invList, invLine, invConnected)
    B = Dict()

    for i in busList
        for k in brList
            # initialize the matrix to be zero matrices
            B[i,k] = zeros(9, 6)
            iInd = parse(Int64, i)
            if i in invList
                L = mpData["branch"][invLine[invConnected[i][1]]]["br_x"] ./ omega0
                if mpData["branch"][k]["f_bus"] == iInd
                    B[i,k][4:6,1:3] = rN * inv(L)
                    B[i,k][7:9,4:6] = rN * inv(L)
                elseif mpData["branch"][k]["t_bus"] == iInd
                    B[i,k][4:6,1:3] = -rN * inv(L)
                    B[i,k][7:9,4:6] = -rN * inv(L)
                end
            end
        end
    end
    return B
end

function obtainC_inverter_global(mpData, rN, omega0, busList, brList, invList, invLine, load_L)
    C = Dict()

    for i in busList
        # initialize the matrix to be zero matrices
        C[i] = zeros(9, 6)
        iInd = parse(Int64, i)
        if (i in invList) & (i in keys(load_L))
            L = load_L[i]
            C[i][4:6,1:3] = rN * inv(L)
            C[i][7:9,4:6] = rN * inv(L)
        end
    end
    return C
end

function obtainD_inverter_global(mpData, rN, omega0, busList, brList, invList, invLine)
    D = Dict()

    for k in brList
        L = mpData["branch"][k]["br_x"] ./ omega0
        for i in busList
            # it is a regular bus
            D[k,i] = zeros(6, 9)
            if i == "$(mpData["branch"][k]["f_bus"])"
                D[k,i][1:3,4:6] = rN * inv(L)
                D[k,i][4:6,7:9] = rN * inv(L)
            end
            if i == "$(mpData["branch"][k]["t_bus"])"
                D[k,i][1:3,4:6] = -rN * inv(L)
                D[k,i][4:6,7:9] = -rN * inv(L)
            end
        end
    end
    return D
end

function obtainE_inverter_global(mpData, rN, omega0, brList)
    E = Dict()

    for k1 in brList
        L = mpData["branch"][k1]["br_x"] ./ omega0
        for k2 in brList
            E[k1,k2] = zeros(6, 6)
            if k1 == k2
                E[k1,k2][1:3,1:3] = -2 * rN * inv(L) - inv(L) * mpData["branch"][k1]["br_r"]
                E[k1,k2][1:3,4:6] = inv(L) * mpData["branch"][k1]["br_x"]
                E[k1,k2][4:6,1:3] = -inv(L) * mpData["branch"][k1]["br_x"]
                E[k1,k2][4:6,4:6] = -2 * rN * inv(L) - inv(L) * mpData["branch"][k1]["br_r"]
            elseif (mpData["branch"][k1]["f_bus"] == mpData["branch"][k2]["f_bus"]) | (mpData["branch"][k1]["t_bus"] == mpData["branch"][k2]["t_bus"])
                E[k1,k2][1:3,1:3] = -rN * inv(L)
                E[k1,k2][4:6,4:6] = -rN * inv(L)
            elseif mpData["branch"][k1]["f_bus"] == mpData["branch"][k2]["t_bus"]
                if mpData["branch"][k1]["t_bus"] == mpData["branch"][k2]["f_bus"]
                    E[k1,k2][1:3,1:3] = 2 * rN * inv(L)
                    E[k1,k2][4:6,4:6] = 2 * rN * inv(L)
                else
                    E[k1,k2][1:3,1:3] = rN * inv(L)
                    E[k1,k2][4:6,4:6] = rN * inv(L)
                end
            end
        end
    end

    return E
end

function obtainF_inverter_global(mpData, rN, omega0, busList, brList, loadList)
    F = Dict()

    for k in brList
        L = mpData["branch"][k]["br_x"] ./ omega0
        for i in busList
            F[k,i] = zeros(6, 6)
            # obtain line L
            if i in keys(loadList)
                if mpData["branch"][k]["f_bus"] == parse(Int64, i)
                    F[k,i][1:3,1:3] = -rN * inv(L)
                    F[k,i][4:6,4:6] = -rN * inv(L)
                elseif mpData["branch"][k]["t_bus"] == parse(Int64, i)
                    F[k,i][1:3,1:3] = rN * inv(L)
                    F[k,i][4:6,4:6] = rN * inv(L)
                end
            end
        end
    end

    return F
end

function obtainG_inverter_global(mpData, rN, omega0, busList, brList, invList, loadList, load_L)
    G = Dict()

    for i in busList
        G[i] = zeros(6, 9)
        if (i in invList) & (i in keys(loadList))
            L = load_L[i]
            G[i][1:3,4:6] = rN * inv(L)
            G[i][4:6,7:9] = rN * inv(L)
        end
    end

    return G
end

function obtainH_inverter_global(mpData, rN, omega0, busList, brList, load_L)
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
                        H[i,k][1:3,1:3] = -rN * inv(L)
                        H[i,k][4:6,4:6] = -rN * inv(L)
                    elseif mpData["branch"][k]["t_bus"] == parse(Int64, i)
                        H[i,k][1:3,1:3] = rN * inv(L)
                        H[i,k][4:6,4:6] = rN * inv(L)
                    end
                catch
                    println("L[$(i)] is not invertable!")
                end
            end
        end
    end

    return H
end

function obtainI_inverter_global(mpData, rN, omega0, busList, brList, loadList, load_L, load_R, load_X)
    I = Dict()

    for i in busList
        I[i] = zeros(6, 6)
        if i in keys(loadList)
            L = load_L[i]
            try
                I[i][1:3,1:3] = -rN * inv(L) - inv(L) * load_R[i]
                I[i][1:3,4:6] = inv(L) * load_X[i]
                I[i][4:6,1:3] = -inv(L) * load_X[i]
                I[i][1:3,4:6] = -rN * inv(L) - inv(L) * load_R[i]
            catch
                println("L[$(i)] is not invertable!")
            end
        end
    end
    return I
end

function combineSub(busList, brList, A, B, C, D, E, F, G, H, I, type)
    n = length(busList)
    m = length(brList)
    if type == 1
        Atot = zeros(15 * n + 6 * m, 15 * n + 6 * m)
    else
        Atot = Array{Any,2}(zeros(15 * n + 6 * m, 15 * n + 6 * m))
    end
    # Fill in the total matrix
    for i in 1:n
        # fill in A
        Atot[9 * (i - 1) + 1:9 * i,9 * (i - 1) + 1:9 * i] = A[busList[i]]
        # fill in C
        Atot[9 * (i - 1) + 1:9 * i,9 * n + 6 * m + 6 * (i - 1) + 1:9 * n + 6 * m + 6 * i] = C[busList[i]]
        # fill in G
        Atot[9 * n + 6 * m + 6 * (i - 1) + 1:9 * n + 6 * m + 6 * i,9 * (i - 1) + 1:9 * i] = G[busList[i]]
        # fill in I
        Atot[9 * n + 6 * m + 6 * (i - 1) + 1:9 * n + 6 * m + 6 * i,9 * n + 6 * m + 6 * (i - 1) + 1:9 * n + 6 * m + 6 * i] = I[busList[i]]
        for j in 1:m
            # fill in B
            Atot[9 * (i - 1) + 1:9 * i,9 * n + 6 * (j - 1) + 1:9 * n + 6 * j] = B[busList[i],brList[j]]
            # fill in D
            Atot[9 * n + 6 * (j - 1) + 1:9 * n + 6 * j,9 * (i - 1) + 1:9 * i] = D[brList[j],busList[i]]
            # fill in F
            Atot[9 * n + 6 * (j - 1) + 1:9 * n + 6 * j,9 * n + 6 * m + 6 * (i - 1) + 1:9 * n + 6 * m + 6 * i] = F[brList[j],busList[i]]
            # fill in H
            Atot[9 * n + 6 * m + 6 * (i - 1) + 1:9 * n + 6 * m + 6 * i,9 * n + 6 * (j - 1) + 1:9 * n + 6 * j] = H[busList[i],brList[j]]
        end
    end
    for j1 in 1:m
        for j2 in 1:m
            # fill in E
            Atot[9 * n + 6 * (j1 - 1) + 1:9 * n + 6 * j1,9 * n + 6 * (j2 - 1) + 1:9 * n + 6 * j2] = E[brList[j1],brList[j2]]
        end
    end
    return Atot
end

"obtain the global matrix of the small-signal stability control"
function obtainGlobal(mpData, opfSol, omega0, rN)
    # preprocessing
    busList, brList, invList, invConnected, invLine, loadList, vnomList = preproc(mpData)
    load_L, load_R, load_X = procLoad(mpData, loadList, vnomList, omega0)

    # obtain A matrix
    Asub = obtainA_inverter_global(mpData, opfSol, rN, omega0, busList, invList, invLine, invConnected)
    Bsub = obtainB_inverter_global(mpData, rN, omega0, busList, brList, invList, invLine, invConnected)
    Csub = obtainC_inverter_global(mpData, rN, omega0, busList, brList, invList, invLine, load_L)
    Dsub = obtainD_inverter_global(mpData, rN, omega0, busList, brList, invList, invLine)
    Esub = obtainE_inverter_global(mpData, rN, omega0, brList)
    Fsub = obtainF_inverter_global(mpData, rN, omega0, busList, brList, loadList)
    Gsub = obtainG_inverter_global(mpData, rN, omega0, busList, brList, invList, loadList, load_L)
    Hsub = obtainH_inverter_global(mpData, rN, omega0, busList, brList, load_L)
    Isub = obtainI_inverter_global(mpData, rN, omega0, busList, brList, loadList, load_L, load_R, load_X)

    Atot = combineSub(busList, brList, Asub, Bsub, Csub, Dsub, Esub, Fsub, Gsub, Hsub, Isub, 1)

    return Atot
end
