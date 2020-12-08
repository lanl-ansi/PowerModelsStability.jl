"obtain the A matrix for diesel generator"
function obtainA_diesel(mpData,opfSol,iBusList,omega0,kD)
    # obtain the voltage solution
    T = 2/3*[1 cos(-2*pi/3) cos(2*pi/3);
             0 -sin(-2*pi/3) -sin(2*pi/3);
             1/2 1/2 1/2]
    vD = Dict()
    vQ = Dict()
    vd = Dict()
    vq = Dict()
    delta0 = Dict()
    v0 = Dict()
    id = Dict()
    iq = Dict()
    iD = Dict()
    iQ = Dict()
    A = Dict()
    errorList = []
    for i in iBusList
        if mpData["bus"][i]["bus_type"] == 3
            delta0[i] = opfSol["solution"]["bus"][i]["va"][1]
            vabcComplex = opfSol["solution"]["bus"][i]["vm"].*cos.(opfSol["solution"]["bus"][i]["va"]) +
                im*opfSol["solution"]["bus"][i]["vm"].*sin.(opfSol["solution"]["bus"][i]["va"])
            vDQComplex = T*vabcComplex
            if abs(vDQComplex[1].re) > 1e-8
                vD[i] = vDQComplex[1].re
            else
                vD[i] = 0
            end
            if abs(vDQComplex[1].im) > 1e-8
                vQ[i] = vDQComplex[1].im
            else
                vQ[i] = 0
            end
            vd[i] = vD[i]*cos(delta0[i]) + vQ[i]*sin(delta0[i])
            vq[i] = vQ[i]*cos(delta0[i]) - vD[i]*sin(delta0[i])
            # if it is a generator bus, search for the matching generator
            pg[i] = 0
            qg[i] = 0
            for g in keys(mpData["gen"])
                if "$(mpData["gen"][g]["gen_bus"])" == i
                    pg[i] += sum(opfSol["solution"]["gen"][g]["pg"])
                    qg[i] += sum(opfSol["solution"]["gen"][g]["qg"])
                    # Question: is this consistent with the calculation from line flow?
                end
            end
            id[i] = (pg[i]*vd[i] + qg[i]*vq[i])/(vd[i]^2 + vq[i]^2)
            iq[i] = (pg[i]*vq[i] - qg[i]*vd[i])/(vd[i]^2 + vq[i]^2)

            A[i] = zeros(2,2)
            A[i][1,2] = omega0
            A[i][2,2] = -kD[i]
            A[i][2,1] = -(2*vd[i]*vq[i])/(omega0*(vd[i]^2+vq[i]^2))*(pg[i] - qg[i])
        else
            push!(errorList,i)
        end
    end

    return A,errorList
end


""
function preproc(mpData)
    # obtain the branches that are between two buses and buses with inverters
    brList = []
    invList = []
    invLine = Dict()
    invConnected = Dict()
    for k in keys(mpData["branch"])
        t_bus = mpData["branch"][k]["t_bus"]
        f_bus = mpData["branch"][k]["f_bus"]
        if mpData["bus"]["$(t_bus)"]["bus_type"] != 3 && mpData["bus"]["$(f_bus)"]["bus_type"] != 3
            push!(brList,k)
        elseif mpData["bus"]["$(t_bus)"]["bus_type"] == 3
            push!(invList,"$(f_bus)")
            invConnected["$(f_bus)"] = "$(t_bus)"
            invLine["$(f_bus)"] = k
        else
            push!(invList,"$(t_bus)")
            invConnected["$(t_bus)"] = "$(f_bus)"
            invLine["$(t_bus)"] = k
        end
    end

    # obtain the list of buses with loads, assuming balanced loads
    loadList = Dict()
    vnomList = Dict()
    for i in keys(mpData["load"])
        println(mpData["load"][i]["load_bus"])
        if !("$(mpData["load"][i]["load_bus"])" in keys(loadList))
            loadList["$(mpData["load"][i]["load_bus"])"] = [[0.0,0.0,0.0],[0.0,0.0,0.0]]
            loadList["$(mpData["load"][i]["load_bus"])"][1] += mpData["load"][i]["pd"]
            loadList["$(mpData["load"][i]["load_bus"])"][2] += mpData["load"][i]["qd"]
            vnomList["$(mpData["load"][i]["load_bus"])"] = [0.0,0.0,0.0]
            vnomList["$(mpData["load"][i]["load_bus"])"] += (mpData["load"][i]["pd"] .!= 0.0)*mpData["load"][i]["vnom_kv"]
        else
            loadList["$(mpData["load"][i]["load_bus"])"][1] += mpData["load"][i]["pd"]
            vnomList["$(mpData["load"][i]["load_bus"])"] += (mpData["load"][i]["pd"] .!= 0.0)*mpData["load"][i]["vnom_kv"]
        end
    end

    # get the bus list
    busList = []
    for i in keys(mpData["bus"])
        if mpData["bus"][i]["bus_type"] != 3
            push!(busList, i)
        end
    end

    return busList, brList, invList, invConnected, invLine, loadList, vnomList
end


""
function procLoad(mpData, loadList, vnomList, omega0)
    θList = [0,2*pi/3,-2*pi/3]
    load_R = Dict()
    load_X = Dict()
    load_L = Dict()
    # Get the load resistance and reactance
    for i in keys(loadList)
        load_R[i] = zeros(3,3)
        load_X[i] = zeros(3,3)
        for j in 1:3
            load_R[i][j,j] = loadList[i][1][j]/vnomList[i][j]^2
            load_X[i][j,j] = loadList[i][2][j]/((vnomList[i][j]*cos(θList[j]))^2 - (vnomList[i][j]*sin(θList[j]))^2)
        end
        load_L[i] = load_X[i]./(omega0)
    end
    return load_L,load_R,load_X
end


"obtain the A matrix for inverter based generator"
function obtainA_inverter(mpData,opfSol,iBusList,omega0,mP,mQ,tau)
    A = Dict()

    errorList = []
    for i in iBusList
        if mpData["bus"][i]["bus_type"] == 3
            # there should be only one branch from the generator
            iInd = parse(Int64,i)
            link = [mpData["branch"][br] for br in keys(mpData["branch"]) if mpData["branch"][br]["f_bus"] == iInd | mpData["branch"][br]["t_bus"] == iInd][1]
            L = link["br_x"]./omega0
            delta0 = opfSol["solution"]["bus"][i]["va"][1]
            vd0 = opfSol["solution"]["bus"][i]["vm"].*cos.(opfSol["solution"]["bus"][i]["va"])
            vq0 = opfSol["solution"]["bus"][i]["vm"].*sin.(opfSol["solution"]["bus"][i]["va"])
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
                id0[j],iq0[j] = inv([vd0[j] vq0[j]vq0[j] -vd0[j]])*[P0/3Q0/3]
            end
            dvd_delta = -vq0
            dvq_delta = vd0
            dvd_v = cos.(opfSol["solution"]["bus"][i]["va"])
            dvq_v = sin.(opfSol["solution"]["bus"][i]["va"])
            dP_delta = dot(dvd_delta,id0) + dot(dvq_delta,iq0)
            dQ_delta = -dot(dvd_delta,iq0) + dot(dvq_delta,id0)
            dP_v = dot(dvd_v,id0) + dot(dvq_v,iq0)
            dQ_v = -dot(dvd_v,iq0) + dot(dvq_v,id0)
            dP_id = vd0
            dP_iq = vq0
            dQ_id = vq0
            dQ_iq = -vd0

            A[i] = zeros(9,9)
            A[i][1,2] = 1
            A[i][2,:] = -1/tau[i]*[mP[i]*dP_delta
                                1
                                mP[i]*dP_v
                                mP[i]*dP_id[1] mP[i]*dP_id[2] mP[i]*dP_id[3]
                                mP[i]*dP_iq[1] mP[i]*dP_iq[2] mP[i]*dP_iq[3]
                                ]
            A[i][3,:] = -1/tau[i]*[mQ[i]*dQ_delta
                                0
                                1+mQ[i]*dQ_v
                                mQ[i]*dQ_id[1] mQ[i]*dQ_id[2] mP[i]*dQ_id[3]
                                mQ[i]*dQ_iq[1] mQ[i]*dQ_iq[2] mP[i]*dQ_iq[3]
                                ]
            A[i][4:6,1] = inv(L)*dvd_delta
            A[i][4:6,3] = inv(L)*dvd_v
            A[i][4:6,4:6] = -inv(L)*link["br_r"]
            A[i][4:6,7:9] = inv(L)*link["br_x"]
            A[i][7:9,1] = inv(L)*dvq_delta
            A[i][7:9,3] = inv(L)*dvq_v
            A[i][7:9,4:6] = -inv(L)*link["br_x"]
            A[i][7:9,7:9] = -inv(L)*link["br_r"]
        else
            push!(errorList,i)
        end
    end
    return A,errorList
end


""
function obtainA_inverter_global(mpData, opfSol, omega0, mP, mQ, tau, rN, busList, invList, invConnected)
    A = Dict()

    for i in busList
        # initialize the matrix to be zero matrices, if it is a regular bus
        A[i] = zeros(9,9)
    end

    for busConnected in invList
        # there should be only one branch from the generator
        i = invConnected[busConnected]
        iInd = parse(Int64,i)
        link = mpData["branch"][invLine[busConnected]]
        L = link["br_x"]./omega0
        delta0 = opfSol["solution"]["bus"][i]["va"][1]
        vd0 = opfSol["solution"]["bus"][i]["vm"].*cos.(opfSol["solution"]["bus"][i]["va"])
        vq0 = opfSol["solution"]["bus"][i]["vm"].*sin.(opfSol["solution"]["bus"][i]["va"])
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
            id0[j],iq0[j] = inv([vd0[j] vq0[j]vq0[j] -vd0[j]])*[P0/3Q0/3]
        end
        dvd_delta = -vq0
        dvq_delta = vd0
        dvd_v = cos.(opfSol["solution"]["bus"][i]["va"])
        dvq_v = sin.(opfSol["solution"]["bus"][i]["va"])
        dP_delta = dot(dvd_delta,id0) + dot(dvq_delta,iq0)
        dQ_delta = -dot(dvd_delta,iq0) + dot(dvq_delta,id0)
        dP_v = dot(dvd_v,id0) + dot(dvq_v,iq0)
        dQ_v = -dot(dvd_v,iq0) + dot(dvq_v,id0)
        dP_id = vd0
        dP_iq = vq0
        dQ_id = vq0
        dQ_iq = -vd0

        A[busConnected][1,2] = 1
        A[busConnected][2,:] = -1/tau[i]*[mP[i]*dP_delta
                            1
                            mP[i]*dP_v
                            mP[i]*dP_id[1] mP[i]*dP_id[2] mP[i]*dP_id[3]
                            mP[i]*dP_iq[1] mP[i]*dP_iq[2] mP[i]*dP_iq[3]
                            ]
        A[busConnected][3,:] = -1/tau[i]*[mQ[i]*dQ_delta
                            0
                            1+mQ[i]*dQ_v
                            mQ[i]*dQ_id[1] mQ[i]*dQ_id[2] mP[i]*dQ_id[3]
                            mQ[i]*dQ_iq[1] mQ[i]*dQ_iq[2] mP[i]*dQ_iq[3]
                            ]
        A[busConnected][4:6,1] = inv(L)*dvd_delta
        A[busConnected][4:6,3] = inv(L)*dvd_v
        A[busConnected][4:6,4:6] = -inv(L)*link["br_r"] - rN*inv(L)
        A[busConnected][4:6,7:9] = inv(L)*link["br_x"]
        A[busConnected][7:9,1] = inv(L)*dvq_delta
        A[busConnected][7:9,3] = inv(L)*dvq_v
        A[busConnected][7:9,4:6] = -inv(L)*link["br_x"]
        A[busConnected][7:9,7:9] = -inv(L)*link["br_r"] - rN*inv(L)
    end
    return A
end


""
function obtainB_inverter_global(mpData, rN, omega0, busList, brList, invList, invLine)
    B = Dict()

    for i in busList
        for k in brList
            # initialize the matrix to be zero matrices
            B[i,k] = zeros(9,6)
            iInd = parse(Int64,i)
            if i in invList
                L = mpData["branch"][invLine[i]]["br_x"]./omega0
                if mpData["branch"][k]["f_bus"] == iInd
                    B[i,k][4:6,1:3] = rN*inv(L)
                    B[i,k][7:9,4:6] = rN*inv(L)
                elseif mpData["branch"][k]["t_bus"] == iInd
                    B[i,k][4:6,1:3] = -rN*inv(L)
                    B[i,k][7:9,4:6] = -rN*inv(L)
                end
            end
        end
    end
    return B
end


""
function obtainC_inverter_global(mpData, rN, omega0, busList, brList, invList, invLine, load_L)
    C = Dict()

    for i in busList
        # initialize the matrix to be zero matrices
        C[i] = zeros(9,6)
        iInd = parse(Int64,i)
        if (i in invList) & (i in keys(load_L))
            L = load_L[i]
            C[i][4:6,1:3] = rN*inv(L)
            C[i][7:9,4:6] = rN*inv(L)
        end
    end
    return C
end


""
function obtainD_inverter_global(mpData,rN, omega0, busList, brList, invList, invLine)
    D = Dict()

    for k in brList
        L = mpData["branch"][k]["br_x"]./omega0
        for i in busList
            # it is a regular bus
            D[k,i] = zeros(6,9)
            if i == "$(mpData["branch"][k]["f_bus"])"
                D[k,i][1:3,4:6] = rN*inv(L)
                D[k,i][4:6,7:9] = rN*inv(L)
            end
            if i == "$(mpData["branch"][k]["t_bus"])"
                D[k,i][1:3,4:6] = -rN*inv(L)
                D[k,i][4:6,7:9] = -rN*inv(L)
            end
        end
    end
    return D
end


""
function obtainE_inverter_global(mpData, rN, omega0, brList, invList, invLine)
    E = Dict()

    for k1 in brList
        L = mpData["branch"][k1]["br_x"]./omega0
        for k2 in brList
            E[k1,k2] = zeros(6,6)
            if k1 == k2
                E[k1,k2][1:3,1:3] = -2*rN*inv(L) - inv(L)*mpData["branch"][k1]["br_r"]
                E[k1,k2][1:3,4:6] = inv(L)*mpData["branch"][k1]["br_x"]
                E[k1,k2][4:6,1:3] = -inv(L)*mpData["branch"][k1]["br_x"]
                E[k1,k2][4:6,4:6] = -2*rN*inv(L) - inv(L)*mpData["branch"][k1]["br_r"]
            elseif (mpdata["branch"][k1]["f_bus"] == mpdata["branch"][k2]["f_bus"]) | (mpdata["branch"][k1]["t_bus"] == mpdata["branch"][k2]["t_bus"])
                E[k1,k2][1:3,1:3] = -rN*inv(L)
                E[k1,k2][4:6,4:6] = -rN*inv(L)
            elseif mpdata["branch"][k1]["f_bus"] == mpdata["branch"][k2]["t_bus"]
                if mpdata["branch"][k1]["t_bus"] == mpdata["branch"][k2]["f_bus"]
                    E[k1,k2][1:3,1:3] = 2*rN*inv(L)
                    E[k1,k2][4:6,4:6] = 2*rN*inv(L)
                else
                    E[k1,k2][1:3,1:3] = rN*inv(L)
                    E[k1,k2][4:6,4:6] = rN*inv(L)
                end
            end
        end
    end

    return E
end


""
function obtainF_inverter_global(mpData, rN, omega0, busList, brList, loadList, invList, invLine)
    F = Dict()

    for k in brList
        L = mpData["branch"][k]["br_x"]./omega0
        for i in busList
            F[k,i] = zeros(6,6)
            # obtain line L
            if i in keys(loadList)
                if mpData["branch"][k]["f_bus"] == parse(Int64,i)
                    F[k,i][1:3,1:3] = -rN*inv(L)
                    F[k,i][4:6,4:6] = -rN*inv(L)
                elseif mpData["branch"][k]["t_bus"] == parse(Int64,i)
                    F[k,i][1:3,1:3] = rN*inv(L)
                    F[k,i][4:6,4:6] = rN*inv(L)
                end
            end
        end
    end

    return F
end


""
function obtainG_inverter_global(mpData, rN, omega0, busList, brList, invList, invLine, loadList, load_L)
    G = Dict()

    for i in busList
        G[i] = zeros(6,9)
        if (i in invList) & (i in keys(loadList))
            L = load_L[i]
            G[1:3,4:6] = rN*inv(L)
            G[4:6,7:9] = rN*inv(L)
        end
    end

    return G
end


""
function obtainH_inverter_global(mpData, rN, omega0, busList, brList, invList, invLine, load_L)
    H = Dict()

    for i in busList
        for k in brList
            H[i,k] = zeros(6,6)
            # obtain line L
            if i in keys(load_L)
                L = load_L[i]
                try
                    # if L != 0
                    if mpData["branch"][k]["f_bus"] == parse(Int64,i)
                        H[i,k][1:3,1:3] = -rN*inv(L)
                        H[i,k][4:6,4:6] = -rN*inv(L)
                    elseif mpData["branch"][k]["t_bus"] == parse(Int64,i)
                        H[i,k][1:3,1:3] = rN*inv(L)
                        H[i,k][4:6,4:6] = rN*inv(L)
                    end
                catch
                    println("L[$(i)] is not invertable!")
                end
            end
        end
    end

    return H
end


""
function obtainI_inverter_global(mpData, rN, omega0, busList, brList, invList, invLine, load_L, load_R, load_X)
    I = Dict()

    for i in busList
        I[i] = zeros(6,6)
        if i in keys(loadList)
            L = load_L[i]
            try
                I[i][1:3,1:3] = -rN*inv(L) - inv(L)*load_R[i]
                I[i][1:3,4:6] = inv(L)*load_X[i]
                I[i][4:6,1:3] = -inv(L)*load_X[i]
                I[i][1:3,4:6] = -rN*inv(L) - inv(L)*load_R[i]
            catch
                println("L[$(i)] is not invertable!")
            end
        end
    end
    return I
end


""
function combineSub(busList, brList, A, B, C, D, E, F, G, H, I)
    n = length(busList)
    m = length(brList)
    Atot = zeros(15*n + 6*m, 15*n + 6*m)
    # Fill in the total matrix
    for i in 1:n
        # fill in A
        Atot[9*(i-1)+1:9*i,9*(i-1)+1:9*i] = A[busList[i]]
        # fill in C
        Atot[9*(i-1)+1:9*i,9*n+6*m+6*(i-1)+1:9*n+6*m+6*i] = C[busList[i]]
        # fill in G
        Atot[9*n+6*m+6*(i-1)+1:9*n+6*m+6*i,9*(i-1)+1:9*i] = G[busList[i]]
        # fill in I
        Atot[9*n+6*m+6*(i-1)+1:9*n+6*m+6*i,9*n+6*m+6*(i-1)+1:9*n+6*m+6*i] = I[busList[i]]
        for j in 1:m
            # fill in B
            Atot[9*(i-1)+1:9*i,9*n+6*(j-1)+1:9*n+6*j] = B[busList[i],brList[j]]
            # fill in D
            Atot[9*n+6*(j-1)+1:9*n+6*j,9*(i-1)+1:9*i] = D[brList[j],busList[i]]
            # fill in F
            Atot[9*n+6*(j-1)+1:9*n+6*j,9*n+6*m+6*(i-1)+1:9*n+6*m+6*i] = F[brList[j],busList[i]]
            # fill in H
            Atot[9*n+6*m+6*(i-1)+1:9*n+6*m+6*i,9*n+6*(j-1)+1:9*n+6*j] = H[busList[i],brList[j]]
        end
    end
    for j1 in 1:m
        for j2 in 1:m
            # fill in E
            Atot[9*n+6*(j1-1)+1:9*n+6*j1,9*n+6*(j2-1)+1:9*n+6*j2] = E[brList[j1],brList[j2]]
        end
    end
    return Atot
end


"obtain the global matrix of the"
function obtainGlobal(mpData,opfSol,omega0,mP,mQ,tau,rN)
    # preprocessing
    busList, brList, invList, invConnected, invLine, loadList, vnomList = preproc(mpData)
    load_L,load_R,load_X = procLoad(mpData, loadList, vnomList, omega0)

    # obtain A matrix
    Asub = obtainA_inverter_global(mpData, opfSol, omega0, mP, mQ, tau, rN, busList, invList, invConnected)
    Bsub = obtainB_inverter_global(mpData, rN, omega0, busList, brList, invList, invLine)
    Csub = obtainC_inverter_global(mpData, rN, omega0, busList, brList, invList, invLine, load_L)
    Dsub = obtainD_inverter_global(mpData, rN, omega0, busList, brList, invList, invLine)
    Esub = obtainE_inverter_global(mpData, rN, omega0, brList, invList, invLine)
    Fsub = obtainF_inverter_global(mpData, rN, omega0, busList, brList, loadList, invList, invLine)
    Gsub = obtainG_inverter_global(mpData, rN, omega0, busList, brList, invList, invLine, loadList, load_L)
    Hsub = obtainH_inverter_global(mpData, rN, omega0, busList, brList, invList, invLine, load_L)
    Isub = obtainI_inverter_global(mpData, rN, omega0, busList, brList, invList, invLine, load_L, load_R, load_X)

    Atot = combineSub(busList, brList, Asub, Bsub, Csub, Dsub, Esub, Fsub, Gsub, Hsub, Isub)

    return Atot
end
