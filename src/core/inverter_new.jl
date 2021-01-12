# inverter ODE construction
const _PMs = PowerModels;
const _PMD = PowerModelsDistribution;

function vectorMulti(mp, C, v)
    # given a vector with NL expression and a constant matrix C, calculate the multiplication Cv
    n = length(v);
    multiExpression = [@NLexpression(mp, sum(C[i,j]*v[j] for j in 1:n)) for i in 1:n];
    return multiExpression;
end

function obtainA_inverter_global(mpData, opfSol, ω0, mP, mQ, τ, rN, busList, invList, invLine, invConnected, inverters, invBusDict)
    A = Dict();

    for invItem in inverters
        # initialize the matrix to be zero matrices, if it is a regular bus
        A[invItem,invItem] = zeros(9,9);
    end

    # for every inverter
    for i in inverters
        # there should be only one branch from the generator
        iInd = parse(Int64,i);
        link = mpData["branch"][invLine[i]];
        L = link["br_x"]./ω0;
        δ0 = opfSol["solution"]["bus"][i]["va"][1];
        vd0 = opfSol["solution"]["bus"][i]["vm"].*cos.(opfSol["solution"]["bus"][i]["va"]);
        vq0 = opfSol["solution"]["bus"][i]["vm"].*sin.(opfSol["solution"]["bus"][i]["va"]);
        P0 = 0;
        Q0 = 0;
        for g in keys(mpData["gen"])
            if "$(mpData["gen"][g]["gen_bus"])" == i
                P0 += sum(opfSol["solution"]["gen"][g]["pg"]);
                Q0 += sum(opfSol["solution"]["gen"][g]["qg"]);
                # Question: is this consistent with the calculation from line flow?
            end
        end
        id0 = zeros(3);
        iq0 = zeros(3);
        for j in 1:3
            id0[j],iq0[j] = inv([vd0[j] vq0[j];vq0[j] -vd0[j]])*[P0/3;Q0/3];
        end
        dvd_δ = -vq0;
        dvq_δ = vd0;
        dvd_v = cos.(opfSol["solution"]["bus"][i]["va"]);
        dvq_v = sin.(opfSol["solution"]["bus"][i]["va"]);
        dP_δ = dot(dvd_δ,id0) + dot(dvq_δ,iq0);
        dQ_δ = -dot(dvd_δ,iq0) + dot(dvq_δ,id0);
        dP_v = dot(dvd_v,id0) + dot(dvq_v,iq0);
        dQ_v = -dot(dvd_v,iq0) + dot(dvq_v,id0);
        dP_id = vd0;
        dP_iq = vq0;
        dQ_id = vq0;
        dQ_iq = -vd0;

        A[i,i][1,2] = 1;
        A[i,i][2,:] = -1/τ[i]*[mP[i]*dP_δ;
                            1;
                            mP[i]*dP_v;
                            mP[i]*dP_id[1]; mP[i]*dP_id[2]; mP[i]*dP_id[3];
                            mP[i]*dP_iq[1]; mP[i]*dP_iq[2]; mP[i]*dP_iq[3]
                            ];
        A[i,i][3,:] = -1/τ[i]*[mQ[i]*dQ_δ;
                            0;
                            1+mQ[i]*dQ_v;
                            mQ[i]*dQ_id[1]; mQ[i]*dQ_id[2]; mP[i]*dQ_id[3];
                            mQ[i]*dQ_iq[1]; mQ[i]*dQ_iq[2]; mP[i]*dQ_iq[3]
                            ];
        A[i,i][4:6,1] = inv(L)*dvd_δ;
        A[i,i][4:6,3] = inv(L)*dvd_v;
        A[i,i][4:6,4:6] = -inv(L)*link["br_r"] - rN*inv(L);
        A[i,i][4:6,7:9] = inv(L)*link["br_x"];
        A[i,i][7:9,1] = inv(L)*dvq_δ;
        A[i,i][7:9,3] = inv(L)*dvq_v;
        A[i,i][7:9,4:6] = -inv(L)*link["br_x"];
        A[i,i][7:9,7:9] = -inv(L)*link["br_r"] - rN*inv(L);
    end

    # for every pair of inverters in the same bus
    for iBus in invList
        if length(invConnected[iBus]) > 1
            for i in 1:length(invConnected[iBus])
                for j in (i+1):length(invConnected[iBus])
                    A[i,j] = zeros(9,9);
                    link = mpData["branch"][invLine[i]];
                    L = link["br_x"]./ω0;
                    A[i,j][4:6,4:6] = -inv(L)*link["br_r"] - rN*inv(L);
                    A[i,j][7:9,7:9] = -inv(L)*link["br_r"] - rN*inv(L);

                    link = mpData["branch"][invLine[j]];
                    L = link["br_x"]./ω0;
                    A[j,i] = zeros(9,9);
                    A[j,i][4:6,4:6] = -inv(L)*link["br_r"] - rN*inv(L);
                    A[j,i][7:9,7:9] = -inv(L)*link["br_r"] - rN*inv(L);
                end
            end
        end
    end

    return A;
end

function obtainA_inverter_global_var(mpData, pm, ω0, mP, mQ, τ, rN, busList, invList, invLine, invConnected, inverters, invBusDict)
    A = Dict();
    for invItem in inverters
        A[invItem,invItem] = Array{Any,2}(zeros(9,9));
    end

    for iBus in inverters
        # there should be only one branch from the generator
        iInd = parse(Int64,iBus);
        link = mpData["branch"][invLine[iBus]];
        L = link["br_x"]./ω0;
        i = mpData["bus"][iBus]["index"];
        δ0 = _PMs.var(pm,0,:va,i)[1];
        vmList = _PMs.var(pm,0,:vm,i);
        vaList = _PMs.var(pm,0,:va,i);
        vd0 = [@NLexpression(pm.model,vmList.data[j]*cos(vaList.data[j])) for j in 1:length(vmList.data)];
        vq0 = [@NLexpression(pm.model,vmList.data[j]*sin(vaList.data[j])) for j in 1:length(vmList.data)];
        P0 = @NLexpression(pm.model,0);
        Q0 = @NLexpression(pm.model,0);
        for g in keys(mpData["gen"])
            if "$(mpData["gen"][g]["gen_bus"])" == i
                P0 += sum(_PMs.var(pm,0,:pg,g));
                Q0 += sum(_PMs.var(pm,0,:qg,g));
            end
        end
        id0 = Array{Any,1}([0,0,0]);
        iq0 = Array{Any,1}([0,0,0]);
        for j in 1:3
            id0[j] = @NLexpression(pm.model,(vd0[j]*P0/3 + vq0[j]*Q0/3)/(vmList.data[j]^2));
            iq0[j] = @NLexpression(pm.model,(vq0[j]*P0/3 - vd0[j]*Q0/3)/(vmList.data[j]^2));
        end
        dvd_δ = [@NLexpression(pm.model,-vq0[j]) for j in 1:length(vq0)];
        dvq_δ = [@NLexpression(pm.model,vd0[j]) for j in 1:length(vd0)];
        dvd_v = [@NLexpression(pm.model,cos(vaList.data[j])) for j in 1:length(vaList.data)];
        dvq_v = [@NLexpression(pm.model,sin(vaList.data[j])) for j in 1:length(vaList.data)];
        dP_δ = @NLexpression(pm.model,sum(dvd_δ[j]*id0[j] for j in 1:length(dvd_δ)) + sum(dvq_δ[j]*iq0[j] for j in 1:length(dvq_δ)));
        dQ_δ = @NLexpression(pm.model,-sum(dvd_δ[j]*iq0[j] for j in 1:length(dvd_δ)) + sum(dvq_δ[j]*id0[j] for j in 1:length(dvq_δ)));
        dP_v = @NLexpression(pm.model,sum(dvd_v[j]*id0[j] for j in 1:length(dvd_v)) + sum(dvq_v[j]*iq0[j] for j in 1:length(dvq_v)));
        dQ_v = @NLexpression(pm.model,-sum(dvd_v[j]*iq0[j] for j in 1:length(dvd_v)) + sum(dvq_v[j]*id0[j] for j in 1:length(dvq_v)));
        dP_id = vd0;
        dP_iq = vq0;
        dQ_id = vq0;
        dQ_iq = [@NLexpression(pm.model,-vd0[i]) for i in 1:length(vd0)];

        A[iBus,iBus][1,2] = 1;
        A[iBus,iBus][2,:] = [@NLexpression(pm.model,-1/τ[iBus]*mP[iBus]*dP_δ);
                            -1/τ[iBus];
                            @NLexpression(pm.model,-1/τ[iBus]*mP[iBus]*dP_v);
                            @NLexpression(pm.model,-1/τ[iBus]*mP[iBus]*dP_id[1]);
                            @NLexpression(pm.model,-1/τ[iBus]*mP[iBus]*dP_id[2]);
                            @NLexpression(pm.model,-1/τ[iBus]*mP[iBus]*dP_id[3]);
                            @NLexpression(pm.model,-1/τ[iBus]*mP[iBus]*dP_iq[1]);
                            @NLexpression(pm.model,-1/τ[iBus]*mP[iBus]*dP_iq[2]);
                            @NLexpression(pm.model,-1/τ[iBus]*mP[iBus]*dP_iq[3])
                            ];
        A[iBus,iBus][3,:] = [@NLexpression(pm.model,-1/τ[iBus]*mQ[iBus]*dQ_δ);
                            0;
                            @NLexpression(pm.model,-1/τ[iBus]*(1+mQ[iBus]*dQ_v));
                            @NLexpression(pm.model,-1/τ[iBus]*(mQ[iBus]*dQ_id[1]));
                            @NLexpression(pm.model,-1/τ[iBus]*(mQ[iBus]*dQ_id[2]));
                            @NLexpression(pm.model,-1/τ[iBus]*(mP[iBus]*dQ_id[3]));
                            @NLexpression(pm.model,-1/τ[iBus]*(mQ[iBus]*dQ_iq[1]));
                            @NLexpression(pm.model,-1/τ[iBus]*(mQ[iBus]*dQ_iq[2]));
                            @NLexpression(pm.model,-1/τ[iBus]*(mP[iBus]*dQ_iq[3]))
                            ];
        A[iBus,iBus][4:6,1] = vectorMulti(pm.model,inv(L),dvd_δ);
        A[iBus,iBus][4:6,3] = vectorMulti(pm.model,inv(L),dvd_v);
        A[iBus,iBus][4:6,4:6] = -inv(L)*link["br_r"] - rN*inv(L);
        A[iBus,iBus][4:6,7:9] = inv(L)*link["br_x"];
        A[iBus,iBus][7:9,1] = vectorMulti(pm.model,inv(L),dvq_δ);
        A[iBus,iBus][7:9,3] = vectorMulti(pm.model,inv(L),dvq_v);
        A[iBus,iBus][7:9,4:6] = -inv(L)*link["br_x"];
        A[iBus,iBus][7:9,7:9] = -inv(L)*link["br_r"] - rN*inv(L);
    end

    # for every pair of inverters in the same bus
    for iBus in invList
        if length(invConnected[iBus]) > 1
            for i in 1:length(invConnected[iBus])
                for j in (i+1):length(invConnected[iBus])
                    A[i,j] = zeros(9,9);
                    link = mpData["branch"][invLine[i]];
                    L = link["br_x"]./ω0;
                    A[i,j][4:6,4:6] = -inv(L)*link["br_r"] - rN*inv(L);
                    A[i,j][7:9,7:9] = -inv(L)*link["br_r"] - rN*inv(L);

                    link = mpData["branch"][invLine[j]];
                    L = link["br_x"]./ω0;
                    A[j,i] = zeros(9,9);
                    A[j,i][4:6,4:6] = -inv(L)*link["br_r"] - rN*inv(L);
                    A[j,i][7:9,7:9] = -inv(L)*link["br_r"] - rN*inv(L);
                end
            end
        end
    end

    return A;
end

function obtainB_inverter_global(mpData, rN, ω0, busList, brList, invList, invLine, invConnected, inverters, invBusDict)
    B = Dict();

    for i in inverters
        iBus = invBusDict[i];
        for k in brList
            # initialize the matrix to be zero matrices
            B[i,k] = zeros(9,6);
            iInd = parse(Int64,iBus);
            L = mpData["branch"][invLine[i]]["br_x"]./ω0;
            if mpData["branch"][k]["f_bus"] == iInd
                B[i,k][4:6,1:3] = rN*inv(L);
                B[i,k][7:9,4:6] = rN*inv(L);
            elseif mpData["branch"][k]["t_bus"] == iInd
                B[i,k][4:6,1:3] = -rN*inv(L);
                B[i,k][7:9,4:6] = -rN*inv(L);
            end
        end
    end
    return B;
end

function obtainC_inverter_global(mpData, rN, ω0, busList, brList, invList, invLine, invConnected, load_L, inverters, invBusDict)
    C = Dict();

    for i in inverters
        # initialize the matrix to be zero matrices
        C[i,invBusDict[i]] = zeros(9,6);
        iInd = parse(Int64,i);
        if invBusDict[i] in keys(load_L)
            L = load_L[invBusDict[i]];
            C[i,invBusDict[i]][4:6,1:3] = rN*inv(L);
            C[i,invBusDict[i]][7:9,4:6] = rN*inv(L);
        end
    end
    return C;
end

function obtainD_inverter_global(mpData,rN, ω0, busList, brList, invList, invLine, inverters, invBusDict)
    D = Dict();

    for k in brList
        L = mpData["branch"][k]["br_x"]./ω0;
        for i in inverters
            iBus = invBusDict[i];
            # it is a regular bus
            D[k,i] = zeros(6,9);
            if iBus == "$(mpData["branch"][k]["f_bus"])"
                D[k,i][1:3,4:6] = rN*inv(L);
                D[k,i][4:6,7:9] = rN*inv(L);
            end
            if iBus == "$(mpData["branch"][k]["t_bus"])"
                D[k,i][1:3,4:6] = -rN*inv(L);
                D[k,i][4:6,7:9] = -rN*inv(L);
            end
        end
    end
    return D;
end

function obtainE_inverter_global(mpData, rN, ω0, brList)
    E = Dict();

    for k1 in brList
        L = mpData["branch"][k1]["br_x"]./ω0;
        for k2 in brList
            E[k1,k2] = zeros(6,6);
            if k1 == k2
                E[k1,k2][1:3,1:3] = -2*rN*inv(L) - inv(L)*mpData["branch"][k1]["br_r"];
                E[k1,k2][1:3,4:6] = inv(L)*mpData["branch"][k1]["br_x"];
                E[k1,k2][4:6,1:3] = -inv(L)*mpData["branch"][k1]["br_x"];
                E[k1,k2][4:6,4:6] = -2*rN*inv(L) - inv(L)*mpData["branch"][k1]["br_r"];
            elseif (mpData["branch"][k1]["f_bus"] == mpData["branch"][k2]["f_bus"]) | (mpData["branch"][k1]["t_bus"] == mpData["branch"][k2]["t_bus"])
                E[k1,k2][1:3,1:3] = -rN*inv(L);
                E[k1,k2][4:6,4:6] = -rN*inv(L);
            elseif mpData["branch"][k1]["f_bus"] == mpData["branch"][k2]["t_bus"]
                if mpData["branch"][k1]["t_bus"] == mpData["branch"][k2]["f_bus"]
                    E[k1,k2][1:3,1:3] = 2*rN*inv(L);
                    E[k1,k2][4:6,4:6] = 2*rN*inv(L);
                else
                    E[k1,k2][1:3,1:3] = rN*inv(L);
                    E[k1,k2][4:6,4:6] = rN*inv(L);
                end
            end
        end
    end

    return E;
end

function obtainF_inverter_global(mpData, rN, ω0, busList, brList, loadList)
    F = Dict();

    for k in brList
        L = mpData["branch"][k]["br_x"]./ω0;
        for i in busList
            F[k,i] = zeros(6,6);
            # obtain line L
            if i in keys(loadList)
                if mpData["branch"][k]["f_bus"] == parse(Int64,i)
                    F[k,i][1:3,1:3] = -rN*inv(L);
                    F[k,i][4:6,4:6] = -rN*inv(L);
                elseif mpData["branch"][k]["t_bus"] == parse(Int64,i)
                    F[k,i][1:3,1:3] = rN*inv(L);
                    F[k,i][4:6,4:6] = rN*inv(L);
                end
            end
        end
    end

    return F;
end

function obtainG_inverter_global(mpData, rN, ω0, busList, brList, invList, loadList, load_L, invLine, inverters, invBusDict)
    G = Dict();

    for i in inverters
        iBus = invBusDict[i];
        G[iBus,i] = zeros(6,9);
        if iBus in keys(loadList)
            L = load_L[iBus];
            G[iBus,i][1:3,4:6] = rN*inv(L);
            G[iBus,i][4:6,7:9] = rN*inv(L);
        end
    end

    return G;
end

function obtainH_inverter_global(mpData, rN, ω0, busList, brList, load_L)
    H = Dict();

    for i in busList
        for k in brList
            H[i,k] = zeros(6,6);
            # obtain line L
            if i in keys(load_L)
                L = load_L[i];
                try
                    # if L != 0
                    if mpData["branch"][k]["f_bus"] == parse(Int64,i)
                        H[i,k][1:3,1:3] = -rN*inv(L);
                        H[i,k][4:6,4:6] = -rN*inv(L);
                    elseif mpData["branch"][k]["t_bus"] == parse(Int64,i)
                        H[i,k][1:3,1:3] = rN*inv(L);
                        H[i,k][4:6,4:6] = rN*inv(L);
                    end
                catch
                    println("L[$(i)] is not invertable!")
                end
            end
        end
    end

    return H;
end

function obtainI_inverter_global(mpData, rN, ω0, busList, brList, loadList, load_L, load_R, load_X)
    I = Dict();

    for i in busList
        I[i] = zeros(6,6);
        if i in keys(loadList)
            L = load_L[i];
            try
                I[i][1:3,1:3] = -rN*inv(L) - inv(L)*load_R[i];
                I[i][1:3,4:6] = inv(L)*load_X[i];
                I[i][4:6,1:3] = -inv(L)*load_X[i];
                I[i][1:3,4:6] = -rN*inv(L) - inv(L)*load_R[i];
            catch
                println("L[$(i)] is not invertable!")
            end
        end
    end
    return I;
end

function combineSub(busList, brList, inverters, invBusDict, A, B, C, D, E, F, G, H, I, type)
    n = length(busList);
    ni = length(inverters);
    m = length(brList);
    if type == 1
        Atot = zeros(9*ni + 6*m + 6*n, 9*ni + 6*m + 6*n);
    else
        Atot = Array{Any,2}(zeros(9*ni + 6*m + 6*n, 9*ni + 6*m + 6*n));
    end
    # Fill in the total matrix
    for i in 1:ni
        # fill in A
        for j in 1:ni
            if (inverters[i],inverters[j]) in keys(A)
                Atot[9*(i-1)+1:9*i,9*(j-1)+1:9*j] = A[inverters[i],inverters[j]];
            end
        end
        for j in 1:m
            # fill in B
            Atot[9*(i-1)+1:9*i,9*ni+6*(j-1)+1:9*ni+6*j] = B[inverters[i],brList[j]];
            # fill in D
            Atot[9*ni+6*(j-1)+1:9*ni+6*j,9*(i-1)+1:9*i] = D[brList[j],inverters[i]];
        end

        for j in 1:n
            if busList[j] == invBusDict[inverters[i]]
                # fill in C
                Atot[9*(i-1)+1:9*i,9*ni+6*m+6*(j-1)+1:9*ni+6*m+6*j] = C[inverters[i],busList[j]];
                # fill in G
                Atot[9*ni+6*m+6*(j-1)+1:9*ni+6*m+6*j,9*(i-1)+1:9*i] = G[busList[j],inverters[i]];
            end
        end
    end

    for i in 1:n
        # fill in I
        Atot[9*ni+6*m+6*(i-1)+1:9*ni+6*m+6*i,9*ni+6*m+6*(i-1)+1:9*ni+6*m+6*i] = I[busList[i]];
        for j in 1:m
            # fill in F
            Atot[9*ni+6*(j-1)+1:9*ni+6*j,9*ni+6*m+6*(i-1)+1:9*ni+6*m+6*i] = F[brList[j],busList[i]];
            # fill in H
            Atot[9*ni+6*m+6*(i-1)+1:9*ni+6*m+6*i,9*ni+6*(j-1)+1:9*ni+6*j] = H[busList[i],brList[j]];
        end
    end

    for j1 in 1:m
        for j2 in 1:m
            # fill in E
            Atot[9*ni+6*(j1-1)+1:9*ni+6*j1,9*ni+6*(j2-1)+1:9*ni+6*j2] = E[brList[j1],brList[j2]];
        end
    end
    return Atot;
end

# obtain the global matrix of the
function obtainGlobal_multi(mpData,opfSol,ω0,mP,mQ,τ,rN)
    # preprocessing
    busList, brList, invList, invConnected, invLine, loadList, vnomList = preproc(mpData);
    load_L,load_R,load_X = procLoad(mpData, loadList, vnomList, ω0);

    inverters = [];
    invBusDict = Dict();

    for i in invList
        currentInv = invConnected[i];
        # obtain a list of inverter buses
        append!(inverters, currentInv);
        for invItem in currentInv
            invBusDict[invItem] = i;
        end
    end

    # obtain A matrix
    Asub = obtainA_inverter_global(mpData, opfSol, ω0, mP, mQ, τ, rN, busList, invList, invLine, invConnected, inverters, invBusDict);
    Bsub = obtainB_inverter_global(mpData, rN, ω0, busList, brList, invList, invLine, invConnected, inverters, invBusDict);
    Csub = obtainC_inverter_global(mpData, rN, ω0, busList, brList, invList, invLine, invConnected, load_L, inverters, invBusDict);
    Dsub = obtainD_inverter_global(mpData, rN, ω0, busList, brList, invList, invLine, inverters, invBusDict);
    Esub = obtainE_inverter_global(mpData, rN, ω0, brList);
    Fsub = obtainF_inverter_global(mpData, rN, ω0, busList, brList, loadList);
    Gsub = obtainG_inverter_global(mpData, rN, ω0, busList, brList, invList, loadList, load_L, invLine, inverters, invBusDict);
    Hsub = obtainH_inverter_global(mpData, rN, ω0, busList, brList, load_L);
    Isub = obtainI_inverter_global(mpData, rN, ω0, busList, brList, loadList, load_L, load_R, load_X);

    Atot = combineSub(busList, brList, inverters, invBusDict, Asub, Bsub, Csub, Dsub, Esub, Fsub, Gsub, Hsub, Isub, 1);

    return Atot;
end
