"Modify the PowerModel opf model to initialize a relaxed model"
function initialize_relaxation(pm::PowerModels.AbstractPowerModel, nw::Int, Amg::Array{Any,2})
    # add a variable for matrix B and P
    n = size(Amg)[1]
    # set up P
    JuMP.@variable(pm.model, P[i in 1:n, j in 1:n], Symmetric)
    # diagonal terms are nonnegative for P
    JuMP.@constraint(pm.model, Ppositive[i in 1:n], P[i,i] >= 0)

    # set up B
    # JuMP.@variable(pm.model, B[1:n,1:n])
    # diagonal terms are nonpositive for B
    # JuMP.@constraint(pm.model, Bnegative[i in 1:n], B[i,i] <= 0)
    # JuMP.@constraint(pm.model, BSymmetric[i in 1:n, j in i:n], B[i,j] == B[j,i])
    Bform = Dict();
    for i in 1:n
        for j in 1:n
            Bform[i,j] = JuMP.@NLexpression(pm.model, sum(Amg[k,i]*pm.model[:P][k,j] + pm.model[:P][i,k]*Amg[k,j] for k in 1:n))
        end
    end
end

function initialize_relaxation(pm::PowerModels.AbstractPowerModel, nw::Int, Amg::Array{Any,2})
    # add a variable for matrix B and P
    n = size(Amg)[1]
    # set up P, only positive diagonal terms
    JuMP.@variable(pm.model, P[i in 1:n] >= 0);

end

"Obtain the P matrix from the OPF solution"
function obtain_P(pm::PowerModels.AbstractPowerModel, nw::Int, n::Int)
    P_value_list = JuMP.value.(pm.model[:P])
    P_value_m = zeros(n,n)
    for i in 1:n
        P_value_m[i,i] = P_value_list[i]
    end
    return P_value_m
end

"Obtain the B matrix from the OPF solution"
function obtain_B(pm::PowerModels.AbstractPowerModel, nw::Int, n::Int, A_value::Array{Any,2}, P_value::Array{Any,2})
    B_value = transpose(A_value)*(P_value) + P_value*A_value
    return B_value
end

"Force the B matrix to be symmetric"
function force_symmetric!(B_value::Array{Any,2})
    n = size(B_value)[1]
    for i in 1:(n-1)
        for j in (i+1):n
            B_value[j,i] = B_value[i,j]
        end
    end
end

"Make elements close to 0 to 0"
function force_zero(mat)
    n,m = size(mat)
    for i in 1:n
        for j in 1:n
            if abs(mat[i,j]) < 1e-6
                mat[i,j] = 0
            end
        end
    end
    return mat
end

"Generate SOCP constraints"
function socp_cons(pm::PowerModels.AbstractPowerModel, nw::Int, n::Int, ni::Int,Amg::Array{Any,2})
    # P nonnegative constraint
    JuMP.@constraint(pm.model, Pcons[i in 1:(n-1)], pm.model[:P][i,i]*pm.model[:P][i+1,i+1] - pm.model[:P][i,i+1]*pm.model[:P][i+1,i] >= 0.0)
    # B nonpositive constraint
    AnonzeroDict = Dict()
    for i in 1:n
        AnonzeroDict[i] = []
        for j in 1:n
            if Amg[j,i] != 0.0
                push!(AnonzeroDict[i],j)
            end
        end
    end
    JuMP.@NLconstraint(pm.model, Bcons[i in 1:(n-1)], (sum(Amg[j,i]*pm.model[:P][j,i] for j in 1:n if j in AnonzeroDict[i]) + sum(pm.model[:P][i,j]*Amg[j,i] for j in 1:n if j in AnonzeroDict[i]))*
        (sum(Amg[j,i+1]*pm.model[:P][j,i+1] for j in 1:n if j in AnonzeroDict[i+1]) + sum(pm.model[:P][i+1,j]*Amg[j,i+1] for j in 1:n if j in AnonzeroDict[i+1])) -
        (sum(Amg[j,i]*pm.model[:P][j,i+1] for j in 1:n if j in AnonzeroDict[i]) + sum(pm.model[:P][i,j]*Amg[j,i+1] for j in 1:n if j in AnonzeroDict[i+1]))*
        (sum(Amg[j,i+1]*pm.model[:P][j,i] for j in 1:n if j in AnonzeroDict[i+1]) + sum(pm.model[:P][i+1,j]*Amg[j,i] for j in 1:n if j in AnonzeroDict[i])) <= 0.0)
end

"Add the B-P cuts to the OPF problem"
function gen_cuts(pm::PowerModels.AbstractPowerModel, nw::Int, n::Int, ni::Int,
        B_value::Array{Any,2}, P_value::Array{Any,2}, Amg::Array{Any,2}, vB = [], vP = [], P_Pos = true)
    # obtain violated eigenvalue for B
    if vB == []
        # obtain the B matrix eigenvector/value
        BeigValList, BeigVectorList = eigen(B_value)
        force_symmetric!(B_value)
        for eigInd in 1:length(BeigValList)
            eig = BeigValList[eigInd]
            if eig > 0
                BeigVec = BeigVectorList[:,eigInd]
                for i in 1:length(BeigVec)
                    if abs(BeigVec[i]) < 1e-6
                        BeigVec[i] = 0
                    end
                end
                push!(vB, BeigVec)
            end
        end
    end
    # generate B cuts
    for item in vB
        nonzeroList = [i for i in 1:length(item) if item[i] != 0.0]
        nonzeroBool = false
        AnonzeroDict = Dict()
        for i in nonzeroList
            AnonzeroDict[i] = []
            for k in 1:n
                if Amg[k,i] != 0.0
                    nonzeroBool = true
                    push!(AnonzeroDict[i],k)
                end
            end
        end
        nonzeroA = []
        for i in nonzeroList
            if length(AnonzeroDict[i]) > 0
                for j in nonzeroList
                    push!(nonzeroA,(i,j))
                    push!(nonzeroA,(j,i))
                end
            end
        end
        nonzeroA = unique(nonzeroA)
        nonzeroA_new = []
        for i in nonzeroList
            for j in nonzeroList
                if Amg[i,j] != 0.0
                    push!(nonzeroA_new,(i,j))
                end
                if Amg[j,i] != 0.0
                    push!(nonzeroA_new,(j,i))
                end
            end
        end
        nonzeroA_new = unique(nonzeroA_new)
        if !(P_Pos)
            if nonzeroBool
                if minimum(nonzeroList) <= 9*ni
                    nlexp = JuMP.@NLexpression(pm.model, sum(item[i]*
                        (sum(Amg[k,i]*pm.model[:P][k,j] for k in AnonzeroDict[i]) + sum(pm.model[:P][i,k]*Amg[k,j] for k in AnonzeroDict[j]))
                        *item[j] for (i,j) in nonzeroA))
                    JuMP.@NLconstraint(pm.model, nlexp <= 0.0)
                else
                    JuMP.@constraint(pm.model, sum(item[i]*
                        (sum(Amg[k,i]*pm.model[:P][k,j] for k in AnonzeroDict[i]) + sum(pm.model[:P][i,k]*Amg[k,j] for k in AnonzeroDict[j]))
                        *item[j] for (i,j) in nonzeroA) <= 0.0)
                end
            end
            # generate P cuts
            if vP == []
                # obtain the P matrix eigenvector/value
                PeigValList, PeigVectorList = eigen(P_value)
                for eigInd in 1:length(PeigValList)
                    eig = PeigValList[eigInd]
                    if eig < 0
                        PeigVec = PeigVectorList[:,eigInd]
                        for i in 1:length(PeigVec)
                            if abs(PeigVec[i]) < 1e-6
                                PeigVec[i] = 0
                            end
                        end
                        append!(vP, PeigVec)
                    end
                end
            end
            # generate P cuts
            for item in vP
                nonzeroList = [i for i in 1:length(item) if item[i] != 0.0]
                JuMP.@constraint(pm.model, sum(sum(item[i]*pm.model[:P][i,j]*item[j] for j in nonzeroList) for i in nonzeroList) >= 0)
            end
        else
            if nonzeroBool
                if minimum(nonzeroList) <= 9*ni
                    nlexp = JuMP.@NLexpression(pm.model, sum(item[i]*(pm.model[:P][i]*Amg[i,j])*item[j] for (i,j) in nonzeroA if (i,j) in nonzeroA_new) +
                        sum(item[i]*(Amg[j,i]*pm.model[:P][j])*item[j] for (i,j) in nonzeroA if (j,i) in nonzeroA_new))
                    JuMP.@NLconstraint(pm.model, nlexp <= 0.0)
                else
                    JuMP.@constraint(pm.model, sum(item[i]*(Amg[j,i]*pm.model[:P][j] + pm.model[:P][i]*Amg[i,j])
                        *item[j] for (i,j) in nonzeroA) <= 0.0)
                end
            end
        end
    end
end
