include("../src/PowerModelsStability.jl")
import JuMP
import PowerModels
import PowerModelsDistribution
const PMD = PowerModelsDistribution
import Ipopt
ipopt_solver = PMD.optimizer_with_attributes(Ipopt.Optimizer, "tol" => 1e-5, "linear_solver" => "ma27")
using LinearAlgebra

# Asub = PowerModelsStability.obtainA_inverter_global(mpData_math, opfSol, rN, omega0, busList, invList, invLine, invConnected, inverters, invBusDict)
# Bsub = PowerModelsStability.obtainB_inverter_global(mpData_math, rN, omega0, busList, brList, invList, invLine, invConnected, inverters, invBusDict)
# Csub = PowerModelsStability.obtainC_inverter_global(mpData_math, rN, omega0, busList, brList, invList, invLine, invConnected, load_L, loadConnections, inverters, invBusDict)
# Dsub = PowerModelsStability.obtainD_inverter_global(mpData_math, rN, omega0, busList, brList, invList, invLine, inverters, invBusDict)
# Esub = PowerModelsStability.obtainE_inverter_global(mpData_math, rN, omega0, brList)
# Fsub = PowerModelsStability.obtainF_inverter_global(mpData_math, rN, omega0, busList, brList, loadList, loadConnections)
# Gsub = PowerModelsStability.obtainG_inverter_global(mpData_math, rN, omega0, busList, brList, invList, loadList, load_L, loadConnections, invLine, inverters, invBusDict)
# Hsub = PowerModelsStability.obtainH_inverter_global(mpData_math, rN, omega0, busList, brList, load_L, loadConnections)
# Isub = PowerModelsStability.obtainI_inverter_global(mpData_math, rN, omega0, busList, brList, loadList, load_L, load_R, load_X, loadConnections)

filePath = "./data/ieee13/IEEE13Nodeckt_David_1.dss"
inverter_data = PowerModelsStability.parse_json("./data/case13_inverters.json")
inverter_data_math = Dict()
inverter_data_math["12"] = Dict()
inverter_data_math["12"]["tau"] = inverter_data["inverters"][1]["tau"]
inverter_data_math["12"]["mp"] = inverter_data["inverters"][1]["mp"]
inverter_data_math["12"]["mq"] = inverter_data["inverters"][1]["mq"]

mpData = PMD.parse_file(filePath)
opfSol,mpData_math = PowerModelsStability.run_mc_opf(mpData, PMD.ACRPowerModel, ipopt_solver; solution_processors=[PMD.sol_data_model!])
mpData_math = PowerModelsStability.convert_gen(mpData_math,inverter_data_math)
opfSol,mpData_math = PowerModelsStability.run_mc_opf(mpData_math, PMD.ACRPowerModel, ipopt_solver; solution_processors=[PMD.sol_data_model!])

pm = PMD.instantiate_mc_model(mpData_math, PMD.ACRPowerModel, PMD.build_mc_opf;
    ref_extensions=[PMD.ref_add_arcs_transformer!, PMD.ref_add_connections!])
opfSol_pm = PowerModels.optimize_model!(pm,optimizer = ipopt_solver;solution_processors=[PMD.sol_data_model!])

# change the virtual bus to be an inverter bus for test
omega0 = inverter_data["omega0"]
rN = inverter_data["rN"]
invData = PowerModelsStability.invData_proc(mpData_math, omega0)

P_Pos = true
ni = 1

Amg = PowerModelsStability.obtainGlobal_var(mpData_math, pm, rN, omega0,invData)
n = size(Amg)[1]
# set up P, only positive diagonal terms
JuMP.@variable(pm.model, P[i in 1:n] >= 0);
JuMP.@constraint(pm.model, consP[i in (ni*9+1):n], P[i] == 0)
#JuMP.@constraint(pm.model, sumP, sum(P[i] for i in 1:n) >= 1e-4)
JuMP.@constraint(pm.model, sumP, sum(P[i] for i in 1:(ni*9)) >= 1e-3)
opfSol_pm = PowerModels.optimize_model!(pm,optimizer = ipopt_solver;solution_processors=[PMD.sol_data_model!])

P_value = obtain_P(pm, 0, n)
A_value = PowerModelsStability.obtainGlobal_multi(mpData_math,opfSol_pm,rN,omega0,invData)
# A_value = force_zero(A_value)
B_value = transpose(A_value)*(P_value) + P_value*A_value
# add a step of cleaning small values of A & B !!!
# force small P elements to take 0 !!!
eigValList, eigVectorList = eigen(A_value)
statusTemp = true
vioList = []
for eigInd in 1:length(eigValList)
    eig = eigValList[eigInd]
    if eig.re > 0
        statusTemp = false
        push!(vioList,eigVectorList[:,eigInd])
    end
end

vB = []
BeigValList, BeigVectorList = eigen(B_value)
for eigInd in 1:length(BeigValList)
    eig = BeigValList[eigInd]
    if eig > 0
        BeigVec = BeigVectorList[:,eigInd]
        for i in 1:length(BeigVec)
            if abs(BeigVec[i]) < 1e-12
                BeigVec[i] = 0
            end
        end
        push!(vB, BeigVec)
    end
end

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
