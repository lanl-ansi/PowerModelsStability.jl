# test the two bus system
using PowerModels,PowerModelsDistribution,InfrastructureModels,LinearAlgebra;
using JuMP,Ipopt,Gurobi;
using JSON;
include("../src/io.jl");
include("../src/inverter.jl");

# load the data
filePath = "./data/case2_diag.dss";
# filePath = "../data/case5_phase_drop.dss";

# load the json file for stability-specific parameters

# solve the opf problem
mpData = PowerModelsDistribution.parse_file(filePath);
mpData = transform_data_model(mpData);
ipopt_solver = with_optimizer(Ipopt.Optimizer, tol=1e-5, print_level=0);
opfSol = PowerModelsDistribution.run_mc_opf(mpData, PowerModels.ACPPowerModel, ipopt_solver);

# change the virtual bus to be an inverter bus for test
mpData["bus"]["3"]["bus_type"] = 4;
mP = Dict();
mQ = Dict();
mP["3"] = 0.3;
mQ["3"] = 0.3;
τ = Dict();
τ["3"] =1e-4;
rN = 1000;
Atot = obtainGlobal(mpData,opfSol,ω0,mP,mQ,τ,rN);
eigValList = eigvals(Atot);
statusTemp = true;
for eig in eigValList
    if eig.re > 0
        statusTemp = false;
    end
end
if statusTemp
    println("The current OPF solution is stable.");
end
