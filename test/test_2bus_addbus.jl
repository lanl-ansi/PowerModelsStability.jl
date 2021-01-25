# test the two bus system
include("../src/io/preprocessing.jl");
include("../src/core/inverter.jl");
include("../src/core/inverter_new.jl");
include("../src/core/constraint.jl");

# load the data
filePath = "./data/case2_diag.dss";
# filePath = "../data/case5_phase_drop.dss";

# load the json file for stability-specific parameters
inverter_data = parse_json("./data/case2_inverters.json");

# solve the opf problem
mpData = PowerModelsDistribution.parse_file(filePath);
# add the inverter buses
add_inverters(inverter_data, mpData);

mpData_math = transform_data_model(mpData);
ipopt_solver = with_optimizer(Ipopt.Optimizer, tol=1e-5, print_level=0);

# obtain the opf solution and the opf model
opfSol = PowerModelsDistribution.run_mc_opf(mpData, PowerModels.ACPPowerModel, ipopt_solver);
pm = PowerModelsDistribution.instantiate_mc_model(mpData,PowerModels.ACPPowerModel,PowerModelsDistribution.build_mc_opf);

# change the virtual bus to be an inverter bus for test
mpData["bus"]["3"]["bus_type"] = 4;
mP = Dict();
mQ = Dict();
mP["3"] = 0.3;
mQ["3"] = 0.3;

τ = Dict();
τ["3"] =1e-4;
rN = 1000;
ω0 = 2*pi*60;
Atot = obtainGlobal(mpData,opfSol,ω0,mP,mQ,τ,rN);
eigValList = eigvals(Atot);
eigVectorList = eigvecs(Atot);
statusTemp = true;
for eigInd in 1:length(eigValList)
    eig = eigValList[eigInd];
    vioList = [];
    if eig.re > 0
        statusTemp = false;
        push!(vioList,eigVectorList[eigInd,:]);
    end
end
if statusTemp
    println("The current OPF solution is stable.");
else
    Amg = obtainGlobal_var(mpData,pm,ω0,mP,mQ,τ,rN);
    constraint_stability(pm, 0, vioList, Amg);
end
