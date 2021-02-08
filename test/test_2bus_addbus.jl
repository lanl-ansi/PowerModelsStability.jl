# test the two bus system
include("../src/PowerModelsStability.jl");
# load the data
filePath = "./data/case2_diag.dss";
# filePath = "../data/case5_phase_drop.dss";
const _PMS = PowerModelsStability;

# load the json file for stability-specific parameters
inverter_data = _PMS.parse_json("./data/case2_inverters.json");

# solve the opf problem
mpData = _PMS.PMD.parse_file(filePath);
# add the inverter buses
mpData = _PMS.add_inverters(inverter_data, mpData);
mpData_math = _PMS.transform_data_model(mpData);
ipopt_solver = with_optimizer(Ipopt.Optimizer, tol=1e-5, print_level=0);

# obtain the opf solution and the opf model
opfSol = _PMS.PMD.run_mc_opf(mpData, PowerModels.ACPPowerModel, ipopt_solver);
pm = _PMS.PMD.instantiate_mc_model(mpData,PowerModels.ACPPowerModel,PowerModelsDistribution.build_mc_opf);

# change the virtual bus to be an inverter bus for test
ω0 = inverter_data["omega0"];
rN = inverter_data["rN"];

Atot = obtainGlobal(mpData,opfSol,ω0,rN);
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
