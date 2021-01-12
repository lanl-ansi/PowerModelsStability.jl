# require to load inverter.jl first

function matrixMulti(mp, A, v)
    # given a matrix A with NL expression and an eigenvector v, calculate the multiplication v^\top A v
    n = length(v);
    multiExpression = @NLexpression(mp, sum(sum(v[i]*A[i,j]*v[j] for j in 1:n) for i in 1:n));
    return multiExpression;
end

# create constraints to append to the power models
function obtainGlobal_var(mpData,pm,ω0,mP,mQ,τ,rN)
    busList, brList, invList, invConnected, invLine, loadList, vnomList = preproc(mpData);
    load_L,load_R,load_X = procLoad(mpData, loadList, vnomList, ω0);

    # obtain A matrix
    Asub = obtainA_inverter_global_var(mpData, pm, ω0, mP, mQ, τ, rN, busList, invList, invLine, invConnected);
    Bsub = obtainB_inverter_global(mpData, rN, ω0, busList, brList, invList, invLine);
    Csub = obtainC_inverter_global(mpData, rN, ω0, busList, brList, invList, invLine, load_L);
    Dsub = obtainD_inverter_global(mpData, rN, ω0, busList, brList, invList, invLine);
    Esub = obtainE_inverter_global(mpData, rN, ω0, brList);
    Fsub = obtainF_inverter_global(mpData, rN, ω0, busList, brList, loadList);
    Gsub = obtainG_inverter_global(mpData, rN, ω0, busList, brList, invList, loadList, load_L);
    Hsub = obtainH_inverter_global(mpData, rN, ω0, busList, brList, load_L);
    Isub = obtainI_inverter_global(mpData, rN, ω0, busList, brList, loadList, load_L, load_R, load_X);

    Atot = combineSub(busList, brList, Asub, Bsub, Csub, Dsub, Esub, Fsub, Gsub, Hsub, Isub, 2);
    return Atot;
end

function constraint_stability(pm::PowerModels.AbstractPowerModel, nw::Int, eigenVectorList::Array{Any,1}, Amg::Array{Any,2})
    # create constraint generate by the eigenvalue analysis
    for ev in eigenVectorList
        evMulti_real = matrixMulti(pm.model, Amg, real(ev));
        evMulti_imag = matrixMulti(pm.model, Amg, imag(ev));
        JuMP.@NLconstraint(pm.model, evMulti_real + evMulti_imag <= 0);
    end
end
