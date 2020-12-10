# traditional generator ODE construction

# obtain the A matrix for diesel generator
# To-do: need to change the bus_type input to an external source
function obtainA_diesel(mpData,opfSol,iBusList,ω0,kD)
    # obtain the voltage solution
    T = 2/3*[1 cos(-2*pi/3) cos(2*pi/3);
             0 -sin(-2*pi/3) -sin(2*pi/3);
             1/2 1/2 1/2];
    vD = Dict();
    vQ = Dict();
    vd = Dict();
    vq = Dict();
    δ0 = Dict();
    v0 = Dict();
    id = Dict();
    iq = Dict();
    iD = Dict();
    iQ = Dict();
    A = Dict();
    errorList = [];
    for i in iBusList
        if mpData["bus"][i]["bus_type"] == 3
            δ0[i] = opfSol["solution"]["bus"][i]["va"][1];
            vabcComplex = opfSol["solution"]["bus"][i]["vm"].*cos.(opfSol["solution"]["bus"][i]["va"]) +
                im*opfSol["solution"]["bus"][i]["vm"].*sin.(opfSol["solution"]["bus"][i]["va"]);
            vDQComplex = T*vabcComplex;
            if abs(vDQComplex[1].re) > 1e-8
                vD[i] = vDQComplex[1].re;
            else
                vD[i] = 0;
            end
            if abs(vDQComplex[1].im) > 1e-8
                vQ[i] = vDQComplex[1].im;
            else
                vQ[i] = 0;
            end
            vd[i] = vD[i]*cos(δ0[i]) + vQ[i]*sin(δ0[i]);
            vq[i] = vQ[i]*cos(δ0[i]) - vD[i]*sin(δ0[i]);
            # if it is a generator bus, search for the matching generator
            pg[i] = 0;
            qg[i] = 0;
            for g in keys(mpData["gen"])
                if "$(mpData["gen"][g]["gen_bus"])" == i
                    pg[i] += sum(opfSol["solution"]["gen"][g]["pg"]);
                    qg[i] += sum(opfSol["solution"]["gen"][g]["qg"]);
                    # Question: is this consistent with the calculation from line flow?
                end
            end
            id[i] = (pg[i]*vd[i] + qg[i]*vq[i])/(vd[i]^2 + vq[i]^2);
            iq[i] = (pg[i]*vq[i] - qg[i]*vd[i])/(vd[i]^2 + vq[i]^2);

            A[i] = zeros(2,2);
            A[i][1,2] = ω0;
            A[i][2,2] = -kD[i];
            A[i][2,1] = -(2*vd[i]*vq[i])/(ω0*(vd[i]^2+vq[i]^2))*(pg[i] - qg[i]);
        else
            push!(errorList,i);
        end
    end

    return A,errorList;
end
