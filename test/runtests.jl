using PowerModelsStability
using Test

import PowerModelsDistribution
import Ipopt
import LinearAlgebra: eigvals, eigvecs

const PMD = PowerModelsDistribution
const PMS = PowerModelsStability

PMD.silence!()

ipopt_solver = PMD.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 3, "tol"=>1e-5)

@testset "PowerModelsStability" begin
    include("test_2bus.jl")
end