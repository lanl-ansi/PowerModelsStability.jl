using PowerModelsStability
using Test

import PowerModelsDistribution
import Ipopt
import LinearAlgebra: eigvals, eigvecs

const PMD = PowerModelsDistribution
const PMS = PowerModelsStability

PMD.silence!()

ipopt_solver = PMD.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 3, "tol"=>1e-5, "max_iter" => Int(1E4), "sb" => "yes")

@testset "PowerModelsStability" begin
    include("test_2bus.jl")
end