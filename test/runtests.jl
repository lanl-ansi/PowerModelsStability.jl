using PowerModelsStability

import PowerModelsDistribution
import Ipopt

const PMD = PowerModelsDistribution
PMD._PM.silence()

import LinearAlgebra: eigvals, eigvecs

using Test

ipopt_solver = PMD.optimizer_with_attributes(Ipopt.Optimizer, "tol" => 1e-5, "print_level" => 0)

@testset "PowerModelsStability" begin
    include("test_2bus.jl")
end
