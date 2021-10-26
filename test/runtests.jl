using PowerModelsStability

import PowerModelsDistribution
import Ipopt

const PMD = PowerModelsDistribution
PMD.silence!()

import LinearAlgebra: eigvals, eigvecs

using Test

ipopt_solver = PMD.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)

@testset "PowerModelsStability" begin
    include("test_2bus.jl")
end
