using JuMP, COPT, Random, COPT
using LinearAlgebra, SpecialMatrices


optimizer = optimizer_with_attributes(
    COPT.Optimizer,
    "FeasTol" => 1e-3,
    "DualTol" => 1e-3
)

n = 1000
Random.seed!(1)


A = Hilbert(n) + parse(Float64, ARGS[1]) * I
c = rand(Float64, n)
b = rand(Float64, n)

model = Model(COPT.Optimizer)
@variable(model, x[1:n] >= 0)
@variable(model, s[1:n] >= 0)
@constraint(model, A * x - s == b)

expr_dual = @expression(model, c' * x)

set_objective_function(model, expr_dual)
# set_objective_sense(model, MIN_SENSE)

# optimize!(model)
write_to_file(model, "hilbert-$n-$(ARGS[1]).mps")
@info LinearAlgebra.cond(A)
