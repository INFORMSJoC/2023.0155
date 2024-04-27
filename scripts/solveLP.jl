#solves the matrix game from the Proposition 4 test

#input B is the payoff matrix 

# NOTE: the row player is the searcher and the column player the hider 

# returns a hider optimal strategy, a searcher optimal strategy and the value

function solveLP(B::Array{Float64,2}, n::Int64)
    model = Model(GLPK.Optimizer)
    @variable(model, x[1:n] >= 0.)
    @objective(model, Min, sum(x))
    @constraint(model, con, B * x .>= 1)
    optimize!(model)
    @assert is_solved_and_feasible(model; dual = true)
    v = 1 / objective_value(model)
    hider_opt = value.(x) .* v
    searcher_opt = -dual.(con) .* v
    return hider_opt, searcher_opt, v
end
