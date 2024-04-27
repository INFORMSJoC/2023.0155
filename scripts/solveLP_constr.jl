#solves the constrained matrix game from Step 3 of Algorithm 10

#input B is the payoff matrix 

# NOTE: the row player is the searcher and the column player the hider 

# returns a hider optimal strategy, a searcher optimal strategy and the value
function solveLP_constr(B::Array{Float64,2}, n::Int64, bounds::Array{Float64,1})
    game = Model(GLPK.Optimizer)
    @variable(game, p[i in 1:n] >= bounds[i])
    @variable(game, v >= 0)
    @objective(game, Max, v)
    @constraint(game, con, v .<= B * p)
    @constraint(game, sum(p) == 1)
    optimize!(game)
    @assert is_solved_and_feasible(game; dual = true)
    hider_opt = value.(p)
    searcher_opt = -dual.(con)
    return(hider_opt, searcher_opt, value(v))
end


