#solves the matrix game from the Proposition 4 test

#input B is the payoff matrix 

# NOTE: the row player is the searcher and the column player the hider 

# returns a hider optimal strategy, a searcher optimal strategy and the value

function solveLP(B::Array{Float64,2}, n::Int64)
    s = size(B,1) #number of searcher strategies
    game2 = Model(GLPK.Optimizer)
    @variable(game2, 0. <= x[1:n])
    @objective(game2, Min, sum(x[i] for i=1:n))
    @constraint(game2, con, 1 .<= B * x)
    @constraint(game2,con2, 0 .<= x)

    optimize!(game2)
    v = 1/objective_value(game2)
    hider_opt = value.(x).*v
    has_duals(game2)
    searcher_opt = [-dual(con[i])*v for i=1:s]
    return(hider_opt, searcher_opt, v)
end

