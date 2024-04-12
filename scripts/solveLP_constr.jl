#solves the constrained matrix game from Step 3 of Algorithm 10

#input B is the payoff matrix 

# NOTE: the row player is the searcher and the column player the hider 

# returns a hider optimal strategy, a searcher optimal strategy and the value
function solveLP_constr(B::Array{Float64,2}, n::Int64, bounds::Array{Float64,1})
    s = size(B,1) #number of searcher strategies
    game = Model(GLPK.Optimizer)
    @variable(game, 0 <= p[1:n])
    @variable(game,0. <= v)
    @objective(game, Max, v)
    @constraint(game, con, v .<= B * p)
    @constraint(game, con2, sum(p) == 1)
    @constraint(game, con3, bounds .<= p)

    optimize!(game)
    hider_opt = value.(p)
    has_duals(game)
    searcher_opt = [-dual(con[i]) for i=1:s]
    return(hider_opt, searcher_opt, value.(v))
end



