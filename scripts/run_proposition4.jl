# for a search game specified by detection probabilities qq and search times cc, 
# determines whether p0 is optimal or not by using the test in Proposition 4
# the test compares two computed values. the input δp0 is used to determine how close these values have to be to be classed as "equal" (due to numerical errors in computation)
# see the main file for more details on the use of δp0

# the function returns 1 if p0 is optimal, 0 otherwise

function run_proposition4(qq::Array{Float64,1}, cc::Array{Float64,1}, n::Int64, accuracy::Int64, δp0::Float64)
    
    #creates the game matrix for the search game G_D (using the notation of Proposition 4)
    #also calculates the expected search time e if the hider plays p0 and the searcher optimally counters
    create = createLPmatp0(qq, cc, n, accuracy, true, false, [[1]])

    #solves the search game G_D
    opts = solveLP(create[1], n)
    
    #to return, indicator if optimal or not
    is_opt = 1

    #set is_opt to 0 if the value of G_D (opts[3]) is sufficiently (determined by user input δp0) from the expected search time e if the hider plays p0 and the searcher optimally counters (create[3])
    if abs(opts[3] - create[3])/opts[3] > δp0
        is_opt = 0
    end
    return(is_opt)
end
