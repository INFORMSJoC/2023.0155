# for a search game specified by detection probabilities qq and search times cc, 
# runs algorithm 10 with the bounds in bounds, and parameters ϵ, δ and max_iters

#returns 7 : value, expected search time if hider plays p0 and searcher optimally counters, suboptimality of p0, # search seqs required for convergence, optimal hiding strategy, p0, number of iterations

function run_algorithm10(qq::Array{Float64,1}, cc::Array{Float64,1}, n::Int64, accuracy::Int64, ϵ::Float64, δ::Float64, bounds::Array{Float64,1}, max_iters::Int64)
    
    p0 = cc./qq
    p0 = p0/sum(p0)

    #the n cycle matrix which initialises algorithm 10
    create = createLPmatp0(qq, cc, n, accuracy, false, true, [[1]]) #this is the n cycle matrix
    A = create[1]

    found=0 #indicator for convergence
    cnt = 0 #iteration counter
    L = 0 #initialise L as is redefined using old value of L (no need to initialise U)

    while found == 0 # the algorithm has not converged yet
        opts_constr=solveLP_constr(A, n, bounds) #the constrained solving for the algorithm
        h_opt = opts_constr[1] #pk: the optimal hiding strategy in the constrained problem

        #sometimes we are still getting 0 from solver (despite the p_i>=bound_i) due to bound_i being so small (10^(-22 sometimes)). 
        #below runs a check and replaces with the bound
        for i in 1:n
            if h_opt[i]<bounds[i]
                h_opt[i]=bounds[i]
            end
        end
        h_opt=h_opt/sum(h_opt) #renorm
        U = opts_constr[3] #vk, the value if the constrained problem which updates the upper bound in step 4
        #to update the lower bound:
        B = find_cond_ests(h_opt, qq, cc, n, accuracy) #use hider opt: creates a new row for the next iteration of the algorithm. Also used to calculate the lower bound
        L = max(L, sum(B.*h_opt)) # this is the EST of an optimal counter to h_opt: the new lower bound
        #now test if converged
        hide_less_bound = count(i->(i<δ),h_opt.-bounds) 

        if hide_less_bound == 0 #none of the constraints p_i≥δ_i are binding: conv crit ii) in Step 5
            if (U/L) - 1 < ϵ #conv crit i) in step 5
                found = 1
                p0_sub = 100*(U-create[3])/U
                if p0_sub < δ #just to tidy up the results, so don't have E-14 scientific notation. The version with the proposition is better for p0 subotimality
                    p0_sub = 0.
                end
                return(U, create[3], p0_sub, size(A,1),h_opt,p0, cnt+1)
            end
        end

        #function will have quit if algorithm converged, if not, add B to A for the next iteration
        A = vcat(A,B') #add the transpose of B to A for the next iteration
        cnt+=1
        
        ## to stop if the algorithm exceeds a certain number of iterations
        if cnt>max_iters
            print(h_opt); @printf "\n"; print(bounds); @printf "\n"; print(qq); @printf "\n"; print(cc); @printf "\n"
            @printf "Lower bound is %.9f and upper bound is %.9f\n" L U
            @printf "Have exceeded %.d iterations.\n" exceedence
            error("Exceeded above iterations of the algorithm, have stopped.")
        end
    end
end

