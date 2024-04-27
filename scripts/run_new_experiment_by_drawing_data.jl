#draws data using (23) in paper to create a new experiment with N games, each with n boxes

#arguments N, n, lowq, upq and upt are parameters determining how (23) is used to create search game data on which to run the experiment
#δ, ϵ, bound_factor and max_iters are parameters used by Algorithm 10 on the search game data
#δp0 is the parameter used by the Proposition 4 test on the search game data
#see main file for an full explanation of these input parameters

#outputs 5 measures: whether p0 is optimal, the suboptimality of p0, the number of search sequences used (not reported in the paper), the number of iterations and the runtime;
#both in a raw output file (with each measure for all N games), and in a summary statistics file which corresponds to the results in Table 2
#which quantiles you wish for the summary statistics is determined by user input quan_vec (use decimals, e.g. [0.5,0.9,0.95,0.99])

#also outputs (for all N games) p0, the optimal p found, and the difference between them. This output was used for the analysis in Section 4.3

#input experiment_name will appear in all output names, so can name the output 

#note that the Proposition 4 test is only run if n<=6, since it involves the computation of n! search sequences, which becomes computationally infeasible for n>6.


function run_new_experiment_by_drawing_data(n::Int64,N::Int64, lowq::Float64, upq::Float64, upt::Float64, accuracy::Int64,
      δ::Float64, ϵ::Float64, δp0::Float64,  bound_factor::Float64, quan_vec::Array{Float64,1}, max_iters::Int64, experiment_name::String)


     if n>6
        @printf "Since n=%d, the optimality test in Proposition 4 will not be run, as it requires the computation of %d!=%d search sequences.\n" n n factorial(n)  
     end

     ##setup storage for raw data
     raw_output = zeros(N, 5) #used to record, for each problem: whether p0 optimal, sub opt of p0, the number of sseqs (and iterations) and the runtime
     raw_output_head = ["p0_opt", "Subopt_p0", "sseqs_used", "number_iters", "runtime"] #header for raw_output

     p0_vs_opt = zeros(N, 3n) #records p0, opt p and the difference between them
     p0_vs_opt_head = vcat(["p0_$(i)" for i=1:n], ["p*_$(i)" for i=1:n], ["p*-p0_$(i)" for i=1:n]) #header


     sprobs = zeros(N, 2n) #used to store search games generated
     sprobs_head = vcat(["q$(i)" for i=1:n], ["t$(i)" for i=1:n]) #header
     
     
     ##draw search games
     qq = rand(Uniform(lowq,upq), n,N)
     cc = rand(Uniform(1.,upt), n,N)
     

     #because the first time you run a function in julia takes longer, for a reliable runtime estimate, we run KL_algo_irrat_no_para_2024_no_prop once OUTSIDE the loop for a random problem
     #the results of this problem are not written anywhere
     qq_init = rand(Uniform(0.1,0.9), n); cc_init = rand(Uniform(1.,5.), n)
     bounds=bound_factor.*calc_lower_bound(qq_init,cc_init,n)
     run_algorithm10(qq_init,cc_init,n,accuracy,ϵ,δ,bounds, max_iters)

     ##run algorithm 10 and proposition 4 for each search game in data
     for i in 1:N
         if i % 200 == 0
             @printf "Run %d " i
         end

         #extract and store detection probabilities and search times for search game i
         qq_vec=qq[:,i]; sprobs[i,1:n] = qq_vec 
         cc_vec=cc[:,i]; sprobs[i,n+1:2n] = cc_vec

         #if n<=6, run Proposition 4 and store binary result
         if n>6
            raw_output[i,1] = NaN
         else
            raw_output[i,1] = run_proposition4(qq_vec, cc_vec, n, accuracy, δp0)  
         end

         #setup for algorithm 10
         bounds=bound_factor.*calc_lower_bound(qq_vec,cc_vec,n) #calculate lower bound for algorithm 10
         t= time() #take runtime
         k = run_algorithm10(qq_vec,cc_vec,n,accuracy,ϵ,δ,bounds, max_iters) #run algorithm 10

         raw_output[i,5] = time() - t #store runtime since time t
         raw_output[i,2] = k[3] #suboptimality of p0
         raw_output[i,3] = k[4] #sseqs used
         raw_output[i,4] = k[7] #number of iterations
         
         p0_vs_opt[i,1:n] = k[6] # this is p0
         p0_vs_opt[i,n+1:2n] = k[5] #this is hider opt
         p0_vs_opt[i,2n+1:3n] = k[5].-k[6] #this is difference hider opt and p0
         
     end

     
     ##write output
     acpri = Int64(round(-log(ϵ)/log(10),digits=0)) #the number of decimal places in ϵ, for writing algorithm output

     #check if folder exists for this box, create if not
     if !isdir("results/new_experiments/Box_$(n)")
        mkdir("results/new_experiments/Box_$(n)")
     end

     ##write search problems
     CSV.write("results/new_experiments/Box_$(n)/sprobs_$(n)_$(N)_$(lowq)_$(upq)_$(upt)_$(experiment_name).csv", DataFrame(sprobs, :auto), header = sprobs_head)

     ##write raw output
     #writes the suboptimality of p0, sseqs used, # iterations, runtimes
     CSV.write("results/new_experiments/Box_$(n)/raw_data_$(n)_$(N)_$(lowq)_$(upq)_$(upt)_$(acpri)_$(experiment_name).csv", DataFrame(raw_output, :auto), header = raw_output_head) 
     #writes p0, opt p and the difference between them
     CSV.write("results/new_experiments/Box_$(n)/p0_vs_opt_$(n)_$(N)_$(lowq)_$(upq)_$(upt)_$(acpri)_$(experiment_name).csv", DataFrame(p0_vs_opt, :auto), header = p0_vs_opt_head)


     ##write summary output
     
     summary = zeros(1 + length(quan_vec), 6) #summary table
     summary_head = vcat("Quantile (0 is mean)", raw_output_head)

     #first row is for the means
     summary[1,2:6] = mean(raw_output, dims = 1)
     
     
     #other rows for the quantiles. We don't calculate quantiles for p0 optimal, as raw data is binary.
     for j in 3:6
        summary[2:(1 + length(quan_vec)),j] = quantile(raw_output[:,j-1], quan_vec)
     end

     #label the quantiles using the first column
     summary[2:(1 + length(quan_vec)),1] = quan_vec

     #fill the second column quantile rows with NA (as no quantiles calculated for p0 optimal) 
     summary[2:(1 + length(quan_vec)),2] = repeat([NaN],length(quan_vec))

     #writes the summary
     CSV.write("results/new_experiments/Box_$(n)/summary_$(n)_$(N)_$(lowq)_$(upq)_$(upt)_$(acpri)_$(experiment_name).csv", DataFrame(summary, :auto), header = summary_head) 
    
end