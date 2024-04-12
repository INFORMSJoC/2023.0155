#reads in data to reproduce a numerical experiment with N games, each with n boxes from data/Box_n

#arguments N, n, lowq, upq and upt correspond to the experiment to be reproduced: are used to generate the filename of the old experiment data (search games) to be read in
#δ, ϵ, bound_factor and max_iters are parameters used by Algorithm 10 on the reproduction of the experiment
#δp0 is the parameter used by the Proposition 4 test on the reproduction of the experiment
#see main file for an full explanation of these input parameters

#outputs 5 measures: whether p0 is optimal, the suboptimality of p0, the number of search sequences used (not reported in the paper), the number of iterations and the runtime;
#both in a raw output file (with each measure for all N games), and in a summary statistics file which corresponds to the results in Table 2
#which quantiles you wish for the summary statistics is determined by user input quan_vec (use decimals, e.g. [0.5,0.9,0.95,0.99])

#also outputs (for all N games) p0, the optimal p found, and the difference between them. This output was used for the analysis in Section 4.3


function reproduce_experiment_from_stored_data(n::Int64,N::Int64, lowq::Float64, upq::Float64, upt::Float64, accuracy::Int64, δ::Float64, ϵ::Float64, δp0::Float64, bound_factor::Float64, quan_vec::Array{Float64,1}, max_iters::Int64)


     ##setup storage for raw data
     raw_output = zeros(N, 5) #used to record, for each problem: whether p0 optimal, sub opt of p0, the number of sseqs (and iterations) and the runtime
     raw_output_head = ["p0_opt", "Subopt_p0", "sseqs_used", "number_iters", "runtime"] #header for raw_output

     p0_vs_opt = zeros(N, 3n) #records p0, opt p and the difference between them
     p0_vs_opt_head = vcat(["p0_$(i)" for i=1:n], ["p*_$(i)" for i=1:n], ["p*-p0_$(i)" for i=1:n]) #header


     ##read in data: use inputs to locate correct file
     #reads csv (if is there) and converts to an array
     read_file = "data/Box_$(n)/sprobs_$(n)_$(N)_$(lowq)_$(upq)_$(upt)_from_paper.csv" 
     if isfile(read_file)
        sprobs = CSV.File("data/Box_$(n)/sprobs_$(n)_$(N)_$(lowq)_$(upq)_$(upt)_from_paper.csv") |> Tables.matrix
     else
        print(read_file); print("\n")
        error("No file with the above path exists. Has a file been moved? Or was this experiment not run in the paper?")
     end
    
     
     #because the first time you run a function in julia takes longer, for a reliable runtime estimate, we run KL_algo_irrat_no_para_2024_no_prop once OUTSIDE the loop for a random problem
     #the results of this problem are not written anywhere
     qq = rand(Uniform(0.1,0.9), n); cc = rand(Uniform(1.,5.), n)
     bounds=bound_factor.*calc_lower_bound(qq,cc,n)
     run_algorithm10(qq,cc,n,accuracy,ϵ,δ,bounds, max_iters)

    ##run algorithm 10 and proposition 4 for each search game in data
     for i in 1:N
         if i % 200 == 0
             @printf "Run %d " i
         end

         #extract detection probabilities and search times for search game i
         qq= sprobs[i,1:n]
         cc= sprobs[i,n+1:2n]

         #run Proposition 4 and store binary result
         raw_output[i,1] = run_proposition4(qq, cc, n, accuracy, δp0)

         #setup for algorithm 10
         bounds=bound_factor.*calc_lower_bound(qq,cc,n) #calculate lower bound for algorithm 10
         t= time() #take runtime
         k = run_algorithm10(qq,cc,n,accuracy,ϵ,δ,bounds, max_iters) #run algorithm 10

         raw_output[i,5] = time() - t #store runtime since time t
         raw_output[i,2] = k[3] #suboptimality of p0
         raw_output[i,3] = k[4] #sseqs used
         raw_output[i,4] = k[7] #number of iterations
         
         p0_vs_opt[i,1:n] = k[6] # this is p0
         p0_vs_opt[i,n+1:2n] = k[5] #this is hider opt
         p0_vs_opt[i,2n+1:3n] = k[5].-k[6] #this is difference hider opt and p0
         
     end
    
     ##write output
     acpri = Int64(round(-log(ϵ)/log(10),digits=0)) #the number of decimal places in ϵ, for writing output

     #check if folder exists for this box, create if not
     if !isdir("results/results_from_paper/Box_$(n)")
        mkdir("results/results_from_paper/Box_$(n)")
     end

     ##write raw output
     #writes the suboptimality of p0, sseqs used, # iterations, runtimes
     CSV.write("results/results_from_paper/Box_$(n)/raw_data_$(n)_$(N)_$(lowq)_$(upq)_$(upt)_$(acpri)_from_paper.csv", DataFrame(raw_output, :auto), header = raw_output_head) 
     #writes p0, opt p and the difference between them
     CSV.write("results/results_from_paper/Box_$(n)/p0_vs_opt_$(n)_$(N)_$(lowq)_$(upq)_$(upt)_$(acpri)_from_paper.csv", DataFrame(p0_vs_opt, :auto), header = p0_vs_opt_head)


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

     #fill the second column quantile rows as missing (as no quantiles calculated for p0 optimal) 
     summary[2:(1 + length(quan_vec)),2] = repeat([NaN],length(quan_vec))

     #writes the summary
     CSV.write("results/results_from_paper/Box_$(n)/summary_$(n)_$(N)_$(lowq)_$(upq)_$(upt)_$(acpri)_from_paper.csv", DataFrame(summary, :auto), header = summary_head) 
    
end