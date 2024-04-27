##for a search game with parameters qq and cc, creates a matrix game where the searcher can choose from Gittins index sequences against p0 determined by a permuation of 1,...,n

##uses the techniques in Section 4.1 and Appendix A to make calculations

##several inputs control which Gittins index sequences against p0 the searcher can choose from

# if allperms is true then looks at all n! GIS determined by a permutation of 1,...,n. Used for Proposition 4 test 
# if genperms is true, includes GIS determined by permutations 1...n, 2...1, and so on (for the 'cycle n'). Used for initialising Algorithm 10
# if more than one of the two above are true then the code errors
# if both are false, then uses the input in 'perms' 

# output is a matrix game 

function createLPmatp0(qq::Array{Float64,1}, cc::Array{Float64,1}, n::Int64, accuracy::Int64, allperms::Bool, genperms::Bool, perms::Array{Array{Int64,1},1})

    # first, we calculate the EST contribution after the first n searches
    # for tail probs, this is (t2-0)(1-q)+(t3-t2)(1-q)^2.... (where t2 is the second time we search the box in question)
    # qt1 is already taken care of in the second part of this code (start is t1+(t2-t1)(1-q)=qt1+t2(1-q))
    e=10.0^(-accuracy) # accuracy to which calculate each EST
    p0 = cc./qq
    p0 = p0/sum(p0)
    pmiss =(1 .-qq) #this is the probability that not been found so far,
    #will change throughout the search, initialise at 1-q since we
    #start after the first n searches (the tie)
    inds = p0.*qq.*pmiss./cc #starting after the first n searches (the tie)
    t=sum(cc) #time counter - we start after the first n searches (the tie)
    mxret = floor(log(minimum(1 .-qq))/log(maximum(1 .-qq))) + 1 #the maximum number of searches of one box that can occur between consec searches of another
    mxrett = mxret*t #an upper bound on the return time (also includes the box in question, which it should not, but only an upper bound)
    condW = zeros(n) #Cond ESTs - lower bound
    rat_UL = ones(n) #the ratio of the extra search time following the lousy policy and condW. (equiv to U/L-1)
    lbt = zeros(n) # the last time we searched that box - can start at 0 as we start with t2(1-q)
    a=0 #area to search
    nits = 0
    while maximum(rat_UL) > e
        #test if any of the indices are the same
        ss = 1
        soinds=sort(inds)
        while ss < n &&  soinds[ss+1] != soinds[ss]
            ss += 1
        end
        if ss<n #if so, need big floats
            inds = BigFloat.(p0).*BigFloat.(qq).*BigFloat.(pmiss)./BigFloat.(cc)
        end
        a=argmax(inds)
        t += cc[a]
        condW[a] += pmiss[a]*(t-lbt[a]) #need to do this and previous line BEFORE update pmiss
        pmiss[a] *= (1-qq[a])
        rat_UL[a] = pmiss[a]*(qq[a]*(t+mxrett)+mxrett*(1-qq[a]))/(qq[a]*condW[a]) #but here we want AFTER update pmiss
        inds[a] *= (1-qq[a])
        lbt[a] = t
    end

    #second, we loop over perms for the first n searches and add them to condW
    if allperms + genperms > 1
        error("At most one of these inputs should be true.")
    end
    perm_vec = collect(1:n)
    if allperms
        perms = collect(permutations(perm_vec))
    end
    if genperms
        perms = [circshift(perm_vec, n-i+1) for i=1:n] # shifts 1...n to 2...n1
    end
    #condW = reshape(condW,1,n) #so repeat `fills by row`
    #conds = repeat(condW,length(perms)) #to store the conditional search times,
    #each row is a perm, each column is a box
    ## The above was to initialise conds with the condW vector,
    #but since add a row vector in the perms loop anyway, easier just to initialise as zeros
    conds = zeros(length(perms),n) # to store the conditional search times

    #now everything is set up, the loop
    #now we calculate qt1 for each perm
    k = 0 #used for which row of conds to populate
    for perm in perms
        # @printf "%d %d %d %d\n" perm[1] perm[2] perm[3] perm[4]
        k=k+1
        times=copy(cc)
        permute!(times,perm) #reorder times according to the permutation
        csum=cumsum(times) #consec search times
        permute!(csum,invperm(perm)) #back to numerical ordering (1,2,..,n)
        csum =  csum.*qq #multiplies by the dp
        conds[k,:] = condW + csum
    end
    vp0 = sum(conds[1,:] .*p0) #this is the expected time to detection, needed for the test. Nothing special about 1, could use any row
    return(conds,perms,vp0)
end
