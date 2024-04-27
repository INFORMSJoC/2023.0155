##for a search game with parameters qq and cc, finds conditional expected search times (u(1,ξ),..,u(n,ξ) in the paper) when the searcher plays a Gittins index sequence ξ against a hiding strategy pp

##uses the techniques in Section 4.1 and Appendix A to make calculations

##returns a vector u(1,ξ),..,u(n,ξ)


function find_cond_ests(pp::Array{Float64,1}, qq::Array{Float64,1}, cc::Array{Float64,1}, n::Int64, accuracy::Int64)
    inds = zeros(n)
    e=10.0^(-accuracy) # accuracy to which calculate each EST
    condW = zeros(n) # to store the conditional search times
    mxret = floor(log(minimum(1 .-qq))/log(maximum(1 .-qq))) + 1 #the maximum number of searches of one box that can occur between consec searches of another
    mxrett = mxret*sum(cc) #an upper bound on the return time (also includes the box in question, which it should not, but only an upper bound)
    pmiss = ones(n) #this is the probability that not been found so far
    inds = pp.*qq./cc
    t=0 #time counter
    rat_UL = ones(n) #the ratio of the extra search time following the lousy policy and condW. (equiv to U/L-1)
    a=0 #area to search
    lbt = zeros(n) # the last time we searched that box
    nits=0
    while maximum(rat_UL) > e #equiv to while U/L -1 > e
        #test if any of the indices are the same
        ss = 1
        soinds=sort(inds)
        while ss < n &&  soinds[ss+1] != soinds[ss]
            ss += 1
        end
        if ss<n #if so, need big floats
            inds = BigFloat.(pp).*BigFloat.(qq).*BigFloat.(pmiss)./BigFloat.(cc)
        end
        a=argmax(inds)
        t += cc[a]
        condW[a] += pmiss[a]*(t-lbt[a]) #need to do this and previous line BEFORE update pmiss
        pmiss[a] *= (1-qq[a])
        rat_UL[a] = pmiss[a]*(qq[a]*(t+mxrett)+mxrett*(1-qq[a]))/(qq[a]*condW[a]) #but here we want AFTER update pmiss
        inds[a] *= (1-qq[a])
        lbt[a] = t
        nits+=1
        if nits >100000
            print(pp); @printf "\n"; print(qq); @printf "\n"; print(cc); @printf "\n"; print(rat_UL); print("\n")
            error("It took too long to calculate the conditional expected search time. Check what is printed above. Does pp have a very small or zero element?")
        end
    end
    return(condW)
end
