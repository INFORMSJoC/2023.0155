##calculates the lower bounds η_i for p*_i from Proposition 8 in the paper

function calc_lower_bound(qq::Array{Float64,1}, cc::Array{Float64,1}, n::Int64)
    k = cc./qq
    upper = sum(k) #is U from the paper
    m = floor.(upper ./cc)
    term = k.*((1 .-qq).^-m) #the individual terms of the sum in the statement
    lower_bound=zeros(n)
    for i in 1:n
        lower_bound[i] = k[i]/(k[i]+sum(term[1:end .!= i]))
    end
    return(lower_bound)
end


