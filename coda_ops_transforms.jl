#testing for julia functions :)

using Statistics
using LinearAlgebra
#using InvertedIndices
#operations:
#const LiefiesMatrix = Union{Matrix{Missing}, Matrix{<:Real}}
LiefiesMatrix{T <: Union{Missing, Real}} = AbstractMatrix{T}
#function f=some_function(x::LiefiesMatrix) :D


# function closure(data::Union{Array{Union{Float64,Missing},2},Array{Union{Int64, Missing},2},Vector{Float64},Vector{Int64}}) #Union{Vector{Regex}, NTuple{N, Regex}}
#::LiefiesMatrix) #Union{Vector{Regex}, NTuple{N, Regex}}
function closure(data::T) where T
    cvec = Array{Union{Missing, Float64}}(undef, size(data)...)
    for i in 1:size(data,1)
      cvec[i,:] = data[i,:] ./ sum(skipmissing(data[i,:])) #work around for sum not accepting skipmissing with dims argument
    end
    return cvec
end

function center(data::T) where T #inherit type T from data :D
    #data must be compositions!
    gvec = Matrix{Union{Missing, Float64}}(undef, 1, size(data,2))
    data = closure(data) #just in case
    data = BigFloat.(data)
    for j in 1:size(gvec,2) #by columns
      gvec[j] = prod(skipmissing((data[:,j]))) ^ (1/size(collect(skipmissing(data[:,j])),1)) #work around for sum not accepting skipmissing with dims argument
    end
    gvec=closure(gvec)
    return gvec
end


function ait_var(data::T) where T
    cen = center(data)
    avar = zeros(Float64, size(data,2), size(data,2))
    #avar = Matrix{Float64}(undef, size(data,2), size(data,2))
    #avar = Matrix{Union{Missing, Float64}}(undef, size(data,2), size(data,2))
    for i=1:size(data,2)
        for j=i:size(data,2)
            n=min(length(collect(skipmissing(data[:,i]))),length(collect(skipmissing(data[:,j]))))
            avar[i,j]=(1/(n-1))*sum(skipmissing((log.(data[:,i]./data[:,j]) .- log(cen[i]/cen[j])).^2))
        end
    end
    #avar = Symmetric(avar) #to make symmetric :)
    return avar
end

function geo(x)
    geo_mean = prod(x) .^ (1 ./ length(x))
    return geo_mean
end

function ait_inner(x,y)
    ai=sum( log.(x ./ geo(x)) .* log.(y ./ geo(y))   )
    return ai

    # ai=0.0        #alternative calculation:
    # for i=1:length(x)
    #     ai = ai + sum((log.(x[i] ./ x) .* (log.(y[i] ./ y))))
    # end
    # ai = 1/(2*length(x)) * ai
    #return ai
end

function ait_norm(x)
    an=0.0
    for i=1:length(x)
        an= an + sum((log.(x[i] ./ x) .^2))
    end
    an = sqrt(1/(2*length(x)) * an)
    return an
end

function square_ait_norm(x)
    an=0.0
    for i=1:length(x)
        an= an + sum((log.(x[i] ./ x) .^2))
    end
    an = (1/length(x)) * 1/(2*length(x)) * an
    return an
end

function ait_dist(x,y)
    if length(x) == length(y)
        ad = 0.0
        for i=1:length(x)
            ad = ad + sum((log.(x[i] ./ x) .- log.(y[i] ./ y)).^2)
        end
        ad = sqrt(1/(2*length(x)) * ad)
        return ad
    else
        println("x and y must be the same length")
    end
end

function perturb(x,y)
    p=closure(x) .* closure(y) #works even if x is mulitple rows!
    p=closure(p)
    return p
end

function perturb_diff(x,y)
    p=closure(x) .* inverse(y) #works even if x is mulitple rows!
    p=closure(p)
    return p
end

function power(x,a)
    p=closure(x) .^ a
    p=closure(p)
    return p
end

function inverse(x)
    ix = closure(1 ./ x)
    return ix
end

function clr(data)
    cvec = Matrix{Union{Missing, Float64}}(undef, size(data))
    for i=1:size(data,1) #by rows
        gmean = prod(skipmissing(data[i,:])) ^ (1/size(collect(skipmissing(data[i,:])),1)) #work around for sum not accepting skipmissing with dims argument
        cvec[i,:]=log.(data[i,:]./gmean)
    end
    return cvec
end

function ilr(x,sbp::Matrix{<:Real})

    bal=zeros(size(sbp)) #create orthonormal basis based on partition
    for i=1:size(sbp,1) #number of rows of sbp
        rr=findall(sbp[i,:] .== 1)
        ss=findall(sbp[i,:] .== -1)
        r=length(rr)
        s=length(ss)
        bal[i,rr] .= (1/r).*sqrt((r*s)/(r+s))
        bal[i,ss] .= -(1/s).*sqrt((r*s)/(r+s))
    end

    ilr_x=log.(x)*bal'
    return ilr_x, bal
end

function invilr(ilr_x,bal)
    #balance refers to contrast matrix (it is not the sbp!)
    x=invclr(ilr_x*bal);
    return x
end

function invclr(x_clr)
    x=closure(exp.(x_clr));
    return x
end

function amalga(x,col1,col2)

    A = zeros(size(x,2),3)
    A[col1,1] = 1
    A[col2,2] = 1
    A[setdiff(collect(1:size(x,2)),[col1;col2]),3] .= 1

    y = x * A
    return y, A

end


function bayes_multi_zero_replace(c,method,prior)
    #c is a count vector!
    #tj is the prior estimate
    #s is the strength parameter

    samples, D = size(c) #number of rows then columns
    n=sum(c,dims=2)

    if prior == "uniform"
        t=1/D .* ones(size(c))
    elseif prior == "modified"
        println("Assuming prior based on ratios found within dataset")
        t=zeros(size(c))
        for q=1:samples
            global t
            to_use = setdiff(1:samples, q)
            t[q,:]=sum(c[to_use,:],dims=1)
        end
        t=t./sum(t,dims=2)
    end

    #common imprecise Dirichlet models considered here: Bayes-Laplace and Square-root
    if method == "Bayes-Laplace"
        s = D * ones(samples,1)
    elseif method == "Square-root"
        s = sqrt.(sum(c,dims=2)) #total trials per sample!
    elseif method == "Geometric"
        #if these are very small numbers, can result in a zero multiplicative product
        println("accounting for very small values here...")
        temp = 10 .^ ((1/D) .* sum(log10.(t),dims=2))
        s = 1 ./ temp    
        #s = 1 ./ prod(closure(t),dims=2).^(1/D)
    elseif method == "Jeffreys"
        s = D/2 * ones(samples,1)
    elseif method == "Perks"
        s = ones(samples,1)
    else
        println("Not a recognized method; please choose from Bayes-Laplace, Square-root, Geomteric, Jeffreys, or Perks")
    end

    println(method)
    #@show s
    #multiplicative part: only adjust rows for which there is a zero
    #fieldnames()...
    x=closure(c)
    r=deepcopy(x)
    global r, x, t, s
    for i=1:samples
        if any(x[i,:] .== 0) #row contains zeros
            zeroind=findall(x[i,:] .== 0) #work on the columns
            nonzero=findall(x[i,:] .!= 0)
            r[i,zeroind] = t[i,zeroind] .* s[i]./(n[i] .+ s[i])
            r[i,nonzero] = x[i,nonzero] .* (1 - sum(r[i,zeroind]))
        end
    end

    return r
    # x=closure(c)
    # r=x
    # zeroind=findall(x .== 0)
    # nonzero=findall(x .!= 0)
    # colind=map(x->x.I[2], zeroind)
    # rowind=map(x->x.I[1], zeroind)
    #
    # r[zeroind] = t[zeroind] .* s[colind]./(n[colind] .+ s[colind])
    # r[nonzero] = x[nonzero] .* (1 - sum(r[zeroind]))

end

#To compare with R package- that one seems to be doing strange thins with BDL values
# g=Matrix{Union{Missing, Float64}}(undef, 1, size(y,2))
# for j in 1:size(y,2) #by columns
#   g[j] = prod(skipmissing(y[:,j])) ^ (1/size(collect(skipmissing(y[:,j])),1)) #work around for sum not accepting skipmissing with dims argument
# end
#
# #(0.3 * x).^(1/2) = nrm
# g0=sum(g[1,[1,2,3]])
# #sum(g[1,[1,2,4]]) + q = nrm
# #nrm = y[1,4] + g0
#
# r_cen=[0.5116703 0.2558351 0.1193897 0.1131049]
# mm=((g[1]./r_cen[1] - g0).^2)./y[1,4]
