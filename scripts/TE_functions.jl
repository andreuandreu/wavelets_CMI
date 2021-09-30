"""this is a stupid function that given a vector A, copies it 
to a vector B but changes  one unique value of that vector 
in the position I, leaving all ther rest the same"""

using Suppressor

function changevector!(A, τ, I)
    B = Vector{Int64}(undef, length(A))#zeros(length(A))
     # LinRange(0, 0, length(A))#A#Array{Tuple{Int64,Int64}}(undef, length(couplings))
    @inbounds for i in eachindex(A)
        if i == I
            #append!( B, τ )
            B[i] = τ
        else
            #append!( B, A[i] )
            B[i] = A[i]
        end
    end
    return B
end


""" measures the mean of the TE for many lags """

function meanLag(x, y, js, τ_range, τs, est )
    mTE = Float64[]
    for t in τ_range
        ts =  changevector!(τs, t, length(τs))
        joint = DelayEmbeddings.genembed(Dataset(x, y),  ts, js )
        e = tranfserentropy(joint, est)  
        ##print(ts, ' ', e, ' ')
        push!(mTE, e) 
    end
    return mTE
end


"""Overly complex function that computes mean tranfer entropy (TE) for a dataset consisting in N sets
of pairs of becots dependent on a correlation coefficient e. For each pair the TE is computed for τ lags
and the mean of these lags is stored. There are  N such pairs and only 1 in n are computed,
by default n = 10. The value of the mean TE is stored in a text file of a given name"""

function TE_for_couplings(test_couplings, τs, js, est, τ_range, root2, name, n = 10)
    begin #@suppress_err
    TE = Vector{Vector{Float64}}(undef, 0)
    open("$root2/$name", "w") do f
    for i in 1:size(test_couplings, 3)-1
        if i%n == 0
            a = test_couplings[:,:,i]
            x = a[:,1]#/maximum(abs.(a[:,1]))
            y = a[:,2]#/maximum(abs.(a[:,2]))
            mTE = meanLag(x, y, js, τ_range, τs, est )
            push!(TE, mTE)
            println(i, " , ", couplings[i], ' ', mean(mTE), "±", StatsBase.std(mTE))
            aux = couplings[i]
            mean_TE = mean(mTE)
            sd_TE = StatsBase.std(mTE)
            write(f, "$aux $mean_TE $sd_TE\n")
        end
    end
    end
    return TE
    end 
end