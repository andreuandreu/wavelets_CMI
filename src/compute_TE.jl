"""Overly complex function that computes mean tranfer entropy (TE) for a dataset consisting in N sets
of pairs of becots dependent on a correlation coefficient e. For each pair the TE is computed for τ lags
and the mean of these lags is stored. There are  N such pairs and only 1 in n are computed,
by default n = 10. The value of the mean TE is stored in a text file of a given name"""

function TE_for_couplings(test_couplings, couplings, τs, js, est, τ_range, root2, name, n = 10)
    TE = Array{Float64}(undef,  length(couplings), 1)
    mean_TE_vector = Array{Float64}(undef,  length(couplings), 1)
    @suppress_err begin
    open("$root2/$name", "w") do f
        for i in 1:size(test_couplings, 3)-1
            println(i)
            if i%n == 3
                mean_TE = 0
                a = test_couplings[:,:,i]
                x = a[:,1]#/maximum(abs.(a[:,1]))
                y = a[:,2]#/maximum(abs.(a[:,2]))
                for t in τ_range
                    print(t)
                    ts =  changevector!(τs, t, length(τs))
                    joint = DelayEmbeddings.genembed(Dataset(x, y),  τs, js )
                    #e = tranfserentropy(joint, est; embdim = 5, α =1.0, base =2)
                    e = tranfserentropy(joint, est,  2)
                    TE[i] = e
                    mean_TE += e/length(τ_range)
                end
                println(i, " , ", couplings[i], ' ', mean_TE)
                aux = couplings[i]
                write(f, "$aux $mean_TE\n")

                mean_TE_vector[i] = mean_TE
            end
        end
    end
    end
    return mean_TE_vector, TE
end

"""
Function to automatically create a string given variables of the code
"""

function name_of_file(est, Bin, maxτ, τs, js)
    strEst = string(est)[1:3]
    prefix = "couplings_TE"
    suffix = "txt"

    tsjs = '_'
    for (i, e) in enumerate(τs)
        tsjs = string(τs[length(τs)-i+1], tsjs, js[i])
    end

    #name = "$prefix-$strEst-$Bin-lag_$lag.$suffix"
    d = (est = strEst,  Bin = Bin, maxlag = maxτ, tsjs = tsjs  )
    name = savename(prefix, d , suffix)
    println(name)
    return name
end

"""this is a stupid function that given a vector A, copies it
to a vector B but changes  one unique value of that vector
in the position I, leaving all the rest the same"""

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