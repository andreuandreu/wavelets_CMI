include("../src/read_bin.jl")
include("../src/basicent.jl")
using Plots
using Suppressor
using DrWatson





function changevector!(A, τ, I)

    """this is a stupid function that given a vector A, copies it 
    to a vector B but changes  one unique value of that vector 
    in the position I, leaving all ther rest the same"""
    
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



"
load data
"

n_points = 131072
couplings = LinRange(0, 0.25, 99)
#root = "../data/exp_raw/binfiles/Rossler_bin_"
#root = "/Users/andreu/Desktop/Dropbox/transfer_inormation_prague/code/binfiles/Rossler_bin_"
root = "../../../TE_EG_project/package_CMI_prague/data/exp_raw/binfiles/Rossler_bin_"
test_couplings = read_bin_couplings(root, n_points, couplings)




"
inizilaize delay embedings
"

Bin = 150
kind_of_ergodicity = Entropies.RectangularBinning(Bin)
root2 = "../data/exp_pro/correlations_TE"
lag = 25
#τs = (0,  -10, -5, 0, lag)
τs = (0,  0, lag)
#js = (2, 1, 1, 1, 1)
#js = (1, 2, 2, 2, 2)
#js = (1, 2, 2)
js = (2, 1, 1)
aux = '_'
for (i, e) in enumerate(τs)
    aux = string(τs[i],aux,js[i])
end
#est = VisitationFrequency(b)
#est = KozachenkoLeonenko(RecBin,0)
est = Kraskov(k = Bin)
#name = "couplings_TE_cbv$RecBin-$aux.txt"
strEst = string(est)[1:3]
prefix = "couplings_TE"
suffix = "txt"
name = "$prefix-$strEst-$RecBin-lag_$lag.$suffix"
println(name )
d = (est = strEst,  Bin = Bin, lag = lag)
name = savename(prefix, d , suffix)

#joint= DelayEmbeddings.genembed(Dataset(x, y),  τs, js )
#joint= DelayEmbeddings.genembed(Dataset(x, y),  (0,0,50, 5, 10), (1,2,2,2,2) )


"
compute one TE
"

TE = Array{Float64}(undef,  length(couplings), 1)
@suppress_err begin
open("$root2/$name", "w") do f
    for i in 1:size(test_couplings, 3)-1
        if i%10 == 3
            a = test_couplings[:,:,i]
            x = a[:,1]
            y = a[:,2]
            joint = DelayEmbeddings.genembed(Dataset(x, y),  τs, js )

            e = tranfserentropy(joint, est)

            println(i, ' ', couplings[i], ' ', e)
            aux = couplings[i]
            TE[i] = e
            write(f, "$aux $e\n")
        end
    end
end
end


"
compute many TE and measure a mean
"

#root2 = "/Users/andreu/Desktop/Dropbox/transfer_inormation_prague/code/correlations_TE"
root2 = "../data/exp_pro/correlations_TE"

maxτ = 50
jumpτ = 5
aux = LinRange(jumpτ, maxτ, 10)
τ_range = round.(Int64, aux)
println(τ_range)
#τs = (0,  10, 5, 0, maxτ)
#js = (2, 1, 1, 1, 1)
#js = (1, 2, 2, 2, 2)

τs = (0, 0,  maxτ ) #RIGhT ORDER FOR THE 3 dimensional embbeding case
#js = (2, 1, 1)
js = (1, 2, 2)

aux = '_'
for (i, e) in enumerate(τs)
    aux = string(τs[length(τs)-i+1],aux,js[i])
end
name = "couplings_meanTE_Knn$RecBin-$aux.txt"

open("$root2/$name", "w") do f
    @suppress_err begin
        for i in 1:size(test_couplings, 3)-1
            mean_TE = 0
            a = test_couplings[:,:,i]
            x = a[:,1]
            y = a[:,2]
            for t in τ_range
                print(t)
                ts =  changevector!(τs, t, 2)
                joint = DelayEmbeddings.genembed(Dataset(x, y),  ts, js )
                e = tranfserentropy(joint, VisitationFrequency(kind_of_ergodicity))
                #e = tranfserentropy(joint, VisitationFrequency(b); embdim = 5, α =1.0, base =2)
                #e = tranfserentropy(joint, KozachenkoLeonenko(1,8),  2)
                mean_TE += e/length(τ_range)  
            end
            aux = couplings[i]
            println(couplings[i], ' ', mean_TE)
            write(f, "$aux $mean_TE\n")
        end
    end
end


