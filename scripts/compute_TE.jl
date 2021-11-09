include("../src/read_bin.jl")
include("../src/basicent.jl")
using Plots
using Suppressor
using DrWatson
using DataFrames
using NPZ






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



function load_data(root)
    "
    load data
    "

    n_points = 131072
    couplings = LinRange(0, 0.25, 99)
    #root = "../data/exp_raw/binfiles/Rossler_bin_"
    #root = "/Users/andreu/Desktop/Dropbox/transfer_inormation_prague/code/binfiles/Rossler_bin_"
    ##root = "../../../TE_EG_project/package_CMI_prague/data/exp_raw/binfiles/Rossler_bin_"
    test_couplings = read_bin_couplings(root, n_points, couplings)
    return test_couplings
end


function load_npy_data(nameX, nameY)

    dataX = npzread(nameX)
    dataY = npzread(nameY)

    return dataX, dataY

end




input_file = ARGS[1]
output_file = ARGS[2]

data_frame = readtable(input_file, separator = ',')
writetable(output_file, data_frame)

function extract_instructions(input_file)

    open(input_file, "r") do in_file
        # Use a for loop to process the rows in the input file one-by-one
        for line in eachline(in_file)

        end
    end
    
end

function inizilaize_one_embeding()
    "
    inizilaize delay embedings
    "

    Bin = 150
    kind_of_ergodicity = Entropies.RectangularBinning(Bin)
    root2 = "../data/exp_pro/correlations_TE"
    lag = 25
    #τ_delays = (0,  -10, -5, 0, lag)
    τ_delays = (0,  0, lag)
    #emb_dim = (2, 1, 1, 1, 1)
    #emb_dim = (1, 2, 2, 2, 2)
    #emb_dim = (1, 2, 2)
    emb_dim = (2, 1, 1)
    aux = '_'
    for (i, e) in enumerate(τ_delays)
        aux = string(τ_delays[i], aux, emb_dim[i])
    end
    #estimator = VisitationFrequency(b)
    #estimator = KozachenkoLeonenko(RecBin,0)
    estimator = Kraskov(k = Bin)
    #name = "couplings_TE_cbv$RecBin-$aux.txt"
    strEst = string(estimator)[1:3]
    prefix = "couplings_TE"
    suffix = "txt"
    name = "$prefix-$strEst-$RecBin-lag_$lag.$suffix"
    println(name)
    dim = (est = strEst,  Bin = Bin, lag = lag)
    name = savename(prefix, dim , suffix)
end

function inizilaize_embedings()
    #root2 = "/Users/andreu/Desktop/Dropbox/transfer_inormation_prague/code/correlations_TE"
    root2 = "../data/exp_pro/correlations_TE"

    Bin = 150
    kind_of_ergodicity = Entropies.RectangularBinning(Bin)
    
    maxτ = 50
    jumpτ = 5
    aux = LinRange(jumpτ, maxτ, 10)
    τ_range = round.(Int64, aux)
    τ_delays = (0, 0,  maxτ ) #RIGhT ORDER FOR THE 3 dimensional embbeding case
    #emb_dim = (2, 1, 1)
    emb_dim = (1, 2, 2)

    aux = '_'
    for (i, e) in enumerate(τ_delays)
        aux = string(τ_delays[length(τ_delays)-i+1], aux, emb_dim[i])
    end
    name = "couplings_meanTE_Knn$RecBin-$aux.txt"
    name_file = "$root2/$name"
    #joint= DelayEmbeddings.genembed(Dataset(x, y),  τ_delays, emb_dim )
    #joint= DelayEmbeddings.genembed(Dataset(x, y),  (0,0,50, 5, 10), (1,2,2,2,2) )
end

function compute_one_TE(name_file, test_couplings, couplings, estimator, τ_delays, emb_dim)
    "
    compute one TE
    "
    TE = Array{Float64}(undef,  length(couplings), 1)
    
    @suppress_err begin
    open("$name_file", "w") do file
        for i in 1:size(test_couplings, 3)-1
            if i%10 == 3
                aux = test_couplings[:,:,i]
                dataX = aux[:,1]
                dataY = aux[:,2]
                joint = DelayEmbeddings.genembed(Dataset(dataX, dataY),  τ_delays, emb_dim )

                entropy = tranfserentropy( joint, estimator)

                println(i, ' ', couplings[i], ' ', entropy)
                aux2 = couplings[i]
                TE[i] = entropy
                write(file, "$aux2 $entropy\n")
            end
        end
    end
    end

end 


function TE_means( output_name, τ_delays, emb_dim, kind_of_ergodicity)

    "
    compute many TE and measure a mean
    "

    open("$root2/$output_name", "w") do f
        @suppress_err begin
            for i in 1:size(test_couplings, 3)-1
                mean_TE = 0
                aux = test_couplings[:,:,i]
                dataX = aux[:,1]
                dataY = aux[:,2]
                for t in τ_range
                    print(t)
                    ts =  changevector!(τ_delays, t, 2)
                    joint = DelayEmbeddings.genembed(Dataset(dataX, dataY),  ts, emb_dim )
                    e = tranfserentropy(joint, VisitationFrequency(kind_of_ergodicity))
                    #e = tranfserentropy(joint, VisitationFrequency(b); embdim = 5, α =1.0, base =2)
                    #e = tranfserentropy(joint, KozachenkoLeonenko(1,8),  2)
                    mean_TE += e/length(τ_range)  
                end
                aux2 = couplings[i]
                println(couplings[i], ' ', mean_TE)
                write(file, "$aux2 $mean_TE\n")
            end
        end
    end

end

