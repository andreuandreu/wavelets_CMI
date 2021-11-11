#using Pkg; Pkg.activate("./")
#Pkg.instantiate()
include("../src/read_bin.jl")
include("../src/basicent.jl")
using Plots
using Suppressor
using DrWatson
using DataFrames
using NPZ
using ConfParser


#Pkg.add("DelayEmbeddings")
using DelayEmbeddings


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
    println()
    println(dataX[1])
    println(typeof(dataX))
    return dataX, dataY

end






function inizilaize_one_embeding()
    "
    inizilaize delay embedings
    "

    Bin = 150
    kind_of_ergodicity = Entropies.RectangularBinning(Bin)
    root = "../data/exp_pro/correlations_TE"
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

function inizilaize_embedings(name_conf_file)

    ep = load_embedings_params(name_conf_file)

    kind_of_ergodicity = Entropies.RectangularBinning(ep.Bin)
    estimator = VisitationFrequency(kind_of_ergodicity)

    aux = LinRange(ep.jumpτ, ep.maxτ, 10)
    τ_range = round.(Int64, aux)
    τ_delays = (0, 0,  ep.maxτ ) #RIGhT ORDER FOR THE 3 dimensional embbeding case
    

    auxs = '_'
    for (i, e) in enumerate(τ_delays)
        aux = string(τ_delays[length(τ_delays)-i+1], auxs, ep.emb_dim[i])
    end
    
    aux_name_tag = ep.name_tag
    aux_bin_str = ep.Bin
    aux_root = ep.root
    name = "$aux_name_tag$aux_bin_str-$aux.txt"
    name_file = "$aux_root/$name"
    
    println()
    println(name)
    println()
    #joint= DelayEmbeddings.genembed(Dataset(x, y),  τ_delays, emb_dim )
    #joint= DelayEmbeddings.genembed(Dataset(x, y),  (0,0,50, 5, 10), (1,2,2,2,2) )
    return estimator, τ_range, τ_delays, name_file, ep.emb_dim
end


function load_embedings_params(name_conf_file)
    "
    load the parameters to define the embeding characteristics
    "

    conf = ConfParse(name_conf_file)
    parse_conf!(conf)

    # get and store config parameter
    #root = "/Users/andreu/Desktop/Dropbox/transfer_inormation_prague/code/correlations_TE"
    #root = "../data/exp_pro/correlations_TE"
    #maxτ = 50
    #jumpτ = 5
    #emb_dim = (2, 1, 1)
    #emb_dim = (1, 2, 2)
    #name_tag = "couplings_meanTE_Knn"

    root = retrieve(conf, "names", "root_folder")
    name_tag = retrieve(conf, "names", "name_tag")

    Bin  = parse(Int64,retrieve(conf, "emb_par", "bins")) 
    maxτ = parse(Int64,retrieve(conf, "emb_par", "max_tau"))
    jumpτ = parse(Int64,retrieve(conf, "emb_par", "jump_tau"))
    emb_dim_str = retrieve(conf, "emb_par", "embeding_dimension")

    emb_dim = []#(for e in collect(emb_dim_str) end)
    for e in collect(emb_dim_str)
        push!(emb_dim, parse(Int64,e))
    end

    emb_par = embedings_params(root, name_tag, Bin, maxτ, jumpτ, emb_dim)
    return emb_par

end


struct embedings_params
    "
    define the parameters to define the embeding characteristics
    "

    root :: String
    name_tag :: String
    Bin :: Int64
    maxτ :: Int64
    jumpτ :: Int64
    emb_dim :: Vector
 
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

function read_rossler()
    n_points = 131072
    couplings = LinRange(0, 0.25, 99)
    #root = "../data/exp_raw/binfiles/Rossler_bin_"
    #root = "/Users/andreu/Desktop/Dropbox/transfer_inormation_prague/code/binfiles/Rossler_bin_"
    root = "../../../TE_EG_project/package_CMI_prague/data/exp_raw/binfiles/Rossler_bin_"
    test_couplings = read_bin_couplings(root, n_points, couplings);
    return test_couplings
end


function TE_means(dataX, dataY,  output_name, τ_range, τ_delays, emb_dim, estimator)

    "
    compute many TE and measure a mean
    "

    open(output_name, "w") do f
        @suppress_err begin
            for i in 1:size(dataX)-1
                mean_TE = 0
                for t in τ_range
                    print(t)
                    ts =  changevector!(τ_delays, t, 2)
                    joint = DelayEmbeddings.genembed(Dataset(dataX, dataY),  ts, emb_dim )
                    e = tranfserentropy(joint, estimator)
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


input_file = ARGS[1]
nameX  = ARGS[2]
nameY = ARGS[3]
println(input_file)
println(nameY)
println(nameX)
estimator, τ_range, τ_delays, name_output_file, emb_dim = inizilaize_embedings(input_file)
#ENSO_manuel_wavelet_vecors_pywt_gaussian_amp.npy
#aux = nameX[:end-7]
#tail = "_pha.npy"
dataX, dataY = load_npy_data(nameX, nameY)
TE_means(dataX[1], dataY[1],  name_output_file, τ_range, τ_delays, emb_dim, estimator)
#output_file = ARGS[2]


