#using Pkg#; Pkg.activate("./")
#Pkg.instantiate()
#Pkg.add("InformationMeasures")
include("../src/read_bin.jl")
include("../src/basicent.jl")
using Suppressor
using DrWatson
using DataFrames
using NPZ
using ConfParser
using InformationMeasures
using Statistics

#Pkg.add("DelayEmbeddings")

#julia --project=. ./scripts/compute_TE.jl ./confs/config_embeding_char.ini

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


"loading functions"

function read_rossler()
    n_points = 131072
    couplings = LinRange(0, 0.25, 99)
    #root = "../data/exp_raw/binfiles/Rossler_bin_"
    #root = "/Users/andreu/Desktop/Dropbox/transfer_inormation_prague/code/binfiles/Rossler_bin_"
    root = "../../../TE_EG_project/package_CMI_prague/data/exp_raw/binfiles/Rossler_bin_"
    test_couplings = read_bin_couplings(root, n_points, couplings)

    return test_couplings
end



function load_rossler_data(root)
    "
    load rossler data
    "

    n_points = 131072
    couplings = LinRange(0, 0.25, 99)
    #root = "../data/exp_raw/binfiles/Rossler_bin_"
    #root = "/Users/andreu/Desktop/Dropbox/transfer_inormation_prague/code/binfiles/Rossler_bin_"
    ##root = "../../../TE_EG_project/package_CMI_prague/data/exp_raw/binfiles/Rossler_bin_"
    test_couplings = read_bin_couplings(root, n_points, couplings)
    return test_couplings
end


function load_npy_data(ep)

    datas = Array{AbstractArray}(undef, length(ep.sufixes), 1)

    for (i, s) in enumerate(ep.sufixes)

        name = ep.data_folder * ep.data_name * s
        datas[i] = npzread(name)

    end

    return datas[1], datas[2], datas[3]

end



function load_surrogates(ep, freq)

    sets_of_surrogates = Array{AbstractArray}(undef, length(freq), 1)

    for (i, f) in enumerate(freq)
        name = "./data/" * "surrogates/" * "surr_" * ep.data_name * f * ".npy"
        println(f, name)
        sets_of_surrogates[i] = npzread(name)

    end

    return sets_of_surrogates

end


"inizailzation functions"


function inizilaize_one_embeding()
    "
    inizilaize delay embedings
    "

    Bin = 150
    kind_of_ergodicity = Entropies.RectangularBinning(Bin)
    root = "../data/output/corr_TE_"
    lag = 25
    #τ_delays = (0,  -10, -5, 0, lag)
    τ_delays = (0, 0, lag)
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
    dim = (est = strEst, Bin = Bin, lag = lag)
    name = savename(prefix, dim, suffix)
end

function inizilaize_embedings(ep)



    kind_of_ergodicity = Entropies.RectangularBinning(ep.Bin)

    if ep.prob_est == "knn"
        estimator = KozachenkoLeonenko(kind_of_ergodicity)
    end
    if ep.prob_est == "VisFreq"
        estimator = VisitationFrequency(kind_of_ergodicity)
    end
    if ep.prob_est == "Kraskov"
        estimator = Kraskov(k = kind_of_ergodicity)
    end

    aux = LinRange(ep.jumpτ, ep.maxτ, Int(ep.maxτ / ep.jumpτ))
    τ_range = round.(Int64, aux)
    if length(ep.emb_dim) == 3
        τ_delays = (0, 0, ep.maxτ) #RIGhT ORDER FOR THE 3 dimensional embbeding case
    end
    if length(ep.emb_dim) == 5
        τ_delays = (0, 0, 0, 0, ep.maxτ) #RIGhT ORDER FOR THE 5 dimensional embbeding case
    end

    auxs = '_'
    for (i, e) in enumerate(τ_delays)
        aux = string(τ_delays[length(τ_delays)-i+1], auxs, ep.emb_dim[i])
    end

    aux_name_tag = ep.name_tag
    aux_bin_str = ep.Bin
    base_name_file = ep.name_tag * "_bin-" * string(ep.Bin) * "_eDim-" * string(length(ep.emb_dim))
    #base_name_file = "$aux_name_tag$aux_bin_str-$aux" * "_eDim-" * length(ep.emb_dim)

    #joint= DelayEmbeddings.genembed(Dataset(x, y),  τ_delays, emb_dim )
    #joint= DelayEmbeddings.genembed(Dataset(x, y),  (0,0,50, 5, 10), (1,2,2,2,2) )
    return estimator, τ_range, τ_delays, ep.emb_dim, base_name_file
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


    name_tag = retrieve(conf, "names", "name_tag")
    data_name = retrieve(conf, "names", "data_name")
    sufixes_str = retrieve(conf, "names", "sufixes")

    root = retrieve(conf, "folders", "root_folder")
    data_folder = retrieve(conf, "folders", "data_folder")
    export_folder = retrieve(conf, "folders", "export_folder")

    Bin = parse(Int64, retrieve(conf, "emb_par", "bins"))
    maxτ = parse(Int64, retrieve(conf, "emb_par", "max_tau"))
    jumpτ = parse(Int64, retrieve(conf, "emb_par", "jump_tau"))
    emb_dim_str = retrieve(conf, "emb_par", "embeding_dimension")

    prob_est = retrieve(conf, "prob_est", "prob_kind")

    sufixes = []#Array{String}(undef, 3, 1)
    for s in collect(sufixes_str)
        push!(sufixes, s)
    end

    emb_dim = []#(for e in collect(emb_dim_str) end)
    for e in collect(emb_dim_str)
        push!(emb_dim, parse(Int64, e))
    end



    emb_par = embedings_params(data_name, name_tag, sufixes, root, data_folder,
        export_folder, Bin, maxτ, jumpτ, emb_dim, prob_est)

    return emb_par

end


struct embedings_params
    "
    define the parameters to define the embeding characteristics
    "

    data_name::String
    name_tag::String
    sufixes::Vector
    root::String
    data_folder::String
    export_folder::String
    Bin::Int64
    maxτ::Int64
    jumpτ::Int64
    emb_dim::Vector
    prob_est::String

end


"compute stuff"

function write_entropies_tau(name, entropies)


    open(name, "w") do file
        for e in entropies

            write(file, "$e\n")
        end
    end

end

function TE_manyY_means(dataX, dataY, output_name, τ_range, τ_delays, emb_dim, estimator)

    "
    compute TE and measure a mean for a set of T timeseries
    "

    open(output_name, "w") do file
        for i = 1:size(dataY, 1)
            #println("doing something? ", i)
            entropies = TE_each_delay(dataX, dataY, output_name, τ_range, τ_delays, emb_dim, estimator)
            file_name = output_name[1:end-4] * "_Yset-" * string(i) * ".txt"
            write_entropies_tau(file_name, entropies)
            #println("now doing", dataChar[i], ' ', mean_TE)
            write(file, "$mean_TE\n")
        end


    end


end


function TE_each_delay(dataX, dataY, output_name, τ_range, τ_delays, emb_dim, estimator)

    "
    compute many TE for each tau delay
    "

    file_name = output_name[1:end-4] * "_each-tau" * ".txt"

    entropies = zeros(0)
    open(file_name, "w") do file
        @suppress_err begin
            for t in τ_range[1:end-1]
                println("tttt  ", τ_delays, t)
                ts = changevector!(τ_delays, t, 2)
                joint = DelayEmbeddings.genembed(Dataset(dataX, dataY), ts, emb_dim)
                entropy = tranfserentropy(joint, estimator)
                #entropy = tranfserentropy(joint, VisitationFrequency(b); embdim = 5, α =1.0, base =2)
                #entropy = tranfserentropy(joint, KozachenkoLeonenko(1,8),  2)
                println("doing something? delay ", t, "  ", entropy)
                write(file, "$entropy\n")
                append!(entropies, entropy)

            end
        end
    end
    return mean(entropies)

end



function MI_each_delay(dataX, dataY, output_name, τ_range, τ_delays)

    "
    compute many MI for each tau delay
    "
    file_name = output_name[1:end-4] * "_MI_each-tau" * ".txt"
    #τ_delays = (0, 0)
    #emb_dim = (1, 2)

    open(file_name, "w") do file
        for t in τ_range - 1
            mi_12 = get_mutual_information(dataX[1:end-t], dataY[t:end])
            ts = changevector!(τ_delays, t, 2)
            joint = DelayEmbeddings.genembed(Dataset(dataX, dataY), ts, emb_dim)
            entropy = tranfserentropy(joint, estimator)
            println("doing something? delay ", t, "  ", mi_12, "  ", entropy)
            write(file, "$mi_12  $entropy\n")
        end
        #println("now doing", dataChar[i], ' ', mean_TE)


    end


end



function data_rows_TE(ep, base_name_output_file, τ_range, τ_delays, emb_dim, estimator)

    "
    compute TE for each row of data in one of the two data files compared to the other
    "
    dataX, dataY, dataChar = load_npy_data(ep)

    if isdir(ep.root * ep.export_folder)
        println("BE AWARE folder already exixts!!! \n")
    else
        mkdir(ep.root * ep.export_folder)
    end


    for (i, char) in enumerate(dataChar)
        name_output_file = ep.root * ep.export_folder * base_name_output_file * "$(@sprintf("%.2f", char))" * ".txt"
        TE_means(dataX[i, :], dataY, name_output_file, τ_range, τ_delays, emb_dim, estimator)
        if i % 20 == 3
            println("entropies have been stored in here ", name_output_file)
        end
    end


end

"compute stuff with surrogates"

function TE_surrogate_set(surrogate_set, dataY, output_name, τ_range, τ_delays, emb_dim, estimator)
    entropies = zeros(0)
    count = 0
    for s in surrogate_set
        count += 1
        ampl = abs.(s[1:end-1])
        TE = TE_each_delay(ampl, dataY, output_name,
            τ_range, τ_delays, emb_dim, estimator)

        append!(entropies, TE)

        println(TE)
    end
    return entropies

end





name_conf_file = ARGS[1]
ep = load_embedings_params(name_conf_file)
estimator, τ_range, τ_delays, emb_dim, base_name_output_file = inizilaize_embedings(ep)

#data_rows_TE(ep, base_name_output_file, τ_range, τ_delays, emb_dim, estimator)

dataX, dataY, dataChar = load_npy_data(ep)
frequencies = ["_f009", "_f010"]
surrogates = load_surrogates(ep, frequencies)

println("sizes of surrogates ", size(surrogates), size(surrogates[1]))
output_name = name_output_file = ep.root * ep.export_folder * base_name_output_file * ".txt"#"./data/output/delay_rossler_phase.txt"
output_surr_nam = "./data" * "/surrogates" * "/surr_" * base_name_output_file * ".txt"
TE_surrogate_set(surrogates[1, :], dataY[1, :], output_surr_nam, τ_range, τ_delays, emb_dim, estimator)
#TE_each_delay(dataX[1, :], dataY[1, :], output_name, τ_range, τ_delays, emb_dim, estimator)
#MI_each_delay(dataX[1, :], dataY[1, :], output_name, τ_range, (0, 1))



