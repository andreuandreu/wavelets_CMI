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

        name = ep.data_folder * ep.input_folder * ep.in_data_tag * s
        datas[i] = npzread(name)

    end

    return datas[1], datas[2], datas[3]

end



function load_surrogates(freq, surr_root)

    sets_of_surrogates = Array{AbstractArray}(undef, length(freq), 1)

    for (i, f) in enumerate(freq)

        name = surr_root * f * ".npy"
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
    lag = 25
    #τ_delays = (Int,  -10, -5, 0, lag)
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

function inizilaize_embedings(ep, serieChar)

    backwardLag = trunc(Int, serieChar / 4.0)

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
        τ_delays = (backwardLag * 3, backwardLag * 2, backwardLag, 0, ep.maxτ) #RIGhT ORDER FOR THE 5 dimensional embbeding case
    end


    #joint= DelayEmbeddings.genembed(Dataset(x, y),  τ_delays, emb_dim )
    #joint= DelayEmbeddings.genembed(Dataset(x, y),  (0,0,50, 5, 10), (1,2,2,2,2) )
    return estimator, τ_range, τ_delays, ep.emb_dim
end

function get_emb_dim_name(ep)

    auxs = ""
    for (i, e) in enumerate(ep.emb_dim)
        auxs = auxs * string(ep.emb_dim[i])
    end

    extra_name_file = "_eDim-" * auxs#string(length(ep.emb_dim))
    return extra_name_file

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



    data_folder = retrieve(conf, "folders", "data_folder")
    input_folder = retrieve(conf, "folders", "input_folder")
    export_folder = retrieve(conf, "folders", "export_folder")

    in_data_tag = retrieve(conf, "names", "in_data_tag")
    sufixes_str = retrieve(conf, "names", "sufixes")
    pha_amp_com_str = retrieve(conf, "names", "pha_amp_com")


    Bin = parse(Int64, retrieve(conf, "emb_par", "bins"))
    maxτ = parse(Int64, retrieve(conf, "emb_par", "max_tau"))
    jumpτ = parse(Int64, retrieve(conf, "emb_par", "jump_tau"))
    emb_dim_str = retrieve(conf, "emb_par", "embeding_dimension")
    period_range_str = retrieve(conf, "emb_par", "period_range")

    name_tag = retrieve(conf, "prob_est", "name_tag")
    prob_est = retrieve(conf, "prob_est", "prob_kind")

    sufixes = []#Array{String}(undef, 3, 1)
    for s in collect(sufixes_str)
        push!(sufixes, s)
    end
    pha_amp_com = []
    for e in collect(pha_amp_com_str)
        push!(pha_amp_com, e)
    end

    emb_dim = []#(for e in collect(emb_dim_str) end)
    for e in collect(emb_dim_str)
        push!(emb_dim, parse(Int64, e))
    end


    number_of_freq = 11
    max_period = parse(Int64, period_range_str[2])
    min_period = parse(Int64, period_range_str[1])
    period_ratio = log10(max_period)
    phase_periods = []
    for i in range(0, number_of_freq, step = 1)
        append!(phase_periods, min_period + i^period_ratio)
    end


    emb_par = embedings_params(data_folder, input_folder, export_folder, in_data_tag,
        sufixes, pha_amp_com, Bin, maxτ, jumpτ, emb_dim, phase_periods, name_tag, prob_est)

    return emb_par

end


struct embedings_params
    "
    define the parameters to define the embeding characteristics
    "
    data_folder::String
    input_folder::String
    export_folder::String
    in_data_tag::String
    sufixes::Vector
    pha_amp_com::Vector
    Bin::Int64
    maxτ::Int64
    jumpτ::Int64
    emb_dim::Vector
    period_range::Vector
    name_tag::String
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


function TE_each_delay(dataX, dataY, ep, serieChar)

    "
    compute many TE for each tau delay
    "

    estimator, τ_range, τ_delays, emb_dim = inizilaize_embedings(ep, serieChar)

    entropies = zeros(0)
    @suppress_err begin
        for t in τ_range[1:end-1]
            ts = changevector!(τ_delays, t, 2)
            joint = DelayEmbeddings.genembed(Dataset(dataX, dataY), ts, emb_dim)
            entropy = tranfserentropy(joint, estimator)
            #println("doing something? delay ", t, "  ", entropy)
            append!(entropies, entropy)

        end
    end
    return mean(entropies)
end



function TE_each_delay(dataX, dataY, output_name, ep, serieChar)

    "
    compute many TE for each tau delay and write it in a file
    "

    estimator, τ_range, τ_delays, emb_dim, extra_name = inizilaize_embedings(ep, serieChar)
    file_name = output_name[1:end-4] * extra_name * "_each-tau" * ".txt"

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
                #println("doing something? delay ", t, "  ", entropy)
                write(file, "$entropy\n")
                append!(entropies, entropy)

            end
        end
    end
    return mean.(entropies)

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
    #mean.()#eachcol(df)

end


function TE_each_column(dataX, dataY, output_name, ep, serieChar)

    "
    compute TE and measure a mean for a set of N timeseries (the Y one)
    "

    open(output_name, "w") do file
        for i = 1:size(dataY, 1)
            mean_entropy = TE_each_delay(dataX, dataY[i, :], ep, serieChar)# output_name, 
            #file_name = output_name  * ".txt"
            #write_entropies_tau(file_name, entropies)
            #println("now doing", dataChar[i], ' ', mean_TE)
            write(file, "$serieChar $mean_entropy\n")
        end
    end
end


function TE_each_row(dataX, dataY, dataChar, ep, base_name_output_file)

    "
    compute TE for each row of data (the x one) in one of the two data files compared to the other
    "

    for (i, char) in enumerate(dataChar)
        name_output_file = ep.data_folder * ep.export_folder *
                           base_name_output_file * "_row-p$(@sprintf("%.2i", char))" * ".txt"

        TE_each_column(dataX[i, :], dataY, name_output_file, ep, char)
        if i % 20 == 3
            println("entropies have been stored in here ", name_output_file)
        end
    end


end





"compute stuff with surrogates"

function TE_surrogate_set(surrogate_set, output_name, dataX, ep, serieChar)

    "compute the TE of a guiven time series with a guiven set of surrogates, 
    measure the mean and store it given a pattern of x y names or characteristics"


    entropies = zeros(0)
    for i = 1:size(surrogate_set)[1]
        s = surrogate_set[i, :]
        if occursin("-Amp", output_name)
            compute = abs.(s[1:end-1])
        elseif occursin("-Pha", output_name)
            compute = angle.(s[1:end-1])
        else
            println("BE AWARE you are missing '_amp' or '_pha' in the name \n")
        end

        TE = TE_each_delay(dataX, compute, ep, serieChar)
        #if i % 44 == 1
        #    println("random eeeeeent  ", i, "  ", TE)
        #end
        append!(entropies, TE)
    end

    return mean(entropies)

end

function TE_all_surrogates(names_surrogates, output_name, dataSerie, serieChar, ep)
    #τ_range, τ_delays, emb_dim, estimator
    "
    for each set of surrogates, compute the TE with a given data series and print it in a file
    "
    extra_name_emb = get_emb_dim_name(ep)
    name_out_file = output_name * extra_name_emb * "_p$(@sprintf("%.2i", serieChar))" * ".txt"


    open(name_out_file, "w") do file
        for n in names_surrogates
            surrogate_set = npzread(ep.data_folder * ep.input_folder * n)
            one_mean = TE_surrogate_set(surrogate_set, output_name, dataSerie, ep, serieChar)

            println("mmmmmmmm   ", n, "   ", one_mean)
            x = n[end-7:end-4]
            write(file, "$x $serieChar $one_mean\n")
        end
    end
end


function TE_data_rows_surrogates(base_name_output_file, dataChar, dataSeries, 
        tagSeries, amp_or_phase_surr,  ep)

    "
    compute TE for each row of data for each set of surrogates with varing frequencies (periods)
    "

    surr_root = "surr_circ_"
    
    #amp_or_phase_surr = "-pha"

    println("NIGHTMARE in ", ep.data_folder * ep.input_folder * surr_root * ep.in_data_tag)
    names_surrogates = frequencies_names(ep.data_folder * ep.input_folder, surr_root * ep.in_data_tag)

    output_surr_root = ep.data_folder * ep.export_folder * surr_root * base_name_output_file
    println("NIGHTMARE out ", output_surr_root)

    for (i, char) in enumerate(dataChar)

        name_output_file = output_surr_root * tagSeries *"_Su"*amp_or_phase_surr
        println("NANANANA ", name_output_file)

        TE_all_surrogates(names_surrogates, name_output_file, dataSeries[i, :], char, ep)

    end


end




function frequencies_names(path, key)

    "
    given a path and a key name, returns an array with all the files on the path with the correspoonding key
    "

    names = filter(x -> occursin(key, x), readdir(path))
    #println("yeah guy", names_surrogates)
    if length(names) < 1
        throw(DomainError("this character is bad"))
    end
    return names

end


@time begin

    name_conf_file = ARGS[1]

    ep = load_embedings_params(name_conf_file)
    base_name_output_file = ep.name_tag * "_bin-" * string(ep.Bin)

    dataX, dataY, dataChar = load_npy_data(ep)

    if isdir(ep.data_folder * ep.export_folder)
        println("BE AWARE folder already exixts!!! \n")
    else
        mkdir(ep.data_folder * ep.export_folder)
        println("making dir ", ep.data_folder * ep.export_folder)
    end

    println(ep.pha_amp_com)


    if ep.pha_amp_com[1] == "_pha" && ep.pha_amp_com[2] == "_pha"
        TE_each_row(dataX, dataX, dataChar, ep, base_name_output_file * "_pha_pha")
        TE_data_rows_surrogates(base_name_output_file, dataChar, dataX, "_SePha", "-Pha", ep)

    elseif ep.pha_amp_com[1] == "_pha" && ep.pha_amp_com[2] == "_amp"
        TE_each_row(dataX, dataY, dataChar, ep, base_name_output_file * "_pha_amp")
        TE_data_rows_surrogates(base_name_output_file, dataChar, dataX, "_SePha", "-Amp", ep)
    
    elseif ep.pha_amp_com[1] == "_amp" && ep.pha_amp_com[2] == "_amp"
        TE_each_row(dataY, dataY, dataChar, ep, base_name_output_file)
        TE_data_rows_surrogates(base_name_output_file, dataChar, dataX, "_SeAmp", "-Amp", ep)
    
    elseif ep.pha_amp_com[1] == "_amp" && ep.pha_amp_com[2] == "_pha"
        TE_each_row(dataY, dataY, dataChar, ep, base_name_output_file)
        TE_data_rows_surrogates(base_name_output_file, dataChar, dataX, "_SeAmp", "-Pha", ep)
    end


    #TE_data_rows_surrogates(base_name_output_file, dataChar, dataX, "_SePha", "-Pha", ep)

end

#MI_each_delay(dataX[1, :], dataY[1, :], output_name, τ_range, (0, 1))



