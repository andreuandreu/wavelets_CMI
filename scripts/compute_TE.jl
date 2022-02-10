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
    pha_amp_more_str = retrieve(conf, "names", "pha_amp_more")

    Bin = parse(Int64, retrieve(conf, "emb_par", "bins"))
    maxτ = parse(Int64, retrieve(conf, "emb_par", "max_tau"))
    jumpτ = parse(Int64, retrieve(conf, "emb_par", "jump_tau"))
    emb_dim_str = retrieve(conf, "emb_par", "embeding_dimension")
    period_range_str = retrieve(conf, "emb_par", "period_range")

    name_tag = retrieve(conf, "prob_est", "name_tag")
    prob_est = retrieve(conf, "prob_est", "prob_kind")

    surr_kind = retrieve(conf, "surrogates", "surr_kind")

    sufixes = []#Array{String}(undef, 3, 1)
    for s in collect(sufixes_str)
        push!(sufixes, s)
    end
    pha_amp_com = []
    for e in collect(pha_amp_com_str)
        push!(pha_amp_com, e)
    end

    pha_amp_more = []
    for e in collect(pha_amp_more_str)
        push!(pha_amp_more, e)
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
        sufixes, pha_amp_com, pha_amp_more, Bin, maxτ, jumpτ, emb_dim, phase_periods,
        name_tag, prob_est, surr_kind)

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
    pha_amp_more::Vector
    Bin::Int64
    maxτ::Int64
    jumpτ::Int64
    emb_dim::Vector
    period_range::Vector
    name_tag::String
    prob_est::String
    surr_kind::String

end


"compute stuff"

function write_entropies_tau(name, entropies)


    open(name, "w") do file
        for e in entropies

            write(file, "$e\n")
        end
    end

end

function MI_each_column(dataX, dataY, file_name, ep, serieChar)

    "
    compute MI 
    "
    estimator, τ_range, τ_delays, emb_dim = inizilaize_embedings(ep, serieChar)
    ts = (1, 1)
    emb_dim = (1, 2)
    #println("oooaooo", dataY)

    open(file_name, "w") do file
        for i = 1:size(dataY, 1)
            #println("this is Y", i, "  ", dataY[i, :])
            #mutual_info = get_mutual_information(dataX, dataY[i, :])
            joint = DelayEmbeddings.genembed(Dataset(dataX, dataY[i, :]), ts, emb_dim)
            mutual_info = mutualInformation(joint, estimator; α = 1.0, base = 2)
            write(file, "$serieChar $mutual_info \n")
            #println("now doing", ' ', mutual_info)
        end
    end


end

function MI_each_row(dataX, dataY, dataChar, ep, base_name_output_file)

    "
    compute MI for each row of data (the x one) in one of the two data files compared to the other
    "

    TEdata_folder = ep.data_folder * ep.export_folder * "Dat_files/"
    make_dir_if_not_exist(TEdata_folder)

    for (i, char) in enumerate(dataChar)
    
        #println("this is X", i, "  ", dataX[i, :])
        strChar = lpad(@sprintf("%.2i", char), 2, "0")
        name_output_file = TEdata_folder * "MI_" * base_name_output_file *
                           "_row-p" * strChar * ".txt"
    
        MI_each_column(dataX[i, :], dataY, name_output_file, ep, char)
    
    
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

    TEdata_folder = ep.data_folder * ep.export_folder * "Dat_files/"
    make_dir_if_not_exist(TEdata_folder)

    for (i, char) in enumerate(dataChar)
        strChar = lpad(@sprintf("%.2i", char), 2, "0")
        name_output_file = TEdata_folder * base_name_output_file * "_row-p" * strChar * ".txt"

        TE_each_column(dataX[i, :], dataY, name_output_file, ep, char)
        if i % 20 == 3
            println("entropies have been stored in here ", name_output_file)
        end
    end


end




"compute stuff with surrogates"

function TE_surrogate_set(surrogate_set, surr_folder, output_name, dataX, ep, serieChar)

    "compute the TE of a guiven time series with a guiven set of surrogates, 
    measure the mean and store it given a pattern of x y names or characteristics"

    kind_of_ergodicity = Entropies.RectangularBinning(ep.Bin)
    estimator = VisitationFrequency(kind_of_ergodicity)

    entropies = zeros(0)
    mutual_infos = zeros(0)
    for i = 1:size(surrogate_set)[1]
        s = surrogate_set[i, :]
        if occursin("-Amp", output_name)
            compute = abs.(s[1:end-1])
        elseif occursin("-Pha", output_name)
            compute = angle.(s[1:end-1])
        else
            println("BE AWARE you are missing '_amp' or '_pha' in the name \n")
        end

        joint = DelayEmbeddings.genembed(Dataset(dataX, dataY[i, :]), (0, 0), (1, 2))
        MI = mutualInformation(joint, estimator; α = 1.0, base = 2)

        TE = TE_each_delay(dataX, compute, ep, serieChar)

        append!(mutual_infos, MI)
        append!(entropies, TE)
    end

    name_surr_arrays = surr_folder * "set-" * output_name[1:end-4] * ".npy"
    name_MI_surr_arrays = surr_folder * "MI_set-" * output_name[1:end-4] * ".npy"

    npzwrite(name_surr_arrays, entropies)
    npzwrite(name_MI_surr_arrays, mutual_infos)

    return mean(entropies), std(entropies), mean(mutual_infos), std(mutual_infos)

end

function TE_all_surrogates(names_surrogates, surr_folder, output_name, dataSerie, serieChar, ep)
    #τ_range, τ_delays, emb_dim, estimator
    "
    for each set of surrogates, compute the TE with a given data series and print it in a file
    "
    extra_name_emb = get_emb_dim_name(ep)
    strChar = lpad(@sprintf("%.2i", serieChar), 2, "0")
    #name_out_file = output_name * extra_name_emb * "_p$(@sprintf("%.2i", serieChar))" * ".txt"
    name_out_file = output_name * extra_name_emb * "_p" * strChar * ".txt"


    println("\n NANANANA ",  name_out_file)

    open( surr_folder * name_out_file, "w") do file
        for (i, n) in enumerate( names_surrogates)
            surrogate_set = npzread(ep.data_folder * ep.input_folder * n)
            one_mean, one_var, MI_mean, MI_var = TE_surrogate_set(surrogate_set, surr_folder,
             name_out_file, dataSerie, ep, serieChar)

            if i % length(names_surrogates)/2 == 1
                println("mmmmmmmm   ", n, " ", one_mean, "+-", one_var)
            end
        
            x = n[end-7:end-4]
            write(file, "$x $serieChar $one_mean $one_var $MI_mean $MI_var\n")
        end
    end
end


function TE_data_rows_surrogates(base_name_output_file, dataChar, dataSeries,
    tagSeries, amp_or_phase_surr, ep)

    "
    compute TE for each row of data for each set of surrogates with varing frequencies (periods)
    "

    surr_root = "surr_circ_"

    #amp_or_phase_surr = "-pha"

    println("NIGHTMARE in ", ep.data_folder * ep.input_folder * surr_root * ep.in_data_tag)
    names_surrogates = frequencies_names(ep.data_folder * ep.input_folder, surr_root * ep.in_data_tag)

    surr_folder = ep.data_folder * ep.export_folder * "Su-" * ep.surr_kind * "_files/"
    make_dir_if_not_exist(surr_folder)

    name_output_file = base_name_output_file * tagSeries * "_Su" * amp_or_phase_surr
    println("NIGHTMARE out ", name_output_file)

    for (i, char) in enumerate(dataChar)


        TE_all_surrogates(names_surrogates, surr_folder, name_output_file, dataSeries[i, :], char, ep)

    end


end

function TE_binary_percentil(TEmatrix, SurrMats)

    for m in SurrMats
    
        diff = TEmatrix - m
        diff[:, findall(diff .> 0)] .= 1
    
        ifelse.(m .> SurrMats, a, m .= 1)
        ifelse.(m .< SurrMats, a, m .= 0)

        mSum = mSum .+ m
    
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


function what_to_correlate(pha_amp)

    if pha_amp[1] == "_pha" && pha_amp[2] == "_pha"
        TE_each_row(dataX, dataX, dataChar, ep, base_name_output_file * "_pha_pha")
        TE_data_rows_surrogates(base_name_output_file, dataChar, dataX, "_SePha", "-Pha", ep)
        MI_each_row(dataX, dataX, dataChar, ep, base_name_output_file * "_pha_pha")
    
    elseif pha_amp[1] == "_pha" && pha_amp[2] == "_amp"
        TE_each_row(dataX, dataY, dataChar, ep, base_name_output_file * "_pha_amp")
        TE_data_rows_surrogates(base_name_output_file, dataChar, dataX, "_SePha", "-Amp", ep)
    
    elseif pha_amp[1] == "_amp" && pha_amp[2] == "_amp"
        TE_each_row(dataY, dataY, dataChar, ep, base_name_output_file)
        TE_data_rows_surrogates(base_name_output_file, dataChar, dataY, "_SeAmp", "-Amp", ep)
    
    elseif pha_amp[1] == "_amp" && pha_amp[2] == "_pha"
        TE_each_row(dataY, dataX, dataChar, ep, base_name_output_file)
        TE_data_rows_surrogates(base_name_output_file, dataChar, dataY, "_SeAmp", "-Pha", ep)
    end


end

function make_dir_if_not_exist(name_dir)
    if isdir(name_dir)
        println("BE AWARE folder already exixts!!! \n")
    else
        mkdir(name_dir)
        println("making dir ", name_dir)
    end
end






@time begin

    name_conf_file = ARGS[1]

    ep = load_embedings_params(name_conf_file)
    base_name_output_file = ep.name_tag * "_bin-" * string(ep.Bin)

    dataX, dataY, dataChar = load_npy_data(ep)

    make_dir_if_not_exist(ep.data_folder * ep.export_folder)

    println(ep.pha_amp_com)
    what_to_correlate(ep.pha_amp_com)

    #println(ep.pha_amp_more)
    #what_to_correlate(ep.pha_amp_more)



    #TE_data_rows_surrogates(base_name_output_file, dataChar, dataX, "_SePha", "-Pha", ep)

end

#MI_each_delay(dataX[1, :], dataY[1, :], output_name, τ_range, (0, 1))



