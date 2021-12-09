using Pkg#; Pkg.activate("./")
Pkg.add("InformationMeasures")


function load_npy_data(ep)

    datas = Array{AbstractArray}(undef, length(ep.sufixes), 1)

    for (i, s) in enumerate(ep.sufixes)

        name = ep.data_folder * ep.data_name * s
        datas[i] = npzread(name)

    end

    return datas[1], datas[2], datas[3]

end
dataX, dataY, dataChar = load_npy_data(ep)
mi_12 = get_mutual_information(data_1, data_2)

output_name = name_output_file = ep.root * ep.export_folder * base_name_output_file * ".txt"#"./data/output/delay_rossler_phase.txt"


TE_each_delay(dataX[1, :], dataY[1, :], output_name, τ_range, τ_delays, emb_dim, estimator)