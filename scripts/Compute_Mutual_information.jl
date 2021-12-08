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