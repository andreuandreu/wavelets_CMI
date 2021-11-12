using Printf

function read_bin_couplings(root = "../data/", n_points = 131072, couplings = LinRange(0, 0.25, 99))
    "
    read  and return data according to a specific coupling factor of the rossler sistem 
    specified by the coupling vector 

    :param:
        root: folder direction
        n_points: dimension of the data
        couplings: vector desciving the strength of the vector for each dataseries
    "
    all = Array{Float64}(undef,  n_points, 2, length(couplings));
    aa = Array{Float64}(undef,  n_points, 2)
    
    for (i, cc) in enumerate(couplings[1:length(couplings)-1])
        aux = @sprintf("%.3f.bin", cc)
        read!("$root$aux", aa)
        all[:,:,i] = aa
        println(i, ' ', cc)
    end

    return all
end
