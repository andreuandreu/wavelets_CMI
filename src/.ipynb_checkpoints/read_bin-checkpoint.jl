using Printf

function read_bin_couplings(root, n_points, couplings)

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