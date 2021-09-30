"""given the basic parameters, returns a string with the name of the file to store the value"""

function name_of_file(est, Bin, maxτ, τs, js)
    strEst = string(est)[1:3]
    prefix = "couplings_TE"
    suffix = "txt"
    
    tsjs = '_'
    for (i, e) in enumerate(τs)
        tsjs = string(τs[length(τs)-i+1], tsjs, js[i])
    end
    
    #name = "$prefix-$strEst-$Bin-lag_$lag.$suffix"
    d = (est = strEst,  Bin = Bin, maxlag = maxτ, tsjs = tsjs  )
    name = savename(prefix, d , suffix)
    println(name)
    return name
end