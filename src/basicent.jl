#using Pkg; Pkg.activate("../")
#Pkg.add("Entropies")
#Pkg.add("DelayEmbeddings")
#Pkg.add(url="https://github.com/JuliaDynamics/Entropies.jl")


using DelayEmbeddings
using Entropies
export tranfserentropy


"""
tranfserentropy(joint, est::Kraskov,  α =1.0, base =2)

Estimate transfer (generalized orther α =1) for a joint contruct using a default Kraskov
nearest neighbours aproach to estimate provabilitues, it is not able to process α != 1

"""



function tranfserentropy(joint, est; α =1.0, base =2)
    if length(joint[1]) == 3
        H3 = Entropies.genentropy(joint, est)
        H2a =Entropies.genentropy(joint[:,[3,2]], est)
        H2b=Entropies.genentropy(joint[:,[1,2]], est)
        H1=Entropies.genentropy(joint[:,[2]], est)
        return -H3+H2a+H2b-H1
    end
    if length(joint[1]) == 5
        H3 = Entropies.genentropy(joint, est)
        H2a = Entropies.genentropy(joint[:,[5,2, 3, 4]], est)
        H2b = Entropies.genentropy(joint[:,[1,2, 3, 4]], est)
        H1 = Entropies.genentropy(joint[:,[2, 3, 4]], est)
    return -H3+H2a+H2b-H1
    end

end



