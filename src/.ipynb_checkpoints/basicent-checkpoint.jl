#using Pkg; Pkg.activate("../")
#Pkg.add("Entropies")
#Pkg.add("DelayEmbeddings")
#Pkg.add(url="https://github.com/JuliaDynamics/Entropies.jl")


using DelayEmbeddings
using Entropies
export tranfserentropy


"""
tranfserentropy(joint, est::VisitationFrequency(b),  α =1.0, base =2)

Estimate transfer (generalized orther α =1) for a joint contruct using a default visitation 
frequency aproach to estimate provabilitues

"""

function tranfserentropy(joint, est::VisitationFrequency;  α =1.0, base =2)
    #genentropy!(p, x, est; α = 1.0, base = Base.MathConstants.e)
    H3 = Entropies.genentropy(joint, est, α=α)
    H2a = Entropies.genentropy(joint[:,[3,2]], est, α=α)
    H2b = Entropies.genentropy(joint[:,[1,2]], est, α=α)
    H1 = Entropies.genentropy(joint[:,[2]], est, α=α)
    return -H3+H2a+H2b-H1
end

"""
tranfserentropy(joint, est::Kraskov,  α =1.0, base =2)

Estimate transfer (generalized orther α =1) for a joint contruct using a default Kraskov
nearest neighbours aproach to estimate provabilitues, it is not able to process α != 1

"""

#function tranfserentropy(joint, est::Kraskov,  base =2)
#    H3 = Entropies.genentropy(joint, est)
#    H2a =Entropies.genentropy(joint[:,[3,2]], est)
#    H2b=Entropies.genentropy(joint[:,[1,2]], est)
#    H1=Entropies.genentropy(joint[:,[2]], est)
#    return -H3+H2a+H2b-H1
#end

"""
tranfserentropy(joint, est::Kraskov,  α =1.0, base =2)

Estimate transfer (generalized orther α =1) for a joint contruct using a default Kraskov
nearest neighbours aproach for a embedded dimensionality 5 
it is not able to process α != 1 it is not able to process α != 1


function tranfserentropy(joint, est::Kraskov, dim::Int = 5, base =2)
    print("kra5")
    H3 = Entropies.genentropy(joint, est)
    H2a = Entropies.genentropy(joint[:,[5,2, 3, 4]], est)
    H2b = Entropies.genentropy(joint[:,[1,2, 3, 4]], est)
    H1 = Entropies.genentropy(joint[:,[2, 3, 4]], est)
    return -H3+H2a+H2b-H1
end

"""



"""
tranfserentropy(joint, est::KozachenkoLeonenko,  α =1.0, base =2)

Estimate transfer (generalized orther α =1) for a joint contruct using a KozachenkoLeonenko 
NearestNeighborEntropyEstimator aproach to estimate provabilitues, it is not able to process α != 1





function tranfserentropy(joint, est::KozachenkoLeonenko,  base =2)
    H3 = Entropies.genentropy(joint, est)
    H2a = Entropies.genentropy(joint[:,[3,2]], est)
    H2b = Entropies.genentropy(joint[:,[1,2]], est)
    H1 = Entropies.genentropy(joint[:,[2]], est)
    return -H3+H2a+H2b-H1
end
"""

"""
tranfserentropy(joint, est::KozachenkoLeonenko,  α =1.0, base =2)

Estimate transfer (generalized orther α =1) for a joint contruct using a KozachenkoLeonenko 
NearestNeighborEntropyEstimator aproach for a embedded dimensionality 5 it is not able to process α != 1





function tranfserentropy(joint, est::KozachenkoLeonenko, embdim = 5, base =2)
    print("  uh ")
    H3 = Entropies.genentropy(joint, est)
    H2a = Entropies.genentropy(joint[:,[5,2, 3, 4]], est)
    H2b = Entropies.genentropy(joint[:,[1,2, 3, 4]], est)
    H1 = Entropies.genentropy(joint[:,[2, 3, 4]], est)
   
    return -H3+H2a+H2b-H1
end

"""

"""
tranfserentropy(joint, est::VisitationFrequency; embdim::5, α =1.0, base =2)

Estimate transfer (generalized orther α =1) for a joint contruct using a VisitationFrequency
for a embedded dimensionality 5

"""

function tranfserentropy(joint, est::VisitationFrequency; embdim = 5, α =1.0, base =2)
    #genentropy!(p, x, est; α = 1.0, base = Base.MathConstants.e)
    H3 = Entropies.genentropy(joint, est, α=α)
    H2a = Entropies.genentropy(joint[:,[5,2, 3, 4]], est, α=α)
    H2b = Entropies.genentropy(joint[:,[1,2, 3, 4]], est, α=α)
    H1 = Entropies.genentropy(joint[:,[2, 3, 4]], est, α=α)
    return -H3+H2a+H2b-H1
end