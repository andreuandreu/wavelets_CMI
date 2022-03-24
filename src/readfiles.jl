using DataFrames
using CSV
using DelimitedFiles

function read_surrogate(name, npts = 1000)
    #name = "/Users/andreu/transfer_entropy_mock_data/arosf11n00eps100raw.dat"
    #x, y = CSV.read("/Users/andreu/transfer_entropy_mock_data/arosf11n00eps100raw.dat", DataFrame, rows = npts)
    df= CSV.File(name; comment="#", delim=' ', 
        ignorerepeated=true, skipto=12, limit=npts, select=[1, 2])#,select=[1, 2]
    b = DataFrame(df)
end 