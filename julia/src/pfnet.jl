module pfnet

# Dependencies
if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
    include("../deps/deps.jl")
else
    error("PFNET is not properly installed. Please run Pkg.build(\"pfnet\")")
end

export all

# Includes
include("net.jl")
include("parser.jl")

end
