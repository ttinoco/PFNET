module pfnet

# Dependencies
if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
    include("../deps/deps.jl")
else
    error("PFNET is not properly installed. Please run Pkg.build(\"pfnet\")")
end

# Includes
include("strings.jl")
include("net.jl")
include("parser.jl")

# Symbols
for name in names(pfnet, true)
    @eval export $(Symbol(name))
end

end
