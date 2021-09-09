module Mesoscale
using Pkg
export add_branch

function add_branch()
    Pkg.add(url = "https://github.com/CliMA/Oceananigans.jl.git", rev = "glw-ncc/overturning-channel-example")
end

end # end of module