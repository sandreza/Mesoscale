# Boiler Plate "using" in seperate file
include(pwd() * "/scripts/dependencies.jl")
include(pwd() * "/scripts/loop_names.jl")

Qlist = [10, 1, 5]
τlist = [0.2, 0.01, 0.05, 0.1]

for Q in Qlist, τ in τlist
    simulationdescription = generate_descriptor(Q, τ)
    println("Currentely doing ", simulationdescription)
    channel_looper(Q, τ)
end

