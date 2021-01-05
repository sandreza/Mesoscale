function ic!(model, filename; ArrayType = Array)
    file = jldopen(filename)
    b = ArrayType(file["tracers"]["b"]["data"])
    u = ArrayType(file["velocities"]["u"]["data"])
    v = ArrayType(file["velocities"]["v"]["data"])
    w = ArrayType(file["velocities"]["w"]["data"])
    model.tracers.b.data[:] .= view(b, : )
    model.velocities.u.data[:] .= view(u, : )
    model.velocities.v.data[:] .= view(v, : )
    model.velocities.w.data[:] .= view(w, : )
    close(file)
    return nothing
end

function archarray(arch::CPU)
    return Array
end
function archarray(arch::GPU)
    return CuArray
end

function getlatest(filename)
    files = readdir()
    bools = [occursin(filename*"_checkpoint", file)*occursin("iteration", file) for file in files]
    ffiles = files[bools]
    bools2 = [filename == ffile[1:length(filename)] for ffile in ffiles]
    ffiles2 = ffiles[bools2]
    iterations = [parse(Int, split(ffile, "iteration")[2][1:end-length(".jld2")]) for ffile in ffiles2]
    p = sortperm(iterations)
    return latest_file = ffiles2[p[end]]
end