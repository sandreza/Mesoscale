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