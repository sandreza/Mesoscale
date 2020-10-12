searchdir(path, key) = filter(x -> occursin(key, x), readdir(path))
checkpoints = searchdir(pwd(), filename_1 * "_checkpoint_iteration")
if length(checkpoints) > 0
    checkpointed = true
    checkpoint_filepath = joinpath(pwd(), checkpoints[end])
    @info "Restoring from checkpoint: $checkpoint_filepath"
    model = restore_from_checkpoint(checkpoint_filepath, boundary_conditions = bcs, forcing = forcings,timestepper = :RungeKutta3, advection = advection_scheme)
else
    checkpointed = false
	model = IncompressibleModel(
	           architecture = arch,
	             float_type = FT,
	                   grid = grid,
	               coriolis = coriolis,
	               buoyancy = buoyancy,
	                closure = closure,
	                tracers = (:b,),
        boundary_conditions = bcs,
                advection = advection_scheme,
                timestepper = timestepping_scheme,
                    forcing = forcings,
	)
end

if !checkpointed
    if !geostrophic_balance
        set!(model, b=B₀)
    else
        set!(model, u = U₀, b=B₀)
    end
end

if !checkpointed
	Δt= 120.0
else
	Δt = 300.0
end
