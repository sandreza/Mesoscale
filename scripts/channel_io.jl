## Progress printer
Ni = 1000

function print_progress(simulation)
    model = simulation.model
    i, t = model.clock.iteration, model.clock.time

    progress = 100 * (model.clock.time / end_time)

    umax = maximum(abs, model.velocities.u.data.parent)
    vmax = maximum(abs, model.velocities.v.data.parent)
    wmax = maximum(abs, model.velocities.w.data.parent)
    meanb = mean(view(model.tracers.b.data, 2:Nx-1, 2:Ny-1, 2:Nz-1))

    @printf("[%05.2f%%] i: %d, t: %.2e days, umax: (%6.3e, %6.3e, %6.3e) m/s, CFL: %6.4e, next Δt: %.1e s\n",
            progress, i, t / day, umax, vmax, wmax, cfl(model), Δt_wizard.Δt)
    println(" ")
    @printf("The mean value of b is %.2e", meanb)
    println(" ")
end

## Diagnostics
if write_slices
fields = Dict(
    "u" => model.velocities.u,
    "v" => model.velocities.v,
    "w" => model.velocities.w,
    "b" => model.tracers.b
)

surface_output_writer =
    NetCDFOutputWriter(model, fields, filepath = filename_1 * "_surface.nc",
			           schedule=TimeInterval(slice_output_interval), field_slicer = FieldSlicer(k = Nz))

middepth_output_writer =
    NetCDFOutputWriter(model, fields, filepath = filename_1 * "_middepth.nc",
                       schedule=TimeInterval(slice_output_interval), field_slicer = FieldSlicer(k = Int(floor(Nz/2))))

zonal_output_writer =
    NetCDFOutputWriter(model, fields, filepath = filename_1 * "_zonal.nc",
                       schedule=TimeInterval(slice_output_interval), field_slicer = FieldSlicer(j = Int(floor(Ny/2))))

meridional_output_writer =
    NetCDFOutputWriter(model, fields, filepath = filename_1 * "_meridional.nc",
                       schedule=TimeInterval(slice_output_interval), field_slicer = FieldSlicer(i = Int(floor(Nx/2))))
function debug_f(debug)
    if debug
        close(surface_output_writer)
        close(middepth_output_writer)
        close(zonal_output_writer)
        close(meridional_output_writer)
    end
    return nothing
end
end

## Zonal and Time Averages
if write_zonal
u, v, w = model.velocities
b = model.tracers.b

u_scratch = XFaceField(model.architecture, model.grid)
v_scratch = YFaceField(model.architecture, model.grid)
w_scratch = ZFaceField(model.architecture, model.grid)
b_scratch =  CellField(model.architecture, model.grid)

zonal_fields = Dict(
    :u => AveragedField( u, dims=(1,)),
    :v => AveragedField( v, dims=(1,)),
    :w => AveragedField( w, dims=(1,)),
    :b => AveragedField( b, dims=(1,)),
    :uu => AveragedField(u * u, dims= 1),
    :vv => AveragedField(v * v, dims= 1),
    :ww => AveragedField(w * w, dims= 1),
    :uv => AveragedField(u * v, dims= 1),
    :vw => AveragedField(w * v, dims= 1),
    :uw => AveragedField(w * u, dims= 1),
    :ub => AveragedField(u * b, dims= 1),
    :vb => AveragedField(v * b, dims= 1),
    :wb => AveragedField(w * b, dims= 1),
)

zonal_statistics = JLD2OutputWriter(model, zonal_fields,
                                    schedule = AveragedTimeInterval(time_avg_window, window=zonal_output_interval, stride=5),
                                    prefix = filename_1 * "_zonal_averages", force = true)
end
## Checkpointer
checkpointer = Checkpointer(model, prefix = filename_1 * "_checkpoint", 
                            schedule = TimeInterval(checkpoint_interval),  
                            force = true)
