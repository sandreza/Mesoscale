# Slices
fields = Dict(
    "u" => model.velocities.u,
    "v" => model.velocities.v,
    "w" => model.velocities.w,
    "b" => model.tracers.b
)

surface_output_writer =
    NetCDFOutputWriter(model, fields, filename= filename_1 * "_surface.nc",
			           time_interval=output_interval, zC=Nz, zF=Nz)

middepth_output_writer =
    NetCDFOutputWriter(model, fields, filename= filename_1 * "_middepth.nc",
                       time_interval = output_interval, zC=Int(Nz/2), zF=Int(Nz/2))

zonal_output_writer =
    NetCDFOutputWriter(model, fields, filename= filename_1 * "_zonal.nc",
                       time_interval=output_interval, yC=Int(Ny/2), yF=Int(Ny/2))

meridional_output_writer =
    NetCDFOutputWriter(model, fields, filename= filename_1 * "_meridional.nc",
                       time_interval=output_interval, xC=Int(Nx/2), xF=Int(Nx/2))

## Checkpointer
checkpointer = Checkpointer(model, prefix = filename_1 * "_checkpoint", time_interval = checkpoint_interval, force = true)

## Zonal and Time Averages

u, v, w = model.velocities
b = model.tracers
u_scratch = XFaceField(model.architecture, model.grid)
v_scratch = YFaceField(model.architecture, model.grid)
w_scratch = ZFaceField(model.architecture, model.grid)
b_scratch =  CellField(model.architecture, model.grid)

zonal_fields = Dict(
    :u => AveragedField(u, dims=(1,), computed_data=u_scratch.data),
    :v => AveragedField(v, dims=(1,), computed_data=v_scratch.data),
    :w => AveragedField(w, dims=(1,), computed_data=w_scratch.data),
    :b => AveragedField(b, dims=(1,), computed_data=b_scratch.data),

    :uu => AveragedField(u * u, dims=(1,), computed_data=u_scratch.data),
    :vv => AveragedField(v * v, dims=(1,), computed_data=v_scratch.data),
    :ww => AveragedField(w * w, dims=(1,), computed_data=w_scratch.data),
    :uv => AveragedField(u * v, dims=(1,), computed_data=u_scratch.data),
    :vw => AveragedField(w * v, dims=(1,), computed_data=w_scratch.data),
    :uw => AveragedField(w * u, dims=(1,), computed_data=w_scratch.data),
    :ub => AveragedField(u * b, dims=(1,), computed_data=b_scratch.data),
    :vb => AveragedField(v * b, dims=(1,), computed_data=b_scratch.data),
    :wb => AveragedField(w * b, dims=(1,), computed_data=b_scratch.data),
)

zonal_statistics = JLD2OutputWriter(model, merge(turbulence_statistics, tke_budget_statistics),
                     time_averaging_window = time_avg_window,
                             time_interval = output_interval,
                                    prefix = filename_1 * "_surface.nc",
                                     force = true)