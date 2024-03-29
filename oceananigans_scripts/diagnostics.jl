# Slices
if write_slices
fields = Dict(
    "u" => model.velocities.u,
    "v" => model.velocities.v,
    "w" => model.velocities.w,
    "b" => model.tracers.b
)

surface_output_writer =
    NetCDFOutputWriter(model, fields, filename= filename_1 * "_surface.nc",
			           time_interval= slice_output_interval, zC=Nz, zF=Nz)

middepth_output_writer =
    NetCDFOutputWriter(model, fields, filename= filename_1 * "_middepth.nc",
                       time_interval = slice_output_interval, zC=Int(Nz/2), zF=Int(Nz/2))

zonal_output_writer =
    NetCDFOutputWriter(model, fields, filename= filename_1 * "_zonal.nc",
                       time_interval=slice_output_interval, yC=Int(Ny/2), yF=Int(Ny/2))

meridional_output_writer =
    NetCDFOutputWriter(model, fields, filename= filename_1 * "_meridional.nc",
                       time_interval=slice_output_interval, xC=Int(Nx/2), xF=Int(Nx/2))
end
## Zonal and Time Averages
if write_zonal
u, v, w = model.velocities
b = getproperty(model.tracers, :b)

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
                     time_averaging_window = time_avg_window,
                             time_interval = zonal_output_interval,
                                    prefix = filename_1 * "_zonal_averages",
                                     force = true)
end

## Checkpointer
checkpointer = Checkpointer(model, prefix = filename_1 * "_checkpoint", time_interval = checkpoint_interval, force = true)