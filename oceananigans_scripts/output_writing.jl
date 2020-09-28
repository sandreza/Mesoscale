if write_slices
    simulation.output_writers[:surface] = surface_output_writer
    simulation.output_writers[:middepth] = middepth_output_writer
    simulation.output_writers[:zonal] = zonal_output_writer
    simulation.output_writers[:meridional] = meridional_output_writer
end
if write_zonal
    simulation.output_writers[:statistics] = zonal_statistics
end
simulation.output_writers[:checkpoint] = checkpointer