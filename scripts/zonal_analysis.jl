files = [
    pwd() * "/Channel_1_zonal_averages.jld2"
]

file = files[1]
zonalstatistics = jldopen(file)
tkeys = keys(zonalstatistics["timeseries"]["t"])
t = [zonalstatistics["timeseries"]["t"][tkey] for tkey in tkeys]
u = [zonalstatistics["timeseries"]["u"][tkey][1,:,:] for tkey in tkeys]
v = [zonalstatistics["timeseries"]["v"][tkey][1,:,:] for tkey in tkeys]
w = [zonalstatistics["timeseries"]["w"][tkey][1,:,:] for tkey in tkeys]
b = [zonalstatistics["timeseries"]["b"][tkey][1,:,:] for tkey in tkeys]
vb = [zonalstatistics["timeseries"]["vb"][tkey][1,:,:] for tkey in tkeys]
vv = [zonalstatistics["timeseries"]["vv"][tkey][1,:,:] for tkey in tkeys]
vw = [zonalstatistics["timeseries"]["vw"][tkey][1,:,:] for tkey in tkeys]
ub = [zonalstatistics["timeseries"]["ub"][tkey][1,:,:] for tkey in tkeys]
close(zonalstatistics)
##
heatmap(vw[1])
