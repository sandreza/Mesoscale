filename = files[end]
states, statenames, units1 = grabstates(filename)
title = "Coarsened " * grabtitle(filename)
grid = getgrid(filename)

filename2 = files[1]
states2, statenames2, units2 = grabstates(filename2)
title2 = "Coarse " * grabtitle(filename2)
coarsesize = size(states2[1])
topo = (Periodic, Bounded, Bounded)
domain = (x=(0, 1e6), y=(0, 1e6), z=(-3000, 0))
newgrid = RegularCartesianGrid(topology=topo, size= coarsesize; domain...)

coarsestates = interpolate(states, grid, newgrid)
scene = visualize(states, states2, statenames = statenames, statenames2 = statenames2, aspect = (1,1, 32/192), statistics = true, title = title, title2 = title2, units1 = units1, units2 = units2)