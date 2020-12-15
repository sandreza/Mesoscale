using NetCDF, Plots, GLMakie, AbstractPlotting
using Printf, Statistics
using ImageTransformations, Colors, JLD2
using AbstractPlotting.MakieLayout

filename = pwd() * "/Weno_20_checkpoint_iteration21030963.jld2"
file = jldopen( filename )
b = file["tracers"]["b"]["data"][2:end-1, 2:end-1, 2:end-1]
u = file["velocities"]["u"]["data"][2:end-1, 2:end-1, 2:end-1]
v = file["velocities"]["v"]["data"][2:end-1, 2:end-1, 2:end-1]
w = file["velocities"]["w"]["data"][2:end-1, 2:end-1, 2:end-1]
close(file)
nx, ny, nz = size(u)
v = (v[:,1:end-1,:] + v[:, 2:end, :]) * 0.5
w = (w[:, :, 1:end-1] + w[:, :, 2:end]) * 0.5
states = [u, v, w, b]
##
Φ = b
x, y, z = size(Φ)
clims = (minimum(Φ), maximum(Φ))
cmap_rgb = to_colormap(:thermometer)
volume(0..x, 0..y, 0..z, Φ, camera = cam3d!, colormap = cmap_rgb, colorrange = clims)

##
function avgx(Φ, n)
    m = size(Φ)[1]
    scale = Int(floor(m/n))
    if ( abs(Int(floor(m/n)) - m/n) > eps(1.0))
        return error
    end
    nx, ny, nz = size(Φ)
    Φ2 = zeros(n, ny, nz)
    for i in 1:n
        Φ2[i,:,:] .= 0
            for j in 1:scale
                Φ2[i,:,:] .+= Φ[scale*(i-1) + j,:,:] / scale
            end
    end
    return Φ2
end
function avgy(Φ, n)
    m = size(Φ)[2]
    scale = Int(floor(m/n))
    if ( abs(Int(floor(m/n)) - m/n) > eps(1.0))
        return error
    end
    nx, ny, nz = size(Φ)
    Φ2 = zeros(nx, n, nz)
    for i in 1:n
        Φ2[:,i,:] .= 0
            for j in 1:scale
                Φ2[:,i,:] .+= Φ[:,scale*(i-1) + j,:] / scale
            end
    end
    return Φ2
end
avgxy(Φ, n) = avgx(avgy(Φ, n),n )
##
# need prime factorization of points in the horizontal direction
coarsenings = sort([2^i * 3^j for i in 0:6 for j in 0:1], rev = true) 
##

stateindex = [1, 2, 3, 4]
statenode = Node(stateindex[4])
statenames = ("u", "v", "w", "b")
state = @lift(states[$statenode])
scene, layout = layoutscene()
lscene = layout[1:4, 2:4] = LScene(scene)
# lscene = layout[3:4, 2:4] = LScene(scene)

coarsening = Node(coarsenings[1])
Φ = @lift(avgxy(states[$statenode], $coarsening))
x, y, z = size(b) #@lift(size(states[$statenode])) #size(state)
clims = @lift((minimum(states[$statenode]), maximum(states[$statenode]))) # (minimum(state), maximum(state))
cmap_rgb = to_colormap(:thermometer)

supertitle = layout[1,2] = LText(scene, " "^10 * " Field and Coarse Graining " * " "^10, textsize = 50, color = :black)

volume!(lscene, 0..x, 0..y, 0..z, Φ, 
        camera = cam3d!, 
        colormap = cmap_rgb, 
        colorrange = clims)

volume!(lscene, x..2x, 0..y, 0..z, state, 
        camera = cam3d!, 
        colormap = cmap_rgb, 
        colorrange = clims)

menu = LMenu(scene, options = zip(coarsenings, coarsenings))
statemenu = LMenu(scene, options = zip(statenames, stateindex))

on(menu.selection) do s
    coarsening[] = s
end

on(statemenu.selection) do s
    statenode[] = s
end

layout[1, 1] = vgrid!(
    LText(scene, "Coarseness", width = nothing),
    menu,
    LText(scene, "State", width = nothing),
    statemenu
)

display(scene)