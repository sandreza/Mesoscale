using JLD2
using GLMakie

tracer_case_j = 1
tracer_case_k = 1
# file output
prefix = "fluxernathy_tracers_iteration9986400.jld2"
prefix = "fluxernathy_tracers_restarted_j1_k1_averages.jld2"
jl_file = jldopen(prefix, "r+")


function get_grid(field, jl_file; ghost = 3)
    #=
    xC = jl_file["grid"]["xC"][ghost+1:end-ghost]
    xF = jl_file["grid"]["yF"][ghost+1:end-ghost]
    yC = jl_file["grid"]["yC"][ghost+1:end-ghost]
    zC = jl_file["grid"]["zC"][ghost+1:end-ghost]
    yF = jl_file["grid"]["yF"][ghost+1:end-ghost]
    zF = jl_file["grid"]["zF"][ghost+1:end-ghost]
    =#
    xC = collect(jl_file["grid"].xC)[ghost+1:end-ghost]
    yC = collect(jl_file["grid"].yC)[ghost+1:end-ghost]
    zC = collect(jl_file["grid"].zC)[ghost+1:end-ghost]

    xF = collect(jl_file["grid"].xF)[ghost+1:end-ghost]
    yF = collect(jl_file["grid"].yF)[ghost+1:end-ghost]
    zF = collect(jl_file["grid"].zF)[ghost+1:end-ghost]


    x = size(field)[1] == size(xC)[1] ? xC : xF
    y = size(field)[2] == size(yC)[1] ? yC : yF
    z = size(field)[3] == size(zC)[1] ? zC : zF

    return x,y,z
end

function get_field(field_string::String, jl_file::JLD2.JLDFile; ghost = 3)
    if field_string == "u" || field_string == "v" || field_string == "w"
        field = jl_file["velocities"][field_string]["data"]
    else
        field = jl_file["tracers"][field_string]["data"]
    end
    return field[ghost+1:end-ghost,ghost+1:end-ghost,ghost+1:end-ghost]
end

field = get_field("c1", jl_file)
x,y,z = get_grid(field, jl_file)

heatmap(x,y, field[:,:,end], colormap = :balance, interpolate = true)

Nx, Ny, Nz = size(field)
Lx = maximum(abs.(x))[end]
Ly = maximum(abs.(y))[end]
Lz = maximum(abs.(z))[end]
xsurf = range(extrema(x)..., length = Nx)
ysurf = range(extrema(y)..., length = Ny)
zsurf = range(extrema(z)..., length = Nz)
ϕsurf = field
clims = extrema(ϕsurf)
zscale = 100
fig = Figure(resolution = (1920, 1080))
ax = fig[1,1] = LScene(fig, title= "Eddying Channel")

# edge 1
ϕedge1 = ϕsurf[:,1,:]
GLMakie.surface!(ax, xsurf, zsurf .* zscale, ϕedge1, transformation = (:xz, 0),  colorrange = clims, colormap = :balance, show_axis=false)

# edge 2
ϕedge2 = ϕsurf[:,end,:]
GLMakie.surface!(ax, xsurf, zsurf .* zscale, ϕedge2, transformation = (:xz, Ly),  colorrange = clims, colormap = :balance)

# edge 3
ϕedge3 = ϕsurf[1,:,:]
GLMakie.surface!(ax, ysurf, zsurf .* zscale, ϕedge3, transformation = (:yz, 0),  colorrange = clims, colormap = :balance)

# edge 4
ϕedge4 = ϕsurf[end,:,:]
GLMakie.surface!(ax, ysurf, zsurf .* zscale, ϕedge4, transformation = (:yz, Lx),  colorrange = clims, colormap = :balance)

# edge 5
ϕedge5 = ϕsurf[:,:,1]
GLMakie.surface!(ax, xsurf, ysurf, ϕedge5, transformation = (:xy, -Lz *  zscale), colorrange = clims, colormap = :balance)


# edge 6
ϕedge6 = ϕsurf[:,:,end]
GLMakie.surface!(ax, xsurf, ysurf, ϕedge6, transformation = (:xy, 0 *  zscale), colorrange = clims, colormap = :balance)

display(fig)

