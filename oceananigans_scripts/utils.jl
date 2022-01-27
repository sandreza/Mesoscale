function get_grid(field, jl_file; ghost = 3)
    
    gridstrings = keys(jl_file["grid"])

    yC_string = "yC" in gridstrings ? "yC" : "yᵃᶜᵃ" 
    zC_string = "zC" in gridstrings ? "zC" : "zᵃᵃᶜ"
    yF_string = "yF" in gridstrings ? "yF" : "Δyᵃᶠᵃ"
    zF_string = "zF" in gridstrings ? "zF" : "zᵃᵃᶠ"

    yC = collect(jl_file["grid"][yC_string])[ghost+1:end-ghost]
    zC = collect(jl_file["grid"][zC_string])[ghost+1:end-ghost]

    yF = collect(jl_file["grid"][yF_string])[ghost+1:end-ghost]
    zF = collect(jl_file["grid"][zF_string])[ghost+1:end-ghost]

    y = size(field)[1] == size(yC)[1] ? yC : yF
    z = size(field)[2] == size(zC)[1] ? zC : zF

    return y, z
end

function get_grid2(field, jl_file; ghost = 3)

    yC = collect(jl_file["grid"]["yᵃᶜᵃ"])[ghost+1:end-ghost]
    zC = collect(jl_file["grid"]["zᵃᵃᶜ"])[ghost+1:end-ghost]

    yF = collect(jl_file["grid"]["yᵃᶠᵃ"])[ghost+1:end-ghost]
    zF = collect(jl_file["grid"]["zᵃᵃᶠ"])[ghost+1:end-ghost]


    y = size(field)[1] == size(yC)[1] ? yC : yF
    z = size(field)[2] == size(zC)[1] ? zC : zF

    return y, z
end

function get_field(field_string::String, time_index::Int, jl_file::JLD2.JLDFile; ghost = 3)
    t_keys = keys(jl_file["timeseries"][field_string])
    field = jl_file["timeseries"][field_string][t_keys[time_index]]
    return field[1,:,:]
end

function ∂y(field, y)
    δf = field[2:end, :] - field[1:end-1, :]
    δy = y[2:end] - y[1:end-1]
    δy = reshape(δy[:], (length(δy), 1))
    return δf ./ δy
end

function ∂z(field,z)
    δf = field[:, 2:end] - field[:, 1:end-1]
    δz = z[2:end] - z[1:end-1]
    δz = reshape(δz[:], (1, length(δz)))
    return δf ./ δz
end

function avg(field)
    return (field[2:end] + field[1:end-1]) .* 0.5
end
function avgy(field)
    return (field[2:end, :] + field[1:end-1, :]) .* 0.5
end
function avgz(field)
    return (field[:, 2:end] + field[:, 1:end-1]) .* 0.5
end
