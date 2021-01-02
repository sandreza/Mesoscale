"""
grabstates(filename)

# Description
Constructs states for visualisation automatically from jld2 file

# Arguments
-  `filename`: a string that corresponds to a filename

# Return states, statenames
- `states`: array of arrays of data
- `statenames`: string names for states


"""
function grabstates(filename)
    file = jldopen(filename)
    grid = file["grid"]
    xS, yS, zS = size(grid)
    btmp = file["tracers"]["b"]["data"]
    ghost = Int((length(btmp[:,1,1])-xS)/2)
    b = file["tracers"]["b"]["data"][ghost:(xS + ghost - 1), ghost:(yS + ghost - 1), ghost:(zS + ghost - 1)]
    u = file["velocities"]["u"]["data"][ghost:(xS + ghost - 1), ghost:(yS + ghost - 1), ghost:(zS + ghost - 1)]
    v = file["velocities"]["v"]["data"][ghost:(xS + ghost - 1), ghost:(yS + 1 + ghost - 1), ghost:(zS + ghost - 1)]
    w = file["velocities"]["w"]["data"][ghost:(xS + ghost - 1), ghost:(yS + ghost - 1), ghost:(zS + 1 + ghost - 1)]
    x = grid.xC[1:xS]
    y = grid.yC[1:yS]
    z = grid.zC[1:zS]
    f = file["coriolis"].f₀
    β = file["coriolis"].β
    # Cell-Centered Values
    v = (v[:, 1:end-1, :] + v[:, 2:end, :]) .* 0.5
    w = (w[:, :, 1:end-1] + w[:,:, 2:end]) .* 0.5
    close(file)
    ## Create Gradient
    include(pwd() * "/analysis_scripts/" * "post_analysis.jl") # Gradients etc. here
    ∂x = PartialDerivative((x, 1))
    ∂y = PartialDerivative((y, 2))
    ∂z = PartialDerivative((z, 3))
    ∇ = [∂x, ∂y, ∂z]
    ω = ∇ × [u, v, w]
    ∇b = [∂x(b), ∂y(b), ∂z(b)]
    coriolis =  f .+ β * reshape(y, (1, length(y) ,1)) 
    fpβ = coriolis[1:end-1] + coriolis[2:end]  # f plus beta
    ω∇b = ω⋅∇b
    pv = ω∇b + fpβ .* ∇b[3]
    alignment = ω∇b ./ ( sqrt.(ω⋅ω) .* sqrt.(∇b⋅∇b) )
    # Hydrostatic Hydrostatic hydrostatic_pressure
    Δz = z[2] - z[1]
    p = cumsum(-b .* Δz, dims = 3) # since evenly spaced
    p = p .- mean(p)
    # Deformation Radius
    relu(x) = x ≥ 0 ? x : 0.0
    Lᵣ = -sqrt.(relu.(∇b[3])) ./ fpβ * abs(z[1]-z[end]) # (NH/f)
    # check Thermal Wind
    ∂xp = ∂x(p) 
    fv = v .* coriolis
    ∂yp = ∂y(p) 
    fu = u .* coriolis
    deviationx = ∂xp .+ (fv[1:end-1,1:end-1,1:end-1] + fv[2:end,2:end,2:end]) * 0.5 # off by a gridcell (should avg)
    deviationy = ∂yp .- (fu[1:end-1,1:end-1,1:end-1] + fu[2:end,2:end,2:end]) * 0.5 
    velocity = [u, v, w]
    speed = sqrt.(velocity ⋅ velocity)
    states = [u, abs.(u), v, w, b, -p, 
              speed, ω[3], ω[1], ω[2], abs.(ω∇b), abs.(pv), 
              ∇b[3], Lᵣ, -∂xp, fv, -∂yp, -fu, 
              deviationx, deviationy]
    statenames =  ["u", "|u|", "v", "w", "b", "-p/ρ₀ (hydrostatic)",
                  "speed", "ω₃", "ω₁", "ω₂", "|ω⋅∇b|",  "|Ertel PV|", 
                  "∂z(b)", "Lᵣ", "-∂x(p/ρ₀)", "fv", "-∂y(p/ρ₀)", "-fu", 
                  "deviationx", "deviationy"]
    units = [" [m/s]", " [m/s]", " [m/s]", " [m/s]", " [m/s²]", " [m²/s²]",
              " [m/s]", " [1/s]", " [1/s]", " [1/s]", " [1/s³]", " [1/s³]",
              " [m/s²]", " [m]", " [m/s²]", " [m/s²]", " [m/s²]", " [m/s²]",
              " [m/s²]", " [m/s²]"]
    return states, statenames, units
end