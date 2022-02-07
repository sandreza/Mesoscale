# obtain relevant flux gradient relations
function c_flux_gradient(jld2_file; tracer_strings = ["c1", "c2", "c3", "c4"], time_index = 22)
    jl_file = jld2_file
    b = get_field("b", time_index, jl_file)
    y, z = get_grid(b, jl_file; ghost = 3)
    y̅ = avg(y)
    z̅ = avg(z)
    bz_p = ∂z(b, z)
    by_p = ∂y(b, y)

    bz = avgy(bz_p)
    by = avgz(by_p)

    v = get_field("v", time_index, jl_file)
    w = get_field("w", time_index, jl_file)
    cs = []
    cys = []
    czs = []
    vcps = []
    wcps = []

    ∇c∇ᵖbs = []  # isopycnal tracer gradients ∇c ⋅ ∇ᵖb
    ∇c∇bs = []  # diapycnal tracer gradients ∇c ⋅ ∇b
    u⃗c∇ᵖbs = [] # isopycnal eddy flux u⃗'c' ⋅ ∇ᵖb
    u⃗c∇bs = [] # diapycnal eddy flux u⃗'c' ⋅ ∇b

    for tracer_string in tracer_strings
        local c = get_field(tracer_string, time_index, jl_file)
        local vc = get_field("v" * tracer_string, time_index, jl_file)
        local wc = get_field("w" * tracer_string, time_index, jl_file)

        local vcp = avgy(avgz(avgy(vc) - avgy(v) .* c))
        local wcp = avgy(avgz(avgz(wc) - avgz(w) .* c))

        local cz = avgy(∂z(c, z))
        local cy = avgz(∂y(c, y))

        local vcp_avg = avgz(avgy(vcp))
        local wcp_avg = avgz(avgy(wcp))

        local ∇c∇b = cy .* by + cz .* bz
        local ∇c∇ᵖb = cy .* bz - cz .* by

        local u⃗c∇b = vcp .* by + wcp .* bz
        local u⃗c∇ᵖb = vcp .* bz - wcp .* by

        push!(cs, c)
        push!(cys, cy)
        push!(czs, cz)
        push!(vcps, vcp_avg)
        push!(wcps, wcp_avg)

        push!(∇c∇bs, ∇c∇b)
        push!(∇c∇ᵖbs, ∇c∇ᵖb)
        push!(u⃗c∇bs, u⃗c∇b)
        push!(u⃗c∇ᵖbs, u⃗c∇ᵖb)
    end
    return cs, cys, czs, vcps, wcps, ∇c∇bs, ∇c∇ᵖbs, u⃗c∇bs, u⃗c∇ᵖbs
end