# obtain relevant flux gradient relations
using Statistics
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
    
        # redo the average
        local cy = avgy(avgz(cy))
        local cz = avgy(avgz(cz))
    
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

"""
loop through cases and grab all the associated flux-gradient
"""
function flux_gradient_cases(cases)
    cs = []
    cys = []
    czs = []
    vcps = []
    wcps = []
    tracer_string_list = []

    for case in cases
        println("currently on case ", case)
        if case[1:5] == "trial"
            i = Meta.parse(case[6:end])
            jlist = [0, 1, 2]
            klist = [4 * i + 0, 4 * i + 1, 4 * i + 2, 4 * i + 3]
        elseif case[1:7] == "attempt"
            i = Meta.parse(case[8:end])
            jlist = [2 * i + 3, 2 * i + 4]
            klist = [1, 2, 3, 4, 5, 6]
        end

        # choose file to loads
        prefix = "relaxation_channel_tracers_restarted_smooth_forcing_case_" * case * "_averages.jld2"

        # pick out tracer string names based on file convention
        tracer_strings = []
        for j in jlist, k in klist
            push!(tracer_strings, "c_j" * string(j) * "_k" * string(k))
        end

        # Grab Tracers for analysis
        jl_file = jldopen(prefix, "r+")

        l_cs, l_cys, l_czs, l_vcps, l_wcps, l_∇c∇bs, l_∇c∇ᵖbs, l_u⃗c∇bs, l_u⃗c∇ᵖbs = c_flux_gradient(jl_file, tracer_strings = tracer_strings)

        for i in eachindex(l_cys)
            push!(cs, l_cs[i])
            push!(cys, l_cys[i])
            push!(czs, l_czs[i])
            push!(vcps, l_vcps[i])
            push!(wcps, l_wcps[i])
            push!(tracer_string_list, tracer_strings)
        end
    end
    return cs, cys, czs, vcps, wcps, tracer_string_list
end

function eddy_kinetic_energy(u, v, w, uu, vv, ww)
    eke = 0.5 .* (uu - u .* u + vv - v .* v + ww - w .* w)
    return eke
end

"""
bouyancy and gradients 
bouyancy fluxes 
eddy-kinetic energy 
zonal velocity
meriodional velocity 
vertical velocity
"""
function physical_fields()
    prefix = "relaxation_channel_tracers_restarted_smooth_forcing_case_" * "trial5" * "_averages.jld2"
    jl_file = jldopen(prefix, "r+")
    time_index = 22
    b = get_field("b", time_index, jl_file)
    y, z = get_grid(b, jl_file; ghost = 3)

    vb = avgy(get_field("vb", time_index, jl_file))
    wb = avgz(get_field("wb", time_index, jl_file))

    bz = ∂z(b, z)
    by = ∂y(b, y)

    v = get_field("v", time_index, jl_file)
    w = get_field("w", time_index, jl_file)

    uu = get_field("uu", time_index, jl_file)
    u = get_field("u", time_index, jl_file)

    vv = avgy(get_field("vv", time_index, jl_file))
    v = avgy(get_field("v", time_index, jl_file))

    ww = avgz(get_field("ww", time_index, jl_file))
    w = avgz(get_field("w", time_index, jl_file))

    vpbp = vb - v .* b # v prime b prime 
    wpbp = wb - w .* b # w prime b prime
    eke = eddy_kinetic_energy(u, v, w, uu, vv, ww)
    return b, by, bz, u, v, w, eke, vpbp, wpbp, y, z
end

function symmetric_quantiles(field, quant; symmetrize = false)
    @assert 0 <= quant <= 1
    quant = quant >= 0.5 ? 1-quant : quant
    vals = quantile.(Ref(field[:]), [quant, 1-quant])
    if symmetrize
        @assert vals[1] <= 0.0
        @assert vals[2] >= 0.0
        maxval = maximum(abs.(vals))
        return [-maxval, maxval]
    end
    return vals
end