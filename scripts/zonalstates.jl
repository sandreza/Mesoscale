
function grabzonalstates(file; ghost = 3, startind = 1)
    # ghost should be 3 WENO
    zonalstatistics = jldopen(file)
    tkeys = keys(zonalstatistics["timeseries"]["t"])
    t = [zonalstatistics["timeseries"]["t"][tkey] for tkey in tkeys]
    yC = zonalstatistics["grid"]["yC"][ghost+1:end-ghost]
    zC = zonalstatistics["grid"]["zC"][ghost+1:end-ghost]
    yF = zonalstatistics["grid"]["yF"][ghost+1:end-ghost]
    zF = zonalstatistics["grid"]["zF"][ghost+1:end-ghost]
    fields = [:u, :v, :w, :b, :vb, :vv, :wb, :ww]
    for field in fields
        label = string(field)
        @eval $field = mean([zonalstatistics["timeseries"][$label][tkeys[i]][1,:,:] for i in $startind:length(tkeys)])
    end
    close(zonalstatistics)

    vpbp =  average1(vb) - average1(v) .* b
    uz = Δ2(u) ./ reshape((zC[2:end] - zC[1:end-1]), (1, length(zC)-1))
    bz = Δ2(b) ./ reshape((zC[2:end] - zC[1:end-1]), (1, length(zC)-1))
    by = Δ1(b) ./ reshape((yC[2:end] - yC[1:end-1]), (length(yC)-1, 1))
    coriolis = -1e-4 .+ 1e-11 * reshape(yC, (length(yC), 1))
    κ =  -average1(vpbp) ./ by
    ν = coriolis .* average2(vpbp) ./ bz ./ uz
    
    #bz = sum(bz, dims = 1) ./ size(bz)[1]
    zCA = (zC[2:end] + zC[1:end-1])/2

    states = [eval(field) for field in fields]
    statenames = [string(field) for field in fields]
    units = ["[m/s]", "[m/s]", "[m/s]", "[m/s²]", "[m²/s³]", "[m²/s²]", "[m²/s³]", "[m²/s²]"]
    # bz
    push!(states, bz)
    push!(statenames, "∂ᶻb")
    push!(units, "[1/s²]")

    # by
    push!(states, by)
    push!(statenames, "∂ʸb")
    push!(units, "[1/s²]")

    # v'b'
    push!(states, vpbp)
    push!(statenames, "v'b'")
    push!(units, "[m²/s³]")
    # κ
    push!(states, κ)
    push!(statenames, "κ ≈  v'b' / ∂ʸb")
    push!(units, "[m²/s]")
    # ν
    push!(states, ν)
    push!(statenames, "ν ≈  f v'b' / ∂ᶻb / ∂ᶻu")
    push!(units, "[m²/s]")
    domain = [yF, zF]
    return states, statenames, units, domain
end