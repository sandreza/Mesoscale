function get_data(file)
    data_dictionary = Dict()
    for key in keys(file["timeseries"])
        println(key)
        tmp = []
        for t in keys(file["timeseries"]["t"])
            if key != "t"
                push!(tmp, file["timeseries"][key][t][1, :,:])
            else
                push!(tmp, file["timeseries"][key][t])
            end
        end
        data_dictionary[key] = tmp
    end
    return data_dictionary
end

filename = "/Hybrid_20_zonal_averages.jld2"
file = jldopen( pwd() * filename)
z = file["grid"]["zC"][2:end-1]
y = file["grid"]["yC"][2:end-1]
yF = file["grid"]["yF"][2:end-1]
zF = file["grid"]["zF"][2:end-1]
data_dictionary = get_data(file)
close(file)
##
function get_y(a, y, yF)
    length(a[:,1]) == length(y) ? y : yF
end
function get_z(a, z, zF)
    length(a[1,:]) == length(z) ? z : zF
end

function plot_field(a)
    py = get_y(a, y, yF) # global scope
    pz = get_z(a, z, zF) # global scope
    clims = (minimum(a), maximum(a))
    p1 = contourf(py, pz, a', 
    color = :thermometer, title = "Zonal and Temporal Average ",
    xlabel = "Meridional [m]", ylabel = "Depth [m]"
    , clims = clims, linewidth = 0)
    display(p1)
end

#plot_field(data_dictionary["b"][1])
a = mean(data_dictionary["b"][2:end])
plot_field(a)








