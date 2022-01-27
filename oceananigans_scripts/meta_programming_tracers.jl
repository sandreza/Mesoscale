# Step 1: Create Symbols for all the tracers
symbol_list = Symbol[]
for j in 1:2, k in 1:2
    push!(symbol_list, Meta.parse("c_k"*string(j) * "_j"*string(k)))
end
t_list = Tuple(symbol_list)

# Step 2: Create a named tuple list with the symbols # for now just set equal to one
nt_list = (; zip(t_list, (1,2,3,4))...)

# Step 3: Create a set of functions that do different things 
function_name_list = Symbol[]
for j in 1:2, k in 1:2
    push!(function_name_list, Meta.parse("forcing_c_k"*string(j) * "_j"*string(k)))
end
f_n_list = Tuple(function_name_list)
f_list = []
for i in eachindex(f_n_list)
    f_n = f_n_list[i]
    @eval begin
        function $f_n(x) 
            return x^2 + $i
        end
        push!($f_list, $f_n)
    end
end

# Step 4: Create a named tuple list with functions (don't need tuples before this step)
tracer_forcing_list = (; zip(t_list, Tuple(f_list))...)
# (; zip(symbol_list, f_list)...) # works


# Step 5: Automatically Construct Diagnostics
struct SimulationProxi{S}
    tracers::S
end

struct AveragedFieldProxi{S,T}
    field::S
    dims::T
end

nt_list = (; zip(t_list, (1,2,3,4))...)
field_list = (w=1, v = 2)
v = field_list.v
simulation = SimulationProxi(nt_list)
c_k1_j1 = simulation.tracers.c_k1_j1
averaged_outputs = Dict(
    :c_k1_j1 => AveragedFieldProxi(c_k1_j1, (1,)),
    :vc_k1_j1 => AveragedFieldProxi(v * c_k1_j1, (1,))
)
mean_outputs = Dict(
    :v => AveragedFieldProxi(v, (1,))
)
for i in eachindex(t_list)
    push!(mean_outputs, t_list[i] => AveragedFieldProxi(getproperty(simulation.tracers, t_list[i]), (1,))) 
end


