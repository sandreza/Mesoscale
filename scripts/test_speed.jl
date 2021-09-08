using BenchmarkTools

add_two_things(x,y) = x+y

add_two_things_2(x,p) = x+p.y

@btime add_two_things(10,10)

parameters = (y = 10,)
z = 10
const conts_parameters = (y=10,)

@btime add_two_things_2(10,conts_parameters)



