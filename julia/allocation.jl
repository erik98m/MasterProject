using JuMP
using HiGHS
using Allocations
using BenchmarkTools
using BenchmarkPlots
using StatsPlots
using Serialization
using Unitful

include("csvReader.jl")

mutable struct Bench_res 
    bench_mark_mine
    bench_mark_julia
end


set_initial_alloc(owners) = function(ctx)
    V = ctx.profile
    A = ctx.alloc_var
    for g in items(V)
        set_start_value(A[owners[g], g], 1)
    end
    return ctx
end

function my_alloc_mm(V, C=nothing; owners, cutoff=nothing, ignored_agents=[], solver=Allocations.conf.MIP_SOLVER, kwds...)
    Allocations.init_mip(V, solver; kwds...) |>
    
    set_initial_alloc(owners) |>
    Allocations.achieve_mm(cutoff) |>
    Allocations.enforce(C) |>
    Allocations.solve_mip |>
    Allocations.mm_result
end

convert_array(arr) = arr .+ 1

# index tells function to choose a given solution from a iteration
function get_best_solution(filename, index)
    matrix = process_csv_to_additive_matrices(filename)
    #for (i, matrix) in enumerate(matrix)
        #println("Additive Matrix for Iteration ", i - 1)
        #println(matrix)
    #end
    m = matrix[index]
    a = vec(m)
    return convert_array(a)
end

function get_matrix(filename, index)
    matrices = process_csv_to_additive_matrices(filename)
    matrix = matrices[index]
    return matrix
end

function benchmark_allocations(iterations::Int, owners_path::String, val_profile_path::String)
    n =  []
    for i in 0:iterations-1 
        solution = get_best_solution(owners_path, i+1)        
        val_profile = get_matrix(val_profile_path, i+1)
        b1 = @benchmark my_alloc_mm($val_profile, owners = $solution) samples = 10 evals = 1 seconds = Inf
        b2 = @benchmark alloc_mm($val_profile) samples = 10 evals = 1 seconds = Inf
        res = Bench_res(b1,b2)
        push!(n, res)
    end

    return n
end


function get_median(l)
    return [median(x.bench_mark_mine.times) for x in l] ./ [median(x.bench_mark_julia.times) for x in l]
end

function get_mean(l)
    return [mean(x.bench_mark_mine.times) for x in l] ./ [mean(x.bench_mark_julia.times) for x in l]
end


n2 = deserialize("Benchmark_bin/eqValprofile/2agents.bin")
#=
n3 = deserialize("Benchmark_bin/eqValprofile/3agents.bin")
n4 = deserialize("Benchmark_bin/eqValprofile/4agents.bin")
n5 = deserialize("Benchmark_bin/eqValprofile/5agents.bin")
n6 = deserialize("Benchmark_bin/eqValprofile/6agents.bin")
n7 = deserialize("Benchmark_bin/eqValprofile/7agents.bin")
n8 = deserialize("Benchmark_bin/eqValprofile/8agents.bin")
n9 = deserialize("Benchmark_bin/eqValprofile/9agents.bin")
n10 = deserialize("Benchmark_bin/eqValprofile/10agents.bin")

a2 = get_mean(n2)
a3 = get_mean(n3)
a4 = get_mean(n4)
a5 = get_mean(n5)
a6 = get_mean(n6)
a7 = get_mean(n7)
a8 = get_mean(n8)
a9 = get_mean(n9)
a10 = get_mean(n10)=#

#println(mean(a2))




#b = @benchmark my_alloc_mm($val, owners = $x) samples = 100 evals = 1 seconds = Inf

#display(b)
#b2 = @benchmark alloc_mm($val) samples = 100 evals = 1 seconds = Inf
#display(b2)



b = n2[1]
s = b.bench_mark_mine
t = b.bench_mark_julia

display(s)
display(t)


