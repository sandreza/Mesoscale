# printing
using Printf, Statistics, JLD2
# Oceananigans
using Oceananigans
using Oceananigans.BoundaryConditions
using Oceananigans.Fields
using Oceananigans.OutputWriters
using Oceananigans.Diagnostics
using Oceananigans.Utils
using Oceananigans.AbstractOperations
using Oceananigans.Advection
# CUDA
using CUDA

include(pwd() * "/scripts/ic.jl")