Base.@kwdef struct MeasuredParameter
    val::Float64
    stat::Float64
    model::Float64
    syst::Float64
end
"""
    MeasuredParameter(val, stat)

Create a `MeasuredParameter` with `val` and `stat` uncertainty. The model uncertainty and systematic uncertainty are set to zero.
"""
MeasuredParameter(val, stat) = MeasuredParameter(val, stat, 0.0, 0.0)

"""
    MeasuredParameter(str::AbstractString)

Parse a string of the form `val ± stat ± syst ± model`, or `val ± stat` into a `MeasuredParameter`.
"""
MeasuredParameter(str::AbstractString) = MeasuredParameter(Meta.parse.(split(str, "±"))...)
