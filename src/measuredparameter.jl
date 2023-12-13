Base.@kwdef struct MeasuredParameter
    val::Float64
    stat::Float64
    model::Float64 = 0.0
    syst::Float64 = 0.0
end
MeasuredParameter(str::AbstractString) =
    MeasuredParameter(Meta.parse.(split(str, "±"))...)
