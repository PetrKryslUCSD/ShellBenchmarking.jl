# runtests.jl


for s in readdir("."; join=true, sort=true)
    f, e = splitext(s)
    if e == ".pdf" || e == ".vtu"
  try
   rm(s)
        catch
            @warn "Failed to remove $s"
        end
    end
end

for s in [
    "scordelis_lo_examples",
    "twisted_beam_examples",
    "LE5_Z_cantilever_examples",
    "clamped_hypar_examples",
    "cos_2t_press_hyperboloid_free_examples",
    "cos_2t_press_hyperboloid_fixed_examples",
    "raasch_examples",
    ]
    include(s * ".jl")
end
