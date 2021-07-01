module Controller

include("controller_config_file.jl")
# include("../system_config_file.jl")
using ..CustomStructs
using LinearAlgebra 

include("detumbler.jl")


export generate_command
export DETUMBLER


# DIMENSIONS OF u?


####################################################################
#                         TRIVIAL                                  #
####################################################################
# struct TRIVIAL_CONTROLLER
#     junk
# end

function generate_command(cont::TRIVIAL, state)
    return zeros(3), true
end


####################################################################
#                       MAGNETORQUER                               #
####################################################################



end