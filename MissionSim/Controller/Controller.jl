module Controller

include("controller_config_file.jl")
# include("../system_config_file.jl")
using ..CustomStructs

export generate_command


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
#                        DETUMBLER                                 #
####################################################################

####################################################################
#                       MAGNETORQUER                               #
####################################################################
end