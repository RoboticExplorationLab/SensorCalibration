module Controller

using ..CustomStructs
using LinearAlgebra 
using PyCall, Infiltrator

include("detumbler.jl")

export generate_command
export DETUMBLER
export __init__

####################################################################
#                         TRIVIAL                                  #
####################################################################

function generate_command(cont::TRIVIAL)
    return zeros(3)
end


####################################################################
#                       MAGNETORQUER                               #
####################################################################

# Eventually add an MPC for orientation control


end