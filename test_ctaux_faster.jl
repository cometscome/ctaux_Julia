include("./ctaux_faster.jl")
using .CTAUX_faster

function main_faster()

    β=10.0
    U = 2.0
    μ=U/2
    K = 1.0
    mqs = 1000000
    ntime = 1024 #number of τs
    mfreq = 1024 #number of ωns
    norbs = 2 #number of orbitals. norbs = 2 for 1-band model.
    V = 1.0 #Strength of the hybridization
    nthermal = 1000
    mkink = 1024
    
    ctaux_main(U,K,β,V,μ;ntime,mfreq,mkink,nthermal,norb,mqs)
end
main_faster()
