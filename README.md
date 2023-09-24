# Continuous-time auxiliary-field Quantum Monte Carlo method for Anderson impurity model with Bethe-lattice-bath electrons

See, E. Gull et al., EPL 82, 57003 (2008)

Yuki Nagai, Ph.D 10/24/2017(MM/DD/YY)

This code was made by Julia 0.6.0.

# Note: This code with Julia 0.6.0 is 6 times slower than that with Fortran90.

If you know how to speedup this code, please tell me.

The code was revised for Julia 0.7. 
# The code in Julia 0.7 is 5 times slower than that with Fortran90. 

# update (2023/09/24)
I wrote new code for Julia 1.9.3. The code is ctaux_faster.jl

Old code in Julia 1.9.3 (ctaux_juliav1.9.jl):
```
 15.258614 seconds (170.55 M allocations: 10.790 GiB, 1.24% gc time, 0.11% compilation time)
```
New code in Julia 1.9.3 (ctaux_faster.jl)
```
 10.448534 seconds (2.05 k allocations: 123.203 KiB)
```