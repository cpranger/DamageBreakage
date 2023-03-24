# DamageBreakage
Julia codes for simulating the Damage-Breakage rheology of Lyakhovsky and Ben-Zion in a continuum.

The main file is `scripts/dbr.jl`, but is currently in development and can't be executed.

All staggered grid operations are placed in the package StaggeredKernels.jl, to be found at `https://github.com/cpranger/StaggeredKernels.jl`.

Algorithms will be collected in the algorithms folder, which currently includes well-tested versions of
- power iterations,
- chebyshev iterations.

These algorithms are tested in
- `scripts/test-powerit.jl`,
- `scripts/test-chebyshev.jl`.
Other tests include
- `scripts/test-de-rham.jl`: check if the staggered grid conserves the vector calculus identities curl(grad(f)) = 0, div(curl(v)) = 0, grad(div(v)) - curl(curl(v)) - div(grad(v)) = 0, for any discrete scalar field f and staggered vector field v,
- `scripts/test-eigenmodes.jl`: check if the analytically computed scalar, divergence-free vectorial, and curl-free vectorial staggered eigenmodes satisfy the corresponding constraints and result in the correct eigenvalues.

any of these tests can be run directly from the shell using
```julia --project scripts/test-____.jl --nb 100,100 [--help for further arguments]```
or from the REPL using
```ARGS = ["--nb", "100,100"]; # add further arguments similarly```
```includet("./scripts/test-____.jl")```
```main()```