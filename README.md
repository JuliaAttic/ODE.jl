# Alternate iterators for ODE.jl

This is a prototype for an implementation of iterators for ODE.jl
slightly different to pwl's PR49.  It is a partial rewrite of pwl's work and
tries to use 20-20 hindsight to come up with a cleaner
implementation.  The idea is to use bits of this to update PR49.

To dive in: run the tests with `include("test.jl")`, and look at the
examples in the files `ex_*.jl`.

### To read the code

- `ODE.jl` just contains the `include`s.
- `core.jl` is the core, I suggest starting there. The core links up
  the internals (solvers+helpers) with the user interface and
  separates the two completely.

### Internals

- `runge_kutta.jl` is an implementation of a fixed-step RK4 method
- `rosenbrock.jl` is the `ode23s` solver ported.  Note that this is
   not using in-place functions yet!
- `fixed_adam.jl` fixed step Adam-B and Adam-M methods (by @obiajulu)
- `dense.jl` dense output wrapper
- `events.jl` does event handling (root finding)

User interface in `ui.jl`.  There is a Julian interface (both
high-level and mid-level) and the old Matlab-style interface.

Tests and examples are in `test.jl` and uses the Julia and Matlab
interface.  There is an example on the iteration interface in
`ex_iter.jl`.

### Performance

The `test.jl` does timings.  In particular the last loop tests the
high-level interface and compares it to ODE.jl.  For instance this
`ode23s` seems 10x faster on the most complicated test.  All in all
performance seems good but needs testing with IVPTestSuite.jl.  One
test of `IVPTestSuite.jl` is in `ex_bruss.jl` and shows that old and
new ode23s are equal.

### ToDo

- [X] UI (to some extent)
- [ ] RK solvers
- [X] ode23s (except in-place functions)
- [ ] refactor dense output to allow for other schemes.
- [ ] add a IVPTestSuite.jl wrapper
- [X] events (needs refactor)
- [ ] return status in iteration
- [ ] wrap a Sundials solver to see how that would work

### To ponder

- [X] where to put Jacobian: into IVP.  --> All math-stuff goes
      into IVP.
- [ ] Think about a modular approach to options.  Maybe options should
      be mostly held in a `Dict{Symbol,Any}` and only later be
      converted into a specific option datatype.  -> Yes, so each
      solver holds its options which can be constructed from the
      options dictionary.  However, I'm not sure yet how to handle
      options which are not meant for the solver but for, say, dense
      output. --> `Problem` holds the option-dict and each part
      (solver, dense, events, etc) parses opts and takes what it
      needs.

### Benchmarks
- [Julia 0.4.6](https://gist.github.com/mauro3/c26836f20f750621aeeb3d03c95cc186#file-benchmarks4-md)
- [Julia 0.5](https://gist.github.com/mauro3/c26836f20f750621aeeb3d03c95cc186#file-benchmarks5-md)
