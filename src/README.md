## Some details on the code

The script "./model_CCE_CCM.jl" contains all the functions that are
called in the examples to run CCE, CCM, as well as the IBVP and CIBVP
independently. These functions are "run_cce", "run_ccm", and
"run_ibvp_cibvp", respectively. The code could be optimized to avoid
repetitions in building these functions, but this is left for future
work.

All the derivatives that appear on the right-hand-side (RHS) of any
equation are approximated with 2nd order accurate, centered, finite
difference operators. These operators can be found in the script
"./operators.jl". The script containts upwind operators as well (for
left- and right-moving variables along the angular direction, for
experimentations). The latter are not called in the module, but
whoever is interested can make the appropraite changes to test
them. In the noisy norm convergence tests, they have resulted in less
clean 2nd order convergence.

The null integration of the intrinsic equations is performed using the
[trapezoidal
rule](https://en.wikipedia.org/wiki/Trapezoidal_rule_(differential_equations)). The
integration in time is performed with the [4th order Runge-Kutta
method](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods),
for both the IBVP and the CIBVP. There is no artificial dissipation.

The Cauchy and characteristic domains touch on their right and left
boundary, respectively. When they are evolved indepentently, boundary
data are prescribed for each on of them. For CCE, the solution of the
right-moving fields of the IBVP serves as boundary data for teh
right-moving fields of the characteristic domain. For CCM, in addition
to the latter, the solution of the left-moving characteristic field
serves as boundary data for the left-moving Cauchy field. The
implementation of the boundary data is the difference between the
functions ""run_ibvp_cibvp", "run_cce", and "run_ccm".

