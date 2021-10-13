# ASCC2022
Consider a nonlinear autonomous system as follows,
$$\dot{x}=f(x)+g(x)+d(x),$$
where $f(x)$ denotes the polynomial term, $g(x)$ denotes the non-polynomial term and $d(x)\sim (0,\sigma^2_n)$ denotes the noise over this system. 

In this repo, we use
- Chebfun Toolbox: To approximate nonlinear terms by Chebyshev Interpolants,
- GPML Toolbox: Expressed the unknown term $d(x)$ above into the polynomial form,
- SOSOPT/SOSTOOL/Mosek: To solve some sum-of-squares programmings in this learned polynomial system.

Note that, please run sosaddpath.m and toolbox/gpml/startup.m at the beginning. 
