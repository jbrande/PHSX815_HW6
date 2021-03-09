# PHSX815_HW6

This program integrates the function e^x on a closed interval, [0,1]. There are two choices of integration methods, Simpson's Rule and Gauss-Legendre Quadrature. These can be specified in the program arguments. When passing an integration method, you must also pass its associated parameters.

-Simpson				By default, uses Simpson's rule to integrate the function.

-Gauss 					Uses Gauss-Legendre quadrature to integrate the function.

-Nint					(integer) If Simpson's rule is used, this must be in the range [1,4]. If Gauss-Legendre is used, this must be in the range [1,5].



The program prints out the analytic solution, the numerical solution, and their difference. Then the program calculates some systematic comparisons of both methods and the analytic integral, and plots the results.

The program requires numpy and matplotlib.
