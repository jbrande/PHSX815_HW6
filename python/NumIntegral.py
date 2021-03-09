import numpy as np
import sys
import matplotlib.pyplot as plt

sys.path.append(".")

# Integrate under the curve e^x from 0 to 1 with either of two integration methods:
# Simpson's Rule or Gauss-Legendre quadrature with varying numbers of sub-intervals (default 2)
# the analytic solution will be e - 1

if __name__ == "__main__":

	#set default number of steps
	Nsteps = 10

	# set default number of sub-intervals
	Nint = 2

	method = 0

	# read the user-provided seed from the command line (if there)
	if '-Nsteps' in sys.argv:
		p = sys.argv.index('-Nsteps')
		Nsteps = int(sys.argv[p+1])
	if '-h' in sys.argv or '--help' in sys.argv:
		print ("Usage: %s -Nsample [number]" % sys.argv[0])
		sys.exit(1)
	
	if '-Simpson' in sys.argv:
		method = 0
	elif '-Gauss' in sys.argv:
		method = 1
	else:
		print("You must specify a method of integration!")
		sys.exit(1)
	
	# specify number of sub-intervals
	if '-Nint' in sys.argv:
		p = sys.argv.index("-Nint")
		Nint = int(sys.argv[p+1])
		if method == 0 and (Nint < 1 or Nint > 4):
			print("For -Simpson, Nint must be in the range [1,4]")
			sys.exit(1)
		if method == 1 and (Nint < 1 or Nint > 5):
			print("For -Gauss, Nint must be in the range [1,5]")


	# function pointer to integrate over, e^x.
	fn = np.exp

	# coefficients for the integration methods
	coeff_NC = [[1.0, 1.0], 
				[1.0, 4.0, 1.0],
				[1.0, 3.0, 3.0, 1.0],
				[7.0, 32.0, 12.0, 32.0, 7.0]]

	# weights, roots taken from Wiki on Gaussian quadrature - https://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss%E2%80%93Legendre_quadrature
	weights_gauss = [[2.0],
					 [1.0, 1.0],
					 [0.555556, 0.888889, 0.555556],
					 [0.347855, 0.652145, 0.652145, 0.347855],
					 [0.236927, 0.478629, 0.568889, 0.478629, 0.236927]]

	roots_gauss = [[0.0],
				   [-0.57735, 0.57735],
				   [-0.774597, 0.0, 0.774597],
				   [-0.861136, -0.339981, 0.339981, 0.861136],
				   [-0.90618, -0.538469, 0.0, 0.538469, 0.90618]]

	# simpson's rule
	def simpson(a, b, n):
		h = (b-a)/n
		prefactor = [h/2.0, h/3.0, 3.0*h/8.0, 2.0*h/45.0][n-1]

		s_int = 0.0
		#print(a,b)
		for i in range(0,n+1):
			xi = a + i*h
			#print("xi, ", xi)
			s_int += prefactor*coeff_NC[n-1][i]*fn(xi)
		return s_int

	# gauss legendre, need to do change of interval - see https://en.wikipedia.org/wiki/Gaussian_quadrature#Change_of_interval
	def gauss(a, b, n):
		# (b-a / 2) * sum(w_i * f((b-a / 2)*xi + (b+a)/2))
		g_int = 0.0
		for i in range(1, n+1):
			g_int += ((b-a)/2.0) * weights_gauss[n-1][i-1] * fn(((b-a)/2.0)*roots_gauss[n-1][i-1] + ((b+a)/2))
		return g_int

	# get intervals on [0, 1]
	intervals = []
	rng = np.linspace(0, 1, Nsteps+1)
	for i in range(Nsteps):
		intervals.append((rng[i], rng[i+1]))

	# calculate integral
	integral = 0.0

	if method == 0:
		for interval in intervals:
			integral += simpson(interval[0], interval[1], Nint)
	elif method == 1:
		for interval in intervals:
			integral += gauss(interval[0], interval[1], Nint)
	else:
		print("Something broke!")
		sys.exit(1)


	analytic = np.e - 1.0

	print("Integrating e^x between 0 and 1.")
	print("Analytic answer: ", analytic)
	print("Numerical answer: ", integral)
	print("Analytic-Numerical: ", analytic-integral)



	# code to do some systematic comparisons
	
	simp = np.zeros(4)
	gaus = np.zeros(5)

	# get simpson as function of Nint
	for Nint in [1,2,3,4]:
		for interval in intervals:
			simp[Nint-1] += simpson(interval[0], interval[1], Nint)

	# get gauss as function of Nint
	for Nint in [1,2,3,4,5]:
		for interval in intervals:
			gaus[Nint-1] += gauss(interval[0], interval[1], Nint)


	# compare numerical methods to analytic
	fig = plt.figure(figsize=(10, 6))
	plt.plot([1,2,3,4], simp, ".-", label="Simpson's rule")
	plt.plot([1,2,3,4,5], gaus, ".-", label="Gauss-Legendre")
	plt.axhline(analytic, linestyle="--", color="k", label="Analytic Solution")
	plt.legend()
	plt.xlabel("Number of Sub-Intervals")
	plt.ylabel(r'$\int e^{x} dx$')
	plt.title("Analytic vs Numerical Integration")
	plt.show()
	fig.savefig("analytic_numerical_integration.jpg", dpi=180)

	# compare numerical methods to each other
	fig = plt.figure(figsize=(10, 6))
	plt.plot([1,2,3,4], 1 + (simp-gaus[:4]), ".-", label="Simpson - Gauss")
	plt.plot([1,2,3,4], 1 + (simp-analytic), ".-", label="Simpson - Analytic")
	plt.plot([1,2,3,4, 5], 1 + (gaus-analytic), ".-", label="Gauss - Analytic")
	plt.legend()
	plt.yscale("log")
	plt.xlabel("Number of Sub-Intervals")
	plt.ylabel("Error")
	plt.title("Error Comparison of the Methods")
	plt.show()
	fig.savefig("comparison.jpg", dpi=180)

	# DISCUSSION:
	# For n = 1, the S-G difference is not a good tracer of the actual errors S-A and G-A. For n = 2 through 4,5, the estimates all converge.