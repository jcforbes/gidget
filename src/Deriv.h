// Take numerical derivative of a one-dimensional function

// Return maximum of a and b.
double maximum(const double a, const double b);

// Compute the numerical derivative of a given function:
//  f - 1d function which requires a set of auxilliary parameters.
//  x0- value of the independent variable at which to evaluate the derivative
// etol- fractional error tolerance
// *p - an object containing the parameters needed by f.
double derivDriver(double (*f)(double,void*),
		   const double x0,
		   const double etol,
		   void* p,double *err);

// Approximate f'(x0) by extrapolating h->0, where
// f'(x0) ~ (f(x0+h) - f(x0-h))/(2h)
double deriv(double (*f)(double,void*),const double x0,
	     const double h0,void *p,double *error);
