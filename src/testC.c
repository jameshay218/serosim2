#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* testing wow */
double testCFunc (double x){
	double y = x*x;
	return(y);
}
