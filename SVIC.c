/* file sirmodess.c */
#include <R.h>
#include <math.h>
static double parms[10];
#define b1 parms[0]
#define b2 parms[1]
#define gamma parms[2]
#define rU parms[3]
#define rL parms[4]
#define rUn parms[5]
#define rLn parms[6]
#define frac_lower parms[7]
#define v parms[8]
#define prop parms[9]



double betafunc (double virulence, double rU_effect, double rL_effect, double frac_lower_effect) {
double result;
result = (1-frac_lower_effect)*b1*pow((1-rU_effect)*(virulence-.0025),b2) + (frac_lower_effect)*b1*pow((1-rL_effect)*(virulence-.0025),b2);
return result; 
}
       
/* initializer */
void initmod(void (* odeparms)(int *, double *))
{
  int N=10;
  odeparms(&N, parms);
}

/* ODEs */
void derivs (int *neq, double *t, double *y, double *ydot,double *yout, double *out, int *ip)
{

/* susceptible */
  
ydot[0] = -y[0] * (y[2]*betafunc(v,0,0,frac_lower) + y[3]*betafunc(v,rU,rL,frac_lower) + y[4]*betafunc(v,rUn,rLn,frac_lower));

/* vaccinated */
  
ydot[1] =  -y[1] * (1-rU) * (y[2]*betafunc(v,0,0,frac_lower) + y[3]*betafunc(v,rU,rL,frac_lower) + y[4]*betafunc(v,rUn,rLn,frac_lower));

/* infected--naive */
 
ydot[2] = y[0] * (y[2]*betafunc(v,0,0,frac_lower) + y[3]*betafunc(v,rU,rL,frac_lower) + y[4]*betafunc(v,rUn,rLn,frac_lower)) - (gamma + v*prop) * y[2];

/* infected--vaccinated */
  
ydot[3] = y[1] * (1-rU) * (y[2]*betafunc(v,0,0,frac_lower) + y[3]*betafunc(v,rU,rL,frac_lower) + y[4]*betafunc(v,rUn,rLn,frac_lower)) - (gamma + v*(1-rL)*prop) * y[3];

/* infected--natural immunity */

ydot[4] = y[5] * (1-rUn) * (y[2]*betafunc(v,0,0,frac_lower) + y[3]*betafunc(v,rU,rL,frac_lower) + y[4]*betafunc(v,rUn,rLn,frac_lower)) - (gamma + v*(1-rLn)*prop) * y[4];

/* convalescent */
  
ydot[5] = gamma * (y[2] + y[3] + y[4]) - y[5]* (1-rUn)*(y[2]*betafunc(v,0,0,frac_lower) + y[3]*betafunc(v,rU,rL,frac_lower) + y[4]*betafunc(v,rUn,rLn,frac_lower));


/* output for Fmat */
  
yout[0] = y[0]*betafunc(v,0,0,frac_lower);
yout[1] = y[0]*betafunc(v,rU,rL,frac_lower);
yout[2] = y[0]*betafunc(v,rUn,rLn,frac_lower);

yout[3] = y[1]*(1-rU)*betafunc(v,0,0,frac_lower);
yout[4] = y[1]*(1-rU)*betafunc(v,rU,rL,frac_lower);
yout[5] = y[1]*(1-rU)*betafunc(v,rUn,rLn,frac_lower);

yout[6] = y[5]*(1-rUn)*betafunc(v,0,0,frac_lower);
yout[7] = y[5]*(1-rUn)*betafunc(v,rU,rL,frac_lower);
yout[8] = y[5]*(1-rUn)*betafunc(v,rUn,rLn,frac_lower);
  
 /* output for Vmat */
  
yout[9] = gamma + v*prop;
yout[10] = 0;
yout[11] = 0;

yout[12] = 0;
yout[13] = gamma + v*prop*(1-rL);
yout[14] = 0;

yout[15] = 0;
yout[16] = 0;
yout[17] = gamma + v*prop*(1-rLn);
 
}

/* END file mymod.c */
