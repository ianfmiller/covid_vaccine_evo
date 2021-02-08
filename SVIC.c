/* file sirmodess.c */
#include <R.h>
#include <math.h>
static double parms[11];
#define b1 parms[0]
#define b2 parms[1]
#define gamma parms[2]
#define sigma parms[3]
#define rUv parms[4]
#define rLv parms[5]
#define rUc parms[6]
#define rLc parms[7]
#define frac_lower parms[8]
#define v parms[9]
#define x parms[10]


double betafunc (double virulence, double rU_effect, double rL_effect, double frac_lower_effect) {
double result;
result = (1-frac_lower_effect)*b1*pow((1-rU_effect)*virulence,b2) + (frac_lower_effect)*b1*pow((1-rL_effect)*virulence,b2);
return result; 
}
       
/* initializer */
void initmod(void (* odeparms)(int *, double *))
{
  int N=9;
  odeparms(&N, parms);
}

/* ODEs */
void derivs (int *neq, double *t, double *y, double *ydot,double *yout, double *out, int *ip)
{

double foi; 
foi = betafunc(v,0,0,frac_lower)*(0.5*y[2]+y[5]) + betafunc(v,rUv,rLv,frac_lower)*(0.5*y[3]+y[6]) + betafunc(v,rUc,rLc,frac_lower)*(0.5*y[4]+y[7]);

/* susceptible */
  
ydot[0] = -y[0] * foi;

/* vaccinated */
  
ydot[1] =  -y[1] * (1-rUv) * foi;

/* infected,pre-symptomatic--naive */

ydot[2] = y[0] * foi - y[2]*sigma;

/* infected,pre-symptomatic--vaccinated */

ydot[3] = y[1] * (1-rUv) * foi - y[3]*sigma;

/* infected,pre-symptomatic--natural immunity */

ydot[4] = y[8] * (1-rUc) * foi - y[4]*sigma;

/* infected,symptomatic--naive */
 
ydot[5] = y[2]*sigma - y[5] * (gamma + x*v);

/* infected,symptomatic--vaccinated */
  
ydot[6] = y[3]*sigma - y[6] * (gamma + (1-rLv)*x*v);

/* infected,symptomatic--natural immunity */

ydot[7] = y[4]*sigma - y[7] * (gamma + (1-rLc)*x*v);

/* convalescent */
  
ydot[8] = gamma * (y[5] + y[6] + y[7]) - y[8]* (1-rUc) * foi ;


/* output for Fmat */

yout[0] = y[0]*0.5*betafunc(v,0,0,frac_lower);
yout[1] = y[0]*0.5*betafunc(v,rUv,rLv,frac_lower);
yout[2] = y[0]*0.5*betafunc(v,rUc,rLc,frac_lower);
yout[3] = y[0]*betafunc(v,0,0,frac_lower);
yout[4] = y[0]*betafunc(v,rUv,rLv,frac_lower);
yout[5] = y[0]*betafunc(v,rUc,rLc,frac_lower);

yout[6] = y[1]*(1-rUv)*0.5*betafunc(v,0,0,frac_lower);
yout[7] = y[1]*(1-rUv)*0.5*betafunc(v,rUv,rLv,frac_lower);
yout[8] = y[1]*(1-rUv)*0.5*betafunc(v,rUc,rLc,frac_lower);
yout[9] = y[1]*(1-rUv)*betafunc(v,0,0,frac_lower);
yout[10] = y[1]*(1-rUv)*betafunc(v,rUv,rLv,frac_lower);
yout[11] = y[1]*(1-rUv)*betafunc(v,rUc,rLc,frac_lower);

yout[12] = y[8]*(1-rUc)*0.5*betafunc(v,0,0,frac_lower);
yout[13] = y[8]*(1-rUc)*0.5*betafunc(v,rUv,rLv,frac_lower);
yout[14] = y[8]*(1-rUc)*0.5*betafunc(v,rUc,rLc,frac_lower);
yout[15] = y[8]*(1-rUc)*betafunc(v,0,0,frac_lower);
yout[16] = y[8]*(1-rUc)*betafunc(v,rUv,rLv,frac_lower);
yout[17] = y[8]*(1-rUc)*betafunc(v,rUc,rLc,frac_lower);

yout[18] = 0;
yout[19] = 0;
yout[20] = 0;
yout[21] = 0;
yout[22] = 0;
yout[23] = 0;

yout[24] = 0;
yout[25] = 0;
yout[26] = 0;
yout[27] = 0;
yout[28] = 0;
yout[29] = 0;

yout[30] = 0;
yout[31] = 0;
yout[32] = 0;
yout[33] = 0;
yout[34] = 0;
yout[35] = 0;

 /* output for Vmat */
  
yout[36] = sigma;
yout[37] = 0;
yout[38] = 0;
yout[39] = 0;
yout[40] = 0;
yout[41] = 0;

yout[42] = 0;
yout[43] = sigma;
yout[44] = 0;
yout[45] = 0;
yout[46] = 0;
yout[47] = 0;

yout[48] = 0;
yout[49] = 0;
yout[50] = sigma;
yout[51] = 0;
yout[52] = 0;
yout[53] = 0;

yout[54] = -1*sigma;
yout[55] = 0;
yout[56] = 0;
yout[57] = gamma + x*v;
yout[58] = 0;
yout[59] = 0;

yout[60] = 0;
yout[61] = -1*sigma;
yout[62] = 0;
yout[63] = 0;
yout[64] = gamma+(1-rLv)*x*v;
yout[65] = 0;

yout[66] = 0;
yout[67] = 0;
yout[68] = -1*sigma;
yout[69] = 0;
yout[70] = 0;
yout[71] = gamma + (1-rLc)*x*v;
 
}

/* END file mymod.c */
