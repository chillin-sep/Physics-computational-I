#include <stdio.h>
#include <math.h>
#include <stdlib.h>

typedef struct phases {double x,v;}ph_t;

typedef struct initial_condition { double x0, v0, T, k, beta, m, w;}i_t;


//-------FUNZIONI PER LE ACCELERAZIONI--------//

double acc_osc(ph_t ph,i_t ci,double t)
{
  double a = -pow(ci.w,2)*ph.x-ci.beta * ph.v;

  return a;
}


//--------METODI DI INTEGRAZIONE--------//

ph_t Eulero(ph_t ph,double a,double dt)
{
    ph_t pips;
    pips.x = ph.x + ph.v * dt;
    pips.v = ph.v + a * dt;
    
    return pips;
}