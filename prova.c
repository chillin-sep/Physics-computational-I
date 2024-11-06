#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "pips.h"


int main(){


 i_t ci ={.x0=5,.beta=0.4,.v0=15,.T=100,.w=3};

 double dt = 0.001;
 double t = 0.;
 ph_t ph;

 ph.x = ci.x0;
 ph.v = ci.v0;

 while(t<ci.T)
 {
   printf("%lf\t%lf\t%lf\n", t, ph.x, ph.v);
   ph=Eulero(ph,acc_osc(ph,ci,t),dt);
   t += dt;
 }


}