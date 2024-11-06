#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define g 9.80352  
#define PI 3.14159265358979323846264338327950288 //41971693993751058209749445923078164062862089986280348253421170679 

typedef struct phasespace
{
    double x, v;
}ps_t;

typedef struct initial_condition
{
    double x0, v0, T, k, m, l, beta, F, wF, w;
}ic_t;

typedef double (*acc_t)(ps_t, ic_t, double);




//--------------- FUNZIONI PER ACCELERAZIONE-----------------


double a_osc( ps_t d1, ic_t ci, double t)
{
    double a = - ( d1.x * ( ci.k / ci.m ) + d1.v * ci.beta);
    return a;
}

double a_pend(ps_t d1, ic_t ci, double t)
{
    double b = - pow(ci.w,2) * sin(d1.x) - ci.beta * d1.v;
    return b;
}

double a_forz(ps_t d1, ic_t ci, double t)
{   
    double c = 0;
    c = (ci.F * cos(ci.wF * t)) - (pow(ci.w,2) * sin(d1.x) + ci.beta * d1.v) ;
    return c; 
}


//-----------METODI DI INTEGRAZIONE-------------------

ps_t integration_euler( ps_t  d1, double a, double dt)
{
    ps_t sf;
    sf.x =  d1.x + d1.v * dt;
    sf.v = d1.v + a * dt;

    return sf;
}

ps_t integration_euler_cromer( ps_t d1, double a , double dt)
{
    ps_t sf;
    sf.v =  d1.v + a * dt;
    sf.x = d1.v + a * dt;
    
    return sf;
}

ps_t integration_velocity_verlet(ps_t d1, ic_t ci, acc_t acc, double t,  double dt)
{
    ps_t sf;
   
    
    sf.x = d1.x + d1.v * dt + 0.5 * acc(d1, ci, t) * pow(dt,2);
    

    sf.v = d1.v + 0.5 * (acc(d1, ci, t) + acc(sf, ci, t)) * dt;

    return sf;
}

ps_t integration_verlet(ps_t d1, double a, double t, double dt)
{
    ps_t sf;


    return sf;
}

ps_t integration_rungekutta2(ps_t d1, ic_t ci, acc_t a, double t, double dt)
 {
    ps_t sf;
    ps_t k1, k2;
    

    
    
    
    k1.x = d1.x + d1.v * dt;
    k1.v = d1.v + a(d1, ci, t) * dt;
    k2.x = dt * ( d1.v + k1.v / 2 ) / 2;
    k2.v = dt * a(k2, ci, t );


    sf.v = d1.v + k2.v;
    sf.x = d1.x + sf.v * dt;

    return sf;
}

ps_t integration_rungekutta4(ps_t d1, ic_t ci, acc_t a, double t, double dt) {
    ps_t k1, k2, k3, k4, ktemp, sf;

    // Calcolo di k1
    k1.x = d1.v * dt;
    k1.v = a(d1, ci, t) * dt;

    // Calcolo di k2
    ktemp.x = d1.x + 0.5 * k1.x;
    ktemp.v = d1.v + 0.5 * k1.v;
    k2.x = ktemp.v * dt;
    k2.v = a(ktemp, ci, t + dt / 2) * dt;

    // Calcolo di k3
    ktemp.x = d1.x + 0.5 * k2.x;
    ktemp.v = d1.v + 0.5 * k2.v;
    k3.x = ktemp.v * dt;
    k3.v = a(ktemp, ci, t + dt / 2) * dt;

    // Calcolo di k4
    ktemp.x = d1.x + k3.x;
    ktemp.v = d1.v + k3.v;
    k4.x = ktemp.v * dt;
    k4.v = a(ktemp, ci, t + dt) * dt;

    // Media pesata dei coefficienti per ottenere il risultato finale
    sf.x = d1.x + (k1.x + 2 * k2.x + 2 * k3.x + k4.x) / 6.0;
    sf.v = d1.v + (k1.v + 2 * k2.v + 2 * k3.v + k4.v) / 6.0;

    return sf;
}



//------------FUNZIONE PER CALCOLARE IL PERIODO-------------

double period(FILE *p, ic_t ci)
{
    double t = 0, x, v, x0, t0 = 0;
    int count = 0;
    double period = 0;

    fscanf(p,"%lf %lf %lf\n", &t, &x, &v);

    while( t < (ci.T - 1) )
    {
        x0 = x;

        fscanf(p,"%lf %lf %lf\n", &t, &x, &v);

        if( x0 * x < 0 )
        {
            count += 1;
            if ( count != 1 ) period += (t - t0 );
            t0 = t;
        }

    }

    period = 2 * period / (count - 1 );

    return period;
}


//-------------FUNZIONE PER L'ENERGIA IN BASE AL PASSO DI INTEGRAZIONE--------

void energy_vs_dt(FILE *nrg)
{
    
}

//-------------FUNZIONE SEGNO----------
int sgn(double x)
{
   int sgn = fabs(x) / x;
    return sgn;
}

//--------FUNZIONE PER FAR RIENTRARE L'ANGOLO IN [-PI,PI]------------- 
double fix_theta(double x)
{
    if( x < -PI || x > PI ) x -= 2 * PI * sgn(x);
    return x;
}
