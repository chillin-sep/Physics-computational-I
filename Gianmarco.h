#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define TOLERANCE 1e-5


char sgn(double x);  // funzione segno
double abso(double x); // funzione valore assoluto

typedef struct d1{
    double x, v, a, E;
}d1_t;

typedef struct d2{
    double theta, v, a, E;
}d2_t;

typedef struct k{
    double x, v;
}k_t;

typedef struct ktheta{
    double theta, v;
}kt_t;

typedef struct f{
    double fzero, omega;
}f_t;

double aosc(double w, double x){
   return - w * w * x;  
}

double apen(double w, double theta, double beta, double v, struct f *f, double T, double m){
    return - w * w * sin(theta) - (v * beta)/m + (f->fzero * cos(f->omega * T))/m;
}


double fix(double theta){
    if(theta > M_PI){
        theta -= 2 * M_PI;
    }

    if(theta < - M_PI){
        theta += 2 * M_PI;
    }
    return theta;
}


void energia(double m, double k, double x, double v, double Tmax);

char sgn(double x){
    char segno='+';
     
    if(x<0){
     segno='-';
    }
 return segno;
}


double abso(double x){
    
    if(x<0){
        x=-x;
    }

return x;
}


void VelVerlet(struct d1 *d1, double m, double k, double dt, double Tmax){
    FILE *fp;
    double T=0.0;
    double anew; 
    double E;

    fp = fopen("VerletVel.dat","w");
    while(T<Tmax){
        d1->a = -k * d1->x / m; // Aggiorna la posizione usando la posizione corrente
        d1->x = d1->x + d1->v * dt + 0.5 * d1->a * dt * dt;  // Aggiorna la posizione usando la velocità correnteù
        anew = -k * d1->x / m; // Aggiorna la accelerazione con la nuova posizione
        d1->v = d1->v + 0.5 * (d1->a + anew) * dt;  // Aggiorna la velocità usando l'accelerazione corrente
        d1->E = 0.5 * m * pow(d1->v, 2) + 0.5 * k * pow(d1->x, 2); // Calcola l'energia totale del sistema
        T += dt;

        fprintf(fp,"%lf %lf %lf %lf %lf\n", T, d1->x, d1->v, d1->a, d1->E);
        
        

    }
    fclose(fp);
}


void energia(double m, double k, double x, double v, double Tmax){
    FILE *fp;
    double T=0.0;
    double anew; 
    double E, Ein, deltaE;
    double a;
    double dt=0.001;
    double xin, vin;

    x=x;
    v=v;
    xin=x;
    vin=v;
    Ein = 0.5 * m * pow(v, 2) + 0.5 * k * pow(x, 2);

    fp = fopen("deltaE.dat","w");

    while(dt<1){
    while(T<Tmax){
        a = -k * x / m; // Aggiorna la posizione usando la posizione corrente
        x = x + v * dt + 0.5 * a * dt * dt;  // Aggiorna la posizione usando la velocità correnteù
        anew = -k * x / m; // Aggiorna la accelerazione con la nuova posizione
        v = v + 0.5 * (a + anew) * dt;  // Aggiorna la velocità usando l'accelerazione corrente
        E = 0.5 * m * pow(v, 2) + 0.5 * k * pow(x, 2); // Calcola l'energia totale del sistema
        T += dt;
     }

    deltaE = E-Ein;
    T=0.0;
    x=xin;
    v=vin;
    fprintf(fp,"%lf %lf %lf \n", dt, deltaE, abso(deltaE));
    dt += 0.001;
    }
    fclose(fp);
}

void VelVerletTheta(struct d2 *d2, double L, double m, double theta, double v, double dt, double Tmax){
    FILE *fp;
    double T=0.0;
    double anew; 
    double E, g = 9.81;

    d2->theta=theta;
    d2->v=v;
    fp = fopen("VerletVelTheta.dat","w");
    while(T<Tmax){
        d2->a = -g * sin(d2->theta) / L ; // Aggiorna la posizione usando la posizione corrente
        d2->theta = d2->theta + d2->v * dt + 0.5 * d2->a * dt * dt;  // Aggiorna la posizione usando la velocità correnteù
        anew = -g * sin(d2->theta) / L; // Aggiorna la accelerazione con la nuova posizione
        d2->v = d2->v + 0.5 * (d2->a + anew) * dt;  // Aggiorna la velocità usando l'accelerazione corrente
        d2->E = 0.5 * m * pow(d2->v, 2) * L + m * g * L * (1 - cos(d2->theta)); // Calcola l'energia totale del sistema
        T += dt;

        fprintf(fp,"%lf %lf %lf %lf %lf\n", T, d2->theta, d2->v, d2->a, d2->E);
        
        

    }
    fclose(fp);
}





void periodo(double L, double g, double m, double v, double dt, double thetamax){
    FILE *fp;
    double T=0.0, anew, thetatemp, E, thetai = 0.01, a,Ttheta[1000000], Tmax=10000.0, dtheta=0.01, Ttot = 0.0, Tcalcolato, Tanalitico, theta;
    int npasso, i = 0, z, x;
    fp = fopen("periodipendolo.dat","w");

    Tanalitico = 2 * M_PI * sqrt(L / g);
    
    while (thetai < thetamax) {
        theta = thetai;
        T = 0.0;
        i = 0;
        Ttot = 0.0;  // Inizializza la somma del periodo
        v=0.0;
        
        while (T < Tmax) {
        
            a = -g * sin(theta) / L;  // Aggiorna la posizione usando la posizione corrente
            thetatemp = theta;
            theta = theta + v * dt + 0.5 * a * dt * dt;  // Aggiorna la posizione
            anew = -g * sin(theta) / L;  // Aggiorna l'accelerazione
            v = v + 0.5 * (a + anew) * dt;  // Aggiorna la velocità
            
            E = 0.5 * m * pow(v, 2) * L + m * g * L * (1 - cos(theta));  // Calcola l'energia totale

            // Rileva il passaggio del moto attraverso lo zero
            if ((theta*thetatemp < 0)) {
                if (i < 1000000) {
                    Ttheta[i] = T - dt/2;
                    i++;
                } else {
                    break;  // Evita overflow dell'array
                }
            }

            T += dt;
        }

        // Calcola la differenza nei tempi rilevati
        for (z = i - 1; z > 0; z--) {
            Ttheta[z] = Ttheta[z] - Ttheta[z - 1];
        }

        // Somma i periodi
        for (x = 1; x < i; x++) {
            Ttot += Ttheta[x];
        }

        // Calcola il periodo medio
            Tcalcolato = Ttot / (i / 2);
       

        fprintf(fp, "%lf %lf %lf %d\n", thetai, Tcalcolato, Tanalitico, i);
        thetai += dtheta;
    }

    fclose(fp);
}


void EuleroOsc(struct d1 *d1, double m, double k, double dt, double Tmax){
    FILE *fp;
    double T=0.;

    d1->a= -k * d1->x / m;

    fp=fopen("Eulero.dat","w");
    while(T<Tmax){
        
        fprintf(fp,"%lf %lf %lf %lf %lf \n", T, d1->x, d1->v, d1->a, d1->E);
        d1->a= -k * d1->x / m;  // Aggiorna la posizione usando la posizione corrente
        d1->x = d1->x + d1->v * dt;  // Aggiorna la posizione usando la velocità corrente
        d1->v = d1->v + d1->a * dt;  // Aggiorna la velocità usando l'accelerazione corrente
        d1->E = 0.5 * m * pow(d1->v, 2) + 0.5 * k * pow(d1->x, 2); // Calcola l'energia totale del sistema
        
        T += dt;
    }
    fclose(fp);
}

void EuleroCromerOsc(struct d1 *d1, double m, double k, double dt, double Tmax){
    FILE *fp;
    double T=0.,w;

    w = sqrt(k / m);

    d1->a= aosc(w , d1->x);

    fp=fopen("EuleroCromer.dat","w");
    while(T<Tmax){
        
        fprintf(fp,"%lf %lf %lf %lf %lf \n", T, d1->x, d1->v, d1->a, d1->E);
        d1->a= aosc(w , d1->x);  // Aggiorna la posizione usando la posizione corrente
        d1->v = d1->v + d1->a * dt;  // Aggiorna la velocità usando l'accelerazione corrente
        d1->x = d1->x + d1->v * dt;  // Aggiorna la posizione usando la velocità corrente
        d1->E = 0.5 * m * pow(d1->v, 2) + 0.5 * k * pow(d1->x, 2); // Calcola l'energia totale del sistema
        
        T += dt;
    }
    fclose(fp);
}





void remove_spaces(char *str) {
    char *i = str;
    char *j = str;
    while (*j != '\0') {
        *i = *j++;
        if (*i != ' ') { 
            i++;
        }
    }
    *i = '\0'; 
}

void Verlet(struct d1 *d1, double m, double k, double dt, double Tmax) {
    double T = 0.0, x0, x1, a, w;
    w= sqrt(k /m);
    FILE *fp;

    a = (-k / m) * d1->x;
    x1 = d1->x + d1->v * dt + 0.5 * a * dt * dt;
    fp = fopen("Verlet.dat", "w");

    d1->E = 0.5 * m * pow(d1->v, 2) + 0.5 * k * pow(d1->x, 2);
    fprintf(fp, "%lf %lf %lf %lf\n", T, d1->x, d1->v, d1->E);

    x0 = d1->x;

    while (T < Tmax) {
        d1->x = 2 * x1 - x0 + a * dt * dt;
        d1->v = (d1->x - x0) / (2 * dt);
        a = aosc(w, d1->x);
        d1->E = 0.5 * m * pow(d1->v, 2) + 0.5 * k * pow(d1->x, 2);

        T += dt;
        fprintf(fp, "%lf %lf %lf %lf\n", T, d1->x, d1->v, d1->E);

        x0 = x1;
        x1 = d1->x;
    }

    fclose(fp);
}



void rk4(struct d1 *d1, double m, double k, double dt, double Tmax){

    k_t k1, k2, k3, k4;
    double T = 0.0, w;
    w= sqrt(k /m);
    FILE *fp;

    fp = fopen("rk4.dat", "w");

    d1->a = aosc(w, d1->x);
    d1->E = 0.5 * m * pow(d1->v, 2) + 0.5 * k * pow(d1->x, 2);

    while (T < Tmax) {

        
        fprintf(fp, "%lf %lf %lf %lf %lf \n", T, d1->x, d1->v, d1->a, d1->E);

        k1.x = d1->v * dt;
        k1.v = aosc(w, d1->x) * dt ;

        k2.x = (d1->v + 0.5 * k1.v) * dt;
        k2.v = aosc(w, (d1->x + 0.5 * k1.x)) * dt;

        k3.x = (d1->v + 0.5 * k2.v) * dt;
        k3.v = aosc(w, (d1->x + 0.5 * k2.x)) * dt;

        k4.x = (d1->v + k3.v) * dt;
        k4.v = aosc(w, (d1->x + k3.x)) * dt;

        d1->x += (1.0/6.0)*(k1.x + 2 * k2.x + 2 * k3.x + k4.x);
        d1->v += (1.0/6.0)*(k1.v + 2 * k2.v + 2 * k3.v + k4.v);

        
        

    
        d1->a = aosc(w, d1->x);
        d1->E = 0.5 * m * pow(d1->v, 2) + 0.5 * k * pow(d1->x, 2);

       
        T += dt;
    }

    fclose(fp);
}








void rk2(struct d1 *d1, double m, double k, double dt, double Tmax) {
    double T = 0.0;
    double w = sqrt(k / m);
    k_t k1, k2;
    FILE *fp;

    fp = fopen("rk2.dat", "w");

    while (T < Tmax) {
        d1->a = aosc(w, d1->x);
        d1->E = 0.5 * m * pow(d1->v, 2) + 0.5 * k * pow(d1->x, 2);
        fprintf(fp, "%lf %lf %lf %lf %lf\n", T, d1->x, d1->v, d1->a, d1->E);

        k1.x = d1->v * dt;
        k1.v = d1->a * dt;

        k2.x = (d1->v + 0.5 * k1.v) * dt;
        k2.v = aosc(w, d1->x + 0.5 * k1.x) * dt;

        d1->x += k2.x;
        d1->v += k2.v;

        T += dt;
    }

    fclose(fp);
}



d2_t pendolobis(struct d2 *d2, double m, double L, double beta, struct f *f, double dt, double T){

    double w, g = 1.;
    kt_t k1, k2, k3, k4;
    


    w = sqrt(g / L);


    k1.theta = d2->v * dt;
    k1.v =  apen(w, d2->theta, beta, d2->v * L, f, T, m) * dt;

    k2.theta = (d2->v + 0.5 * k1.v) * dt;
    k2.v =  apen(w, (d2->theta + 0.5 * k1.theta), beta, (d2->v + 0.5 * k1.v) * L, f, T + 0.5 * dt, m) * dt;

    k3.theta = (d2->v + 0.5 * k2.v) * dt;
    k3.v = apen(w, (d2->theta + 0.5 * k2.theta), beta, (d2->v + 0.5 * k2.v) * L, f, T + 0.5 * dt, m) * dt;

    k4.theta =  (d2->v + k3.v) * dt;
    k4.v = apen(w, (d2->theta + k3.theta), beta, (d2->v + k3.v) * L, f, T + dt, m) * dt;

    d2->theta += (1.0/6.0)*(k1.theta + 2 * k2.theta + 2 * k3.theta + k4.theta);
    d2->v += (1.0/6.0)*(k1.v + 2 * k2.v + 2 * k3.v + k4.v);

    //d2->theta = fix(d2->theta);



    d2->a = apen(w, d2->theta, beta, d2->v * L, f, T, m);
    d2->E = 0.5 * m * pow(d2->v, 2) * L + m * g * L * (1 - cos(d2->theta));

    return *d2;

}





void poincare(FILE *fp, double T, double dt, double Tmax){

    int j=0;
    double i=0.,a,b,c,d,e,f=5;
    FILE *fpc = fopen("poincare.dat","w");

    
    while(i<Tmax){
        fscanf(fp,"%lf %lf %lf %lf %lf\n",&a, &b, &c, &d, &e);
        j++;
        
        if(j % 1000 == 1){
            
            fprintf(fpc,"%lf %lf %lf\n", a, b, c);

        }
        
            
        i += dt;
    } 
    fclose(fpc);
    fclose(fp);
}



d2_t resett(d2_t pen, double *T, d2_t CI){
    pen.theta = CI.theta;
    pen.v = CI.v;
    *T = 0;
    return pen;
}

