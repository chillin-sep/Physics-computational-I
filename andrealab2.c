#include <stdio.h>
#include <math.h>
#include <stdlib.h>
 
typedef struct phases{ double T,x,xp,v,a,E,k,m,dt;int value;}phases_t;


struct phases reset(struct phases *p, FILE *a){
fscanf(a,"%lf  %lf  %lf  %lf  %lf ",&p->T,&p->x,&p->v,&p->k,&p->m);
p->a=(-p->k/p->m) *p->x; 
p->E=0.5 * p->m * p->v * p->v + 0.5 * p->k * p->x * p->x;

fclose(a);

return *p;
};

struct phases Eulero(struct phases *p){


  p->a = (-p->k/p->m) *p->x; 
  p->x = p->x + p->v * p->dt;                                                                                                 
  p->v = p->v + p->a * p->dt;    
  p->E= 0.5 * p->m * p->v * p->v + 0.5 * p->k * p->x * p->x;                                                                     

  p->value = 0;  

  return *p;
};

struct phases Cromer(struct phases *p){
 
  p->a = (-p->k/p->m) *p->x; 
  p->v = p->v + p->a * p->dt; 
  p->x = p->x + p->v * p->dt;                                                                                                   
  p->E= 0.5 * p->m * p->v * p->v + 0.5 * p->k * p->x * p->x;     

  p->value = 0; 

  return *p ;                                                                   
};

struct phases Verlet(struct phases *p){

  double xn;

  p->a = (-p->k/p->m) * p->x; // calcola la nuova accelerazione
  xn = 2 * p->x - p->xp + p->a * p->dt * p->dt; // posizione successiva
  p->v=(xn - p->xp)/2*p->dt;
  p->E= 0.5 * p->m * p->v * p->v + 0.5 * p->k * p->x * p->x;     

  p->xp=p->x;
  p->x=xn;

  p->value = 1;

  return *p;
};

struct phases Velocity_Verlet(struct phases *p){

 double ap;

 ap = p->a;

 p->x = p->x + p->v * p->dt + 0.5 * p->a * p->dt * p->dt;
 p->a = (-p->k/p->m) * p->x;
 p->v = p->v + (p->a + ap) * p->dt/2;
 p->E= 0.5 * p->m * p->v * p->v + 0.5 * p->k * p->x * p->x;

 p->value = 0;

 return *p;

};

double periodo(struct phases *p,FILE *f){

 int i,j=0;
 int counts = (int) (p->T/p->dt);
 double periodi[counts];
 double P;
 double x,v,a,E,t;
 double x0,t0;
 
 fscanf(f,"%lf %lf %lf %lf %lf\n",&x,&v,&a,&E,&t);

 for(i=0;i<counts;i++)
 {
  x0=x;
  t0=t;
  fscanf(f,"%lf %lf %lf %lf %lf\n",&x,&v,&a,&E,&t);

  if(x0*x<0 && x0>x) //controllo dei punti in cui la posizione cambia segno sempre nello stesso verso
  {
    if(t0!=0){
      periodi[j]=t0;
      j++;
    }
  }          
 }


 for(i=1;i<j;i++)
 {
   P = P + (periodi[i]-periodi[i-1]);
 }

 P= P/(j-1);

 fclose(f);

 return P;
};

void dEvsdt(struct phases *p, struct phases (*integration)(struct phases *s), FILE *a, FILE *b) {
    double t = 0.0;
    double E0 = 0.5 * p->m * p->v * p->v + 0.5 * p->k * p->x * p->x;
    double dt_i=0.001;
    int i=1000;
    double deltaE;
    
    rewind(b);
    reset(p, b);
    integration(p); 
    rewind(b);
    reset(p, b); 

   p->dt=dt_i;

  if(p->value == 0){
    while (p->dt <= 1. )
     {
        t = 0.;
      
        rewind(b);  
        reset(p, b); 
     
        while (t < p->T) 
        {
            integration(p);
            t += p->dt;
        }
       
       deltaE = (p->E) - E0;
        
        
        if (deltaE == INFINITY) break;
        fprintf(a, "%lf\t%lf\n", deltaE, p->dt ); 
        p->dt += 0.001;
        
      
       
    };
  };
 
  if(p->value == 1)
  {
   while (p->dt <= 1.0 )
   {
    t = 0.0;
      
        rewind(b);  
        reset(p, b); 
        
        Cromer(p);  
        t=t+(p->dt);
        while (t < p->T) {
            integration(p);
            t += p->dt;
        }

       
        double deltaE = (p->E) - E0;
        
        
        if (deltaE == INFINITY) break;
         fprintf(a, "%lf\t%lf\n", deltaE, p->dt  ) ;   
        p->dt += 0.001;
       
      
                                                                                                                     
  };
 };

  
 fclose(a);
 fclose(b);
};

void dEvsdt1(struct phases *p, struct phases (*integration)(struct phases *s), FILE *a, FILE *b) {
    double t = 0.0;
    double E0 = 0.5 * p->m * p->v * p->v + 0.5 * p->k * p->x * p->x;
    double E, deltaE;
    double dt = 0.001;

    // Esegui una chiamata iniziale a integration per impostare p->value
    rewind(b);
    reset(p,b);
    integration(p);
    rewind(b);
    reset(p,b);

    if (p->value == 0) {  
        while (dt <= 1.0) {
            t = 0.0;

            rewind(b);
            reset(p, b);

        
            while (t < p->T) {
                p->a = (-p->k / p->m) * p->x;
                p->x = p->x + p->v * dt;
                p->v = p->v + p->a * dt;
                t += dt;
            }

          
            E = 0.5 * p->m * p->v * p->v + 0.5 * p->k * p->x * p->x;

            
            if (E == INFINITY) break;

       
            deltaE = E - E0;
            
            
            fprintf(a, "%lf\t%lf\n", deltaE, dt);
            
            
            dt += 0.001;
        }
    }

    fclose(a);
    fclose(b);
};



int main(){

double t=0.;

FILE *fp,*fr;

phases_t ph;



/*******************SETTING DEI DATI INIZIALI***************************/
fp=fopen("dati.dat","r");
reset(&ph,fp);

ph.dt=0.001;
t=0.;

/******************************EULERO**********************************/
fp=fopen("eu.dat","w");
while(t<=ph.T)
{
 fprintf(fp,"%f\t%f\t%f\t%f\t%f\n",ph.x,ph.v,ph.a,ph.E,t);
 Eulero(&ph);
 t=t+ph.dt;
}
fclose(fp);

/*****************************RESET DEI DATI**************************/
fp=fopen("dati.dat","r");
reset(&ph,fp);
t=0.;
ph.dt=0.001;

/***************************EULERO-CROMER******************************/
fp=fopen("cromer.dat","w");
while(t<=ph.T)
{
 fprintf(fp,"%f\t%f\t%f\t%f\t%f\n",ph.x,ph.v,ph.a,ph.E,t);
 Cromer(&ph);
 t=t+ph.dt;
}
fclose(fp);

/*****************************RESET DEI DATI**************************/
fp=fopen("dati.dat","r");
reset(&ph,fp);
t=0.;
ph.dt=0.001;

/*******************************VERLET********************************/
fp=fopen("verlet.dat","w");
ph.xp=ph.x;
fprintf(fp,"%f\t%f\t%f\t%f\t%f\n",ph.x,ph.v,ph.a,ph.E,t);
Cromer(&ph);
t=t+ph.dt;
while(t<=ph.T)
 {
  fprintf(fp,"%f\t%f\t%f\t%f\t%f\n",ph.x,ph.v,ph.a,ph.E,t);
  Verlet(&ph);
  t=t+ph.dt;
  
 }
fclose(fp);

/*****************************RESET DEI DATI**************************/
fp=fopen("dati.dat","r");
reset(&ph,fp);
t=0.;
ph.dt=0.001;

/*******************************VELOCITY-VERLET********************************/
fp=fopen("velocity.dat","w");

while(t<=ph.T){
  fprintf(fp,"%f\t%f\t%f\t%f\t%f\n",ph.x,ph.v,ph.a,ph.E,t);
  Velocity_Verlet(&ph);

  t=t+ph.dt;
}

fclose(fp);
/*****************************RESET DEI DATI**************************/
fp=fopen("dati.dat","r");
reset(&ph,fp);
t=0.;
ph.dt=0.001;

/******************************* PERIODI ********************************/
fp=fopen("velocity.dat","r");
printf("%f\n",periodo(&ph,fp));

/*****************************RESET DEI DATI**************************/
fp=fopen("dati.dat","r");
reset(&ph,fp);
t=0.;
ph.dt=0.001;

/************** ENERGIA IN FUNZIONE DEL PASSO DI INTEGRAZIONE *******************/

fp=fopen("energia_e.dat","w");
fr=fopen("dati.dat","r");
dEvsdt(&ph,Eulero,fp,fr);



}
