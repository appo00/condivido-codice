// PRIMA DEL MAIN: librerie, variabili globali, costanti, funzioni
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define N 10000 // numero di punti nel mesh
#define epsilon 5.99 // parametro di LJ, in meV
#define sigma 3.18 // parametro di LJ, in Angstrom, anche se bisogna cambiarlo
//#define prefactor 0.47646 // 2*mu/(h_bar)^2 in (meV * A°^2)^-1  (VALORE VERO 0.47646)

double v_LJ(double r){  // definisco il potenziale di Lennard_Jones in funzione di una posizione (punto del mesh)
    return 4*epsilon*((pow(sigma/r,12))-(pow(sigma/r,6))); // valore di LJ per ogni punto del mesh
}

double kappa_squared(double r, double E, int l, double prefactor){  // definisco la funzione k_E(r) che compare nell'equazione di Shrodinger, in unità atomiche, ovvero m = h_bar = 1
    return prefactor*(E - v_LJ(r) - l*(l+1)/(r*r*prefactor));
}

double j(int l, double x)
{
    double j_l_minus_1 = cos(x) / x;
    double j_l = sin(x) / x;
    double j_l_plus_1;
    int i;
    j_l_plus_1 = j_l;
/*     if (l == 0)
    {
        return j_l;
    }
    else
    { */
        for (i = 0; i < l; i++)
        {
            j_l_plus_1 = (-j_l_minus_1 + ((2*i+1)/x)*j_l);
            j_l_minus_1 = j_l;
            j_l = j_l_plus_1;
        }
        return j_l_plus_1;
    //}
}

double n(int l, double x)
{
    double n_l_minus_1 = sin(x) / x;
    double n_l = -cos(x) / x;
    double n_l_plus_1;
    int i;
    n_l_plus_1 = n_l;
    /* if (l == 0)
    {
        return n_l;
    }
    else
    { */
        for (i = 0; i < l; i++)
        {
            n_l_plus_1 = (-n_l_minus_1 + ((2*i+1)/x)*n_l);
            n_l_minus_1 = n_l;
            n_l = n_l_plus_1;
        }
        return n_l_plus_1;
    //}
}




// ORA C'E' IL MAIN
int main() {
    int n_sigma = 10;
    double sigma_v[n_sigma+1];
    double E = 1.5; // fissiamo l'energia dello stato, in meV
    int l, l_max = 6; // valore massimo del momento angolare orbitale
    double m_H = 1.6736*pow(10,-27);   //kg
    double m_Kr = 1.3915*pow(10,-25); //kg
    double mu = (m_H*m_Kr)/(m_H + m_Kr);
    double h_bar = 1.0546*pow(10,-34); //J*s
    double conversion = 1.6022*pow(10,-42); //per convertire in meV e amstrong
    double prefactor = (2*mu*conversion)/(pow(h_bar,2));
    // definisco l'algoritmo che mi calcola l'equazione di Shrodinger radiale date le condizioni iniziali, per un fissato valore di l

    double r_max = 100.; // estremo superiore del mesh, in Angstrom
    double h = r_max/N; // ampiezza del meshino
    double r[N], v[N], k[N]; // array di punti del mesh, potenziale e k^2
    double u_l[N]; // array delle soluzioni radiali ad ogni punto, fissato l    
    int pos1 = (10*sigma - (sigma/2)) / h;
    int pos2 = ((10*sigma+0.5) -(sigma/2)) / h;
    // 1° FILE: DATI AL VARIARE DI L

    FILE *output_es_1; // definisco un file dove vado a salvare i dati
    output_es_1 = fopen("Dati_es_1.txt","w");  // apro il file e ci scrivo

    //int i; // indice che scorre sui punti del mesh
    double b =  pow(2, 1./6)* sigma;  // costante che appare nell'esponente dell'espressione asintotica per u(r)
    printf("%lf \n", b);

    // 2° FILE: VALORI DI BETA AL VARIARE DI L, PER LA PHASE SHIFT

    
    double valori_r1[5]; //vettore per valori r1
    double valori_r2[5]; //vettore per valori r2
    double valori_ul1[5];
    double valori_ul2[5];
    int m;
    double beta[l_max+1][5]; // vettore dei valori di beta 
    double pi = 4*atan(1);
    double tgdelta[l_max+1][5], delta[l_max+1][5];
    double k_wave = sqrt(prefactor*E);

    FILE *output5_es_1; // definisco un file dove vado a salvare i dati
    output5_es_1 = fopen("Dati_es_1_phase_shift.txt","w");  // apro il file e ci scrivo

    for(l=0; l<=l_max; l++){

        r[0] = sigma/2; 
        r[1] = r[0] + h;
        v[0] = v_LJ(r[0]);
        v[1] = v_LJ(r[1]);
        k[0] = kappa_squared(r[0], E, l, prefactor);
        k[1] = kappa_squared(r[1], E, l, prefactor);        
        u_l[0] = pow(r[0],2)*exp(pow((-b/r[0]),5)); // condizioni iniziali date dal prof (?)
        u_l[1] = pow(r[1],2)*exp(pow((-b/r[1]),5));

        for(int i=2; i<N; i++){  // ciclo per implementare Numerov
            r[i] = r[0] + h*i; // punto del mesh, come distanza dall'origine
            v[i] = v_LJ(r[i]);
            k[i] = kappa_squared(r[i], E, l, prefactor);
            u_l[i] = ( (u_l[i-1] * (2 - 5./6 * h*h * k[i-1])) - (u_l[i-2] * (1 + h*h*k[i-2]/12)) ) / (1 + h*h*k[i]/12); // algoritmo di Numerov            
        };

        for(int i=0; i<N; i++){  // ciclo per stampare tutto su file, compresi i primi due punti
            fprintf(output_es_1, "%3.4e %3.4e \n", r[i], u_l[i]);
        };
    
        
        
        for(int i=0; i<N; i++){ 
            if(i == pos1){
                for(m=0; m<5; m++){
                valori_r1[m] = r[i];
                valori_r2[m] = r[i+(10*m+10)]; //abbiamo modificato la distanza degli r da i+4 a i+10 
                valori_ul1[m] = u_l[i];
                valori_ul2[m] = u_l[i+(10*m+10)];
                };
            };
        };
    
        for(m=0; m<5; m++){
            beta[l][m] = ((valori_ul1[m])*valori_r2[m])/((valori_ul2[m])*valori_r1[m]);
            tgdelta[l][m] = (beta[l][m]*j(l,k_wave*valori_r2[m])-j(l,k_wave*valori_r1[m]))/(beta[l][m]*n(l,k_wave*valori_r2[m])-n(l,k_wave*valori_r1[m]));
            delta[l][m] = atan(tgdelta[l][m]);
            fprintf(output5_es_1, "%3.4e ", delta[l][m]);            
        };
            fprintf(output5_es_1, "\n");
    };
    fclose(output_es_1);
    fclose(output5_es_1);
   


// STAMPIAMO LE BESSEL E NEUMANN FUNCTION E POI LE CONFRONTEREMO CON I GRAFICI CHE TROVIAMO ONLINE

FILE *output7_es_1; // definisco un file dove vado a salvare i dati
output7_es_1 = fopen("Dati_es_1_Bessel.txt","w");  // apro il file e ci scrivo

FILE *output8_es_1; // definisco un file dove vado a salvare i dati
output8_es_1 = fopen("Dati_es_1_Neuman.txt","w");  // apro il file e ci scrivo

int N3 = 2000; // numero di punti che scegliamo per plottare le Bessel

double x[N3];
for(int w=0; w<N3; w++){
    x[w] = w*h;
    fprintf(output7_es_1, "%3.4e ", x[w]);
    for(l=0; l<=l_max; l++){
        fprintf(output7_es_1, "%3.4e ", j(l,x[w]));
    };
    fprintf(output7_es_1, "\n");
};
fclose(output7_es_1);

double xx[N3];
for(int w=0; w<N3; w++){
    xx[w] = w*h;
    fprintf(output8_es_1, "%3.4e ", xx[w]);
    for(l=0; l<=l_max; l++){
        fprintf(output8_es_1, "%3.4e ", n(l,xx[w]));
    };
    fprintf(output8_es_1, "\n");
};
fclose(output8_es_1);


// STEP 3

int N2 = 350; // valori di energia
double Emax = 3.5; // energia massima in meV (la minima è zero)
double h_en = Emax / N2; // intervallo fra i valori di energia
double r1 = valori_r1[4];  // fissiamo i valori su cui calcoliamo beta (e la phase shift) agli ultimi punti dei vettori di sopra
double r2 = valori_r2[4];
double E2[N2][n_sigma+1]; // vettore dei 351 valori di energie
int s; // indice che scorre sulle energie
double k2[N][n_sigma+1], u_l2[N][n_sigma+1]; // vettore dei nuovi k^2 e funzioni d'onda
double u_l2_1[N2][l_max+1][n_sigma+1]; // vettore dei valori della funzione d'onda nel punto r1 al variare di l
double u_l2_2[N2][l_max+1][n_sigma+1]; // vettore dei valori della funzione d'onda nel punto r2 al variare di l
double beta2[N2][l_max+1][n_sigma+1]; // vettore dei valori di beta 
double tgdelta2[N2][l_max+1][n_sigma+1], delta2[N2][l_max+1][n_sigma+1];
double k_wave2[N2];
double argument[N2][n_sigma+1];
double cross_section[N2][n_sigma+1];
double r_v[N][n_sigma+1], v_v[N][n_sigma+1], k_v[N][n_sigma+1];


double E2_peak1, E2_peak2, E2_peak3;  // valori di energia corrispondenti ai 3 massimi
int a1, a2, a3; // posizioni dei 3 massimi di energia nell'array
double peak1, peak2, peak3; // 3 massimi di cross section
double b_v[n_sigma+1];
//int sigma_v1[n_sigma+1];

FILE *output6_es_1; // definisco un file dove vado a salvare i dati
output6_es_1 = fopen("Dati_es_1_cross_section.txt","w");  // apro il file e ci scrivo
FILE *output20_es_1; // definisco un file dove vado a salvare i dati
output20_es_1 = fopen("Beta.txt","w");  // apro il file e ci scrivo
FILE *output9_es_1; // definisco un file dove vado a salvare i dati
output9_es_1 = fopen("Delta.txt","w");  // apro il file e ci scrivo

int pos1_2[n_sigma +1];
int pos2_2[n_sigma +1];

for(int c = 0; c<=n_sigma; c++ ){
    argument[s = 0][c] = 0;
    sigma_v[c] = 3.10 + 0.01*c; 
    pos1_2[c] = (10*sigma_v[c] - (sigma_v[c]/2)) / h;
    pos2_2[c] = ((10*sigma_v[c] + 0.5) - (sigma_v[c]/2)) / h;
    //sigma_v1[c] = 10*sigma_v[c];

b_v[c] =  pow(2, 1./6)* sigma_v[c];
//printf("%f \n", b_v[c]);

for(s=1; s<=N2; s++){

    E2[s][c] = s*h_en; // l'intervallo di energie alla fine è lo stesso dell'intervallo nel mesh sulle distanze
    k_wave2[s] = sqrt(prefactor*E2[s][c]);

    for(l=0; l<=l_max; l++){
        
        r_v[0][c] = sigma_v[c]/2; 
        r_v[1][c] = r_v[0][c] + h;
        v_v[0][c] = v_LJ(r_v[0][c]);
        v_v[1][c] = v_LJ(r_v[1][c]);
        k2[0][c] = kappa_squared(r_v[0][c], E2[s][c], l, prefactor);
        k2[1][c] = kappa_squared(r_v[1][c], E2[s][c], l, prefactor);        
        u_l2[0][c] = pow(r_v[0][c],2)*exp(pow((-b_v[c]/r_v[0][c]),5)); // condizioni iniziali date dal prof (?)
        u_l2[1][c] = pow(r_v[1][c],2)*exp(pow((-b_v[c]/r_v[1][c]),5));
        
        

        for(int i=2; i<N; i++){  // ciclo per implementare Numerov
            r_v[i][c] = r_v[0][c] + h*i; // punto del mesh, come distanza dall'origine
            v_v[i][c] = v_LJ(r_v[i][c]);
            k2[i][c] = kappa_squared(r_v[i][c], E2[s][c], l, prefactor);
            u_l2[i][c] = ( (u_l2[i-1][c] * (2. - 5./6 * h*h * k2[i-1][c])) - (u_l2[i-2][c] * (1. + h*h*k2[i-2][c]/12)) ) / (1. + h*h*k2[i][c]/12); // algoritmo di Numerov     
            //printf("%3.4e \n ", u_l2[i][c]);
        };
        
        for(int i=0; i<N; i++){
            if(  i == pos1_2[c] ){
                u_l2_1[s][l][c] = u_l2[i][c]; // riempio il vettore dei valori della funzione d'onda a r1 al variare di l
            }
            else if(  i == pos2_2[c] ){
                u_l2_2[s][l][c] = u_l2[i][c]; // riempio il vettore dei valori della funzione d'onda a r2 al variare di l
            };

        };

        beta2[s][l][c] = ((u_l2_1[s][l][c])*(10*sigma_v[c] + 0.5))/((u_l2_2[s][l][c])*(10*sigma_v[c]));
        tgdelta2[s][l][c] = (beta2[s][l][c]*j(l,k_wave2[s]*(10*sigma_v[c] + 0.5))-j(l,k_wave2[s]*(10*sigma_v[c])))/(beta2[s][l][c]*n(l,k_wave2[s]*(10*sigma_v[c] + 0.5))-n(l,k_wave2[s]*(10*sigma_v[c])));
        delta2[s][l][c] = atan(tgdelta2[s][l][c]);

        fprintf(output20_es_1, "%3.4e ", beta2[s][l][c]);

    };

        fprintf(output20_es_1, "\n");

    
    for(l=0; l<=l_max; l++){
        argument[s][c] += (2*l+1) * pow(sin(delta2[s][l][c]),2);
    };

    cross_section[s][c] = (4*pi*argument[s][c]) / (pow(k_wave2[s],2));
    fprintf(output6_es_1, "%3.4e %3.4e \n", E2[s][c], cross_section[s][c]);

};
fprintf(output20_es_1, "\n");

};




// STEP 4: TROVIAMO I PICCHI
double maxRelativi[N2][n_sigma+1]; // massimi relativi
int numMaxRelativi; // numero di massimi relativi trovati 
double E2_1[3][n_sigma+1];
int counter;
// Trova i massimi relativi
for(int c=0; c<=n_sigma; c++){
    numMaxRelativi = 0;
    counter = 0;
for (int i = 1; i <= N2; i++) {
    if (cross_section[i][c] > cross_section[i-1][c] && cross_section[i][c] > cross_section[i+1][c]) {
        maxRelativi[numMaxRelativi][c] = cross_section[i][c];
        //printf("Energy of the relative peak = %3.4e  relative peak = %3.4e \n", E2[i], cross_section[i][c]);
        E2_1[counter][c] = E2[i][c];
        numMaxRelativi++;
        counter++;
    }   
}
};


double delta_squared[n_sigma+1];
double energie_exp[3] = {0.50, 1.59, 2.94}; 
double dE_exp[3] = {0.02, 0.06, 0.12};
for(int c=0; c<=n_sigma; c++){
delta_squared[c]=0;
for(int i=0; i<=2; i++){
    delta_squared[c] += (pow((E2_1[i][c]-energie_exp[i]),2)); // / pow(dE_exp[i],2));
    printf("energie_max = %3.4e \n", E2_1[i][c]);
    printf("energie_exp = %3.4e \n", energie_exp[i]);
            };
fprintf(output9_es_1, "%3.4e \n", delta_squared[o]);
        };
/* fclose(output6_es_1);
fclose(output9_es_1) */

    return 0;
};
