// PRIMA DEL MAIN: librerie, costanti, funzioni
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define N 10000 // numero di punti nel mesh
#define epsilon 5.99 // parametro di LJ, in meV

double v_LJ(double r, double sigma){  // definisco il potenziale di Lennard_Jones in funzione di una posizione (punto del mesh)
    return 4*epsilon*((pow(sigma/r,12))-(pow(sigma/r,6))); // valore di LJ per ogni punto del mesh
}

double kappa_squared(double r, double E, int l, double prefactor, double sigma){  // definisco la funzione k_{E}^{2}(r) che compare nell'equazione di Shrodinger
    return prefactor*(E - v_LJ(r,sigma) - l*(l+1)/(r*r*prefactor));
}

double j(int l, double x){  // routine che implementa le funzioni di Bessel sferiche
    double j_l_minus_1 = cos(x) / x;
    double j_l = sin(x) / x;
    double j_l_plus_1;
    int i;
    j_l_plus_1 = j_l;
    for (i = 0; i < l; i++){
        j_l_plus_1 = (-j_l_minus_1 + ((2*i+1)/x)*j_l);
        j_l_minus_1 = j_l;
        j_l = j_l_plus_1;
    };
    return j_l_plus_1;
}

double n(int l, double x){  // routine che implementa le funzioni di Neumann sferiche
    double n_l_minus_1 = sin(x) / x;
    double n_l = -cos(x) / x;
    double n_l_plus_1;
    int i;
    n_l_plus_1 = n_l;
    for (i = 0; i < l; i++){
        n_l_plus_1 = (-n_l_minus_1 + ((2*i+1)/x)*n_l);
        n_l_minus_1 = n_l;
        n_l = n_l_plus_1;
    };
    return n_l_plus_1;
}


// ORA C'E' IL MAIN
int main() {

    double sigma; // parametro del potenziale di LJ
    double E = 1.5; // fissiamo l'energia dello stato, in meV
    int l, l_max; // numero quantico angolare e suo valore massimo
    printf("Enter the maximum value of the angular momentum quantum number: "); // come input da tastiera il valore massimo di momento angolare
    scanf("%d", &l_max);
    double m_H = 1.6736*pow(10,-27);  // massa dell'idrogeno in kg
    double m_Kr = 1.3915*pow(10,-25);  // massa del kripton in kg
    double mu = (m_H*m_Kr)/(m_H + m_Kr);  // massa ridotta del sistema H-Kr
    double h_bar = 1.0546*pow(10,-34);  // h_tagliato in J*s
    double conversion = 1.6022*pow(10,-42);  // fattore di conversione in meV e Angstrom
    double prefactor = (2*mu*conversion)/(pow(h_bar,2));  // fattore che compare nell'espressione di k_{E}^{2}(r)
    double h = 0.01;  // distanza tra punti adiacenti del mesh
    double r[N], v[N], k[N], u_l[N]; // array di punti del mesh, potenziale, k_{E}^{2}(r) e funzioni d'onda, fissato l            
    double valori_r1[5], valori_r2[5], valori_ul1[5], valori_ul2[5];  // array dei punti (r1, r2), e corrispondenti funzioni d'onda, in cui calcoliamo le phase shift
    double beta[l_max+1][5], tgdelta[l_max+1][5], delta[l_max+1][5];  // matrici per il fattore beta, tangente della phase shift e phase shift stessa, in funzione di l e punti (r1,r2)
    double k_wave = sqrt(prefactor*E);  // vettore d'onda
    double pi = 4*atan(1);  // pi greco
    double b;  // parametro per inizializzare le funzioni d'onda
    int pos1, pos2;  // posizioni, negli array, dei punti per calcolare la phase shift al variare dell'energia
    int dim = 10;  // numero di parametri sigma di LJ
    double delta_mat[dim+1];  // vettore per i valori Delta^2(sigma) 

    FILE *output_es_1; // funzioni d'onda radiali al variare di sigma
    output_es_1 = fopen("Dati_es_1","w");  

    FILE *output5_es_1; // phase shift al variare di sigma
    output5_es_1 = fopen("Dati_es_1_phase_shift","w");  

    FILE *output10_es_1; // cross section al variare di sigma per l_max = 6
    output10_es_1 = fopen("Dati_es_1_cross_section_per_ogni_l","w");
    
    FILE *output6_es_1; // cross section al variare di sigma per l_max = 6
    output6_es_1 = fopen("Dati_es_1_cross_section6","w");  

    FILE *output68_es_1; // cross section al variare di sigma per l_max = 6
    output68_es_1 = fopen("Dati_es_1_cross_section8","w");  

    FILE *output20_es_1; // varianza Delta^2 al variare di sigma
    output20_es_1 = fopen("Dati_es_1_Delta_squared","w");  

    FILE *output7_es_1;  // funzioni di Bessel al variare di l
    output7_es_1 = fopen("Dati_es_1_Bessel","w");  

    FILE *output8_es_1; // funzioni di Neumann al variare di l
    output8_es_1 = fopen("Dati_es_1_Neumann","w");  

    int N3 = 2000; // numero di punti che scegliamo per stampare le funzioni di Bessel e Neumann
    double x[N3];  // array degli argomenti in cui valutare le funzioni

    for(int w=0; w<N3; w++){  // ciclo per stampare le funzioni al variare di l
        x[w] = w*h;
        fprintf(output7_es_1, "%3.4e ", x[w]);
        fprintf(output8_es_1, "%3.4e ", x[w]);
        for(l=0; l<=l_max; l++){
            fprintf(output7_es_1, "%3.4e ", j(l,x[w]));
            fprintf(output8_es_1, "%3.4e ", n(l,x[w]));
        };
        fprintf(output7_es_1, "\n");
        fprintf(output8_es_1, "\n");
    };


// !!! DI SEGUITO, I VARI STEP DELL'ESERCIZIO, ANCHE SE VARIAMO SIGMA FIN DA SUBITO !!!

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// STEP 1-2: SOLUZIONI DELL'EQUAZIONE DI SCRHODINGER E PHASE SHIFT 

    
for(int o=0; o<=dim; o++){  // ciclo più esterno che scorre sui possibili valori di sigma (si chiude poco prima della fine del main)
        
    sigma = 3.14 + o*0.01;  // valori di sigma in LJ
    pos1 = (10*sigma - (sigma/2)) / h;  // posizione del punto che corrisponde a 10*sigma, avendo sigma/2 come primo punto del mesh
    pos2 = ((10*sigma+0.5) -(sigma/2)) / h;  // posizione del punto che corrisponde a 10*sigma+0.5
    b = pow(prefactor * (2./5) * epsilon, 1./10) * pow(sigma, 6./5);  // definisco b a partire dal limite r --> 0
    printf("b = %lf \n", b); // stampiamo b per verificare il suo valore


    for(l=0; l<=l_max; l++){  // ciclo su l per ottenere le soluzioni u_{l}(r)

        r[0] = sigma/2;  // primo punto del mesh 
        r[1] = r[0] + h;  // punto successivo
        v[0] = v_LJ(r[0], sigma);
        v[1] = v_LJ(r[1], sigma);
        k[0] = kappa_squared(r[0], E, l, prefactor, sigma);
        k[1] = kappa_squared(r[1], E, l, prefactor, sigma);        
        u_l[0] = pow(r[0],2)*exp(pow((-b/r[0]),5)); // condizioni iniziali per i primi due punti
        u_l[1] = pow(r[1],2)*exp(pow((-b/r[1]),5));
        for(int i=2; i<N; i++){  // ciclo sul mesh per ottenere le precedenti variabili nei punti successivi
            r[i] = r[0] + h*i;
            v[i] = v_LJ(r[i], sigma);
            k[i] = kappa_squared(r[i], E, l, prefactor, sigma);
            u_l[i] = ( (u_l[i-1] * (2 - 5./6 * h*h * k[i-1])) - (u_l[i-2] * (1 + h*h*k[i-2]/12)) ) / (1 + h*h*k[i]/12); // algoritmo di Numerov            
        };
        
        for(int i=0; i<N; i++){  // ciclo per stampare posizioni e funzioni d'onda, compresi i primi due punti
            fprintf(output_es_1, "%3.4e %3.4e \n", r[i], u_l[i]);
        };
        
        for(int i=0; i<N; i++){  // ciclo per individuare le 5 coppie di punti (r1,r2) in cui calcolare le phase shifts
            if(i == pos1){
                for(int m=0; m<5; m++){
                valori_r1[m] = r[i];  // il punto r1 è fisso nelle 5 coppie 
                valori_r2[m] = r[i+(10*m+10)];  // il punto r2 aumenta di 0.1 Angstrom ogni volta
                valori_ul1[m] = u_l[i];  // corrispondenti funzioni d'onda
                valori_ul2[m] = u_l[i+(10*m+10)];
                };
            };
        };
    
        for(int m=0; m<5; m++){  // ciclo per calcolare beta, tan(delta) e delta nelle 5 coppie di punti
            beta[l][m] = ((valori_ul1[m])*valori_r2[m])/((valori_ul2[m])*valori_r1[m]);
            tgdelta[l][m] = (beta[l][m]*j(l,k_wave*valori_r2[m])-j(l,k_wave*valori_r1[m]))/(beta[l][m]*n(l,k_wave*valori_r2[m])-n(l,k_wave*valori_r1[m]));
            delta[l][m] = atan(tgdelta[l][m]);         
        };

    }; // si chiude il ciclo su l

    for(int m=0; m<5; m++){  // ciclo per stampare i valori di r2 e delle corrispondenti phase shift, al variare di l
        fprintf(output5_es_1, "%3.4e ", valori_r2[m]);
        for(l=0; l<=l_max; l++){
            fprintf(output5_es_1, "%3.4e ", delta[l][m]);
        };
        fprintf(output5_es_1, "\n");
    };


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// STEP 3: CROSS SECTION AL VARIARE DELL'ENERGIA

    int N2 = 350; // valori di energia
    double Emax = 3.5; // energia massima in meV (la minima è zero)
    double h_en = Emax / N2; // intervallo fra i valori di energia
    double r1 = valori_r1[4], r2 = valori_r2[4];  // fissiamo i valori su cui calcoliamo beta e le phase shift agli ultimi punti dei vettori di cui sopra
    double E2[N2], k_wave2[N2]; // array dei valori di energie e vettori d'onda
    double k2[N], u_l2[N]; // array dei nuovi k_{E}^{2}(r) e funzioni d'onda, al variare dell'energia
    double u_l2_1[N2][l_max+1], u_l2_2[N2][l_max+1]; // matrici dei valori della funzione d'onda in r1 e r2 al variare di energia e momento angolare l
    double beta2[N2][l_max+1], tgdelta2[N2][l_max+1], delta2[N2][l_max+1]; // matrici dei valori di beta, tan(delta) e delta, al variare di energia e l
    double argument[N2];  // array della sommatoria che compare nell'espressione della cross section, al variare dell'energia
    double argument2[N2];
    double cross_section[N2];  // array dei valori di cross section per ogni energia
    double cross_section2[N2];


    for(int s=1; s<=N2; s++){  // ciclo sui valori di energia

        E2[s] = s*h_en; // riempiamo gli array di energie e vettori d'onda
        k_wave2[s] = sqrt(prefactor*E2[s]);

        for(l=0; l<=l_max; l++){  // di nuovo il ciclo su l per trovare u_{l}(r) e phase shift, come per step 1, ma ora per ogni energia
        
            r[0] = sigma/2; 
            r[1] = r[0] + h;
            v[0] = v_LJ(r[0], sigma);
            v[1] = v_LJ(r[1], sigma);
            k2[0] = kappa_squared(r[0], E2[s], l, prefactor, sigma);
            k2[1] = kappa_squared(r[1], E2[s], l, prefactor, sigma);        
            u_l2[0] = pow(r[0],2)*exp(pow((-b/r[0]),5));
            u_l2[1] = pow(r[1],2)*exp(pow((-b/r[1]),5));

            for(int i=2; i<N; i++){
                r[i] = r[0] + h*i;
                v[i] = v_LJ(r[i], sigma);
                k2[i] = kappa_squared(r[i], E2[s], l, prefactor, sigma);
                u_l2[i] = ( (u_l2[i-1] * (2 - 5./6 * h*h * k2[i-1])) - (u_l2[i-2] * (1 + h*h*k2[i-2]/12)) ) / (1 + h*h*k2[i]/12);          
            };
            
            for(int i=0; i<N; i++){  // stavolta, per calcolare le phase shift, scelgo solo le distanze (r1,r2) localizzate alle posizioni (pos1,pos2)
                if(i == pos1){
                    u_l2_1[s][l] = u_l2[i];
                } else if(i == pos2) {
                    u_l2_2[s][l] = u_l2[i];
                };
            };

            beta2[s][l] = ((u_l2_1[s][l])*r2)/((u_l2_2[s][l])*r1);
            tgdelta2[s][l] = (beta2[s][l]*j(l,k_wave2[s]*r2)-j(l,k_wave2[s]*r1))/(beta2[s][l]*n(l,k_wave2[s]*r2)-n(l,k_wave2[s]*r1));
            delta2[s][l] = atan(tgdelta2[s][l]);  // calcolo le phase shift in funzione di l, fissati (r1,r2), per ogni energia

        };  // si chiude il ciclo su l

        argument[s] = 0;  // inizializzo a zero la sommatoria dell'espressione per la cross section


        for(l=0; l<=l_max; l++){  // per ogni energia, calcolo la sommatoria aggiungendo ogni termine corrispondente ai vari l
            argument[s] += (2*l+1) * pow(sin(delta2[s][l]),2);
        };

        cross_section[s] = (4*pi*argument[s]) / (pow(k_wave2[s],2));  // cross section per ogni energia

        if(l_max == 6){  // salvo su file diversi a seconda del valore inserito per l_max
            fprintf(output6_es_1, "%3.4e %3.4e \n", E2[s], cross_section[s]);
        } else {
            fprintf(output68_es_1, "%3.4e %3.4e \n", E2[s], cross_section[s]);
        };    

    };  // si chiude il ciclo sulle energie

     
    for(l=0; l<=l_max; l++){ 
        for(int s=1; s<=N2; s++){   
            E2[s] = s*h_en; // riempiamo gli array di energie e vettori d'onda
            k_wave2[s] = sqrt(prefactor*E2[s]);
            argument2[s] = 0;
            argument2[s] = (2*l+1) * pow(sin(delta2[s][l]),2);
            cross_section2[s] = (4*pi*argument2[s]) / (pow(k_wave2[s],2));
            fprintf(output10_es_1, "%3.4e %3.4e \n", E2[s], cross_section2[s]);
        };
    };


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// STEP 4: PICCHI E DELTA^{2}(SIGMA)
    
    double E2_1[3]; // array dei 3 valori di energia ai massimi 
    double maxRelativi[3]; // array del massimo assoluto e dei 2 massimi relativi 
    int counter = 0; // indice per trovare le energie ai massimi
    int numMaxRelativi = 0; // indice per trovare i massimi     

    for (int i = 0; i < N2; i++) {  // ciclo per trovare i massimi, stampo punti di massimo e massimi
        if (cross_section[i] > cross_section[i-1] && cross_section[i] > cross_section[i+1]) {
            maxRelativi[numMaxRelativi] = cross_section[i];
            E2_1[counter] = E2[i];
            printf("Energy of the relative peak = %3.4e  relative peak = %3.4e \n", E2_1[counter], maxRelativi[numMaxRelativi]);
            numMaxRelativi++;
            counter++;
        };   
    };

    double energie_exp[3] = {0.50, 1.59, 2.94}; // valori sperimentali di energie ai picchi
    double dE_exp[3] = {0.02, 0.06, 0.12}; // errori sperimentali sulle precedenti energie
    double delta_squared = 0; // inizializzo il Delta^{2}(sigma)

    for(int i=0; i<=2; i++){  // per ogni sigma, calcolo la sommatoria aggiungendo ogni termine corrispondente ai vari picchi
        delta_squared += (pow((E2_1[i]-energie_exp[i]),2) / pow(dE_exp[i],2));
    };

    delta_mat[o] = delta_squared;  // array in cui salvo i vari Delta^{2} per ogni sigma

    fprintf(output20_es_1, "%3.4e %3.4e \n", sigma, delta_mat[o]);  // stampo i valori di sigma e i corrispondenti Delta^{2}

}; // si chiude il ciclo sui valori di sigma 


fclose(output_es_1);  // chiudo tutti i file
fclose(output5_es_1);
fclose(output6_es_1);
fclose(output68_es_1);
fclose(output20_es_1);
fclose(output7_es_1);
fclose(output8_es_1);
fclose(output10_es_1);


return 0;


}; // si chiude il main