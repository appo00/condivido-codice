#include <stdio.h>

main()
{
    int I, J;
    
    /* per I che va da 1 a 10*/
    for( I = 1 ; I <= 10 ; I = J + 1 ){
        /* per J che va da 1 a 10*/
        for( J = 1 ; J <= 10 ; J = J + 1 )
            /* stampare I * J
               stampare uno spazio*/
            printf("%3d ", I*J);
        /* andare a capo*/
        printf("\n");
    }
}
