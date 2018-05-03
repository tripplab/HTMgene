/* HTMGene Angélica Alejandra Serrano Rubio, Mauricio Carrillo-Tripp */
/* Biomolecular Diversity Lab, http://tripplab.com */

/*Declaración de bibliotecas*/
#include "memory.h"

/*Almacena memoria en un solo bloque*/
int **createMatrixBloque(int nr, int nc)
{
    int i;
    int **Matrix=NULL;
    //Reservar memoria
    Matrix=(int **) malloc((nr)*sizeof(int *));
    /*Si no se pudo reservar memoria*/
    if(Matrix==NULL)
        return(NULL);
    /*Almacena memoria en un bloque*/
    Matrix[0]=(int *)malloc(nr*nc*sizeof(int));
    if(Matrix[0]==NULL)
        return(NULL);
    /*Genera un apuntador a cada inicio de fila*/
    for(i=1; i<=nr; i++)
        Matrix[i]=Matrix[i-1]+nc;
    return(Matrix);
}

/*Almacena memoria en varios blaques*/
int **createMatrixBloques(int nr, int nc)
{
    int i;
    int **Matrix=NULL;
    // Reservamos memoria
    Matrix = (int **) malloc( nr*sizeof(int *));
    if(Matrix==NULL)
        return(NULL);
    /*Genera distintos bloques de memoria*/
    for(i=0; i<nr; ++i)
        Matrix[i] = (int *) malloc( nc*sizeof(int));
    return(Matrix);
}

/*Almacena memoria para un vector*/
int *createVector(int n)
{
    int *Vector=NULL;
    /*Aloca memoria para un vector de n elementos de tipo int*/
    Vector=(int*)malloc(n*sizeof(int));
    if(Vector==NULL)
        return(NULL);
    return(Vector);
}

/* Libera memoria para una matriz de un solo bloque*/
void freeMatrix(int **mat)
{
    free(mat[0]);
    free(mat);
}

/* Libera memoria para una matriz de varios bloques*/
void freeMatrixBloques(int **mat, int nr)
{
    int i=0;
    for(i=0; i<nr; i++)
        free(mat[i]);
    free(mat);
}
