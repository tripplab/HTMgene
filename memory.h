/* HTMGene Angélica Alejandra Serrano Rubio, Mauricio Carrillo-Tripp */
/* Biomolecular Diversity Lab, http://tripplab.com */

#ifndef MEMORY_H
#define MEMORY_H

/*Declaracion de bibliotecas*/
#include <stdio.h>
#include <stdlib.h>
#include "data.h"

/* Aloca memoria para una matriz por medio de un solo bloque*/
int **createMatrixBloque(int,int);

/* Aloca memoria para una matriz por medio de bloques */
int **createMatrixBloques(int,int);

/* Aloca memoria para un vector */
int *createVector(int);

/* Libera memoria para una matriz*/
void freeMatrix(int **);

/* Libera memoria para una matriz de varios bloques*/
void freeMatrixBloques(int **, int);

#endif // MEMORY_H
