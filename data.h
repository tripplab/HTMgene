/*HTMGene  Angélica Alejandra Serrano Rubio*/

#ifndef DATA_H
#define DATA_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*SP*/
typedef struct edge
{
	int id;
	int origen;
	int destino;
	int score;
	int win;
	int *cols;
	double probabilidad;
	struct edge *siguiente; 
}s_edge;

/*Almacena los nodos: filas o columnas*/
typedef struct nodo
{
	int pos;
	struct nodo *siguiente;
}s_nodo;

/*Almacena un bicluster*/
typedef struct group
{
	int id;
	int condiciones;
	int genes;
	int *cols;
	s_nodo *rows;
	float pvalue;
	struct group *siguiente;
}s_group;


#endif // DATA_H
