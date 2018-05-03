/* HTMGene Angélica Alejandra Serrano Rubio, Mauricio Carrillo-Tripp */
/* Biomolecular Diversity Lab, http://tripplab.com */

#ifndef HTMGENE_H
#define HTMGENE_H

/*Declaración de bibliotecas*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "data.h"
#define pi 3.1416

/*Algoritmo HTMGene*/
int HTMGene(int **,int, int, double, int, char *, int, int);
/*Overlap SP*/
s_edge *OverlapSP(int **, s_edge *, int, int);
/*Calcula la interseccion de dos filas*/
int IntersectCols(int **, int, int, int);
/*Llena estructura*/
s_edge *LlenarEstructura(s_edge *, s_edge *);
/*Imprime los datos de la estructura*/
void ImprimirEstructura(s_edge *);
/*Funcion recursiva para Limpiar la estructura*/
void LimpiaDatos(s_edge *);
/*Borra un nodo de la estructura*/
s_edge *BorrarDato(s_edge *);
/*Intercambia los datos de la estructura*/
void Swap(s_edge **);
/*Funcion que inserta un nodo de acuerdo al ordenamiento del número aleatorio de menor a mayor*/
void insertar(struct edge **, struct edge *);
/*Calculo de la media aritmética*/
double media(s_edge *);
/*Calculo de la desviacion estandar poblacional*/
double desviacion(s_edge *);
/*InferenceSP*/
s_edge *InferenceSP(s_edge *);
/*Learning*/
s_edge *Preprocessing(s_edge *, double, int);
/*Lee las columnas con valores activos*/
s_edge *LearningSP(int **, int, s_edge*);
/*Agrupar por condiciones*/
s_group *GroupingTP(int **, s_group *, s_edge *, int, int);
/*Llena estructura*/
s_group *LlenarBicluster(int **, s_group *, s_group *, int, int);
/*Borra un nodo de la estructura*/
s_group *BorrarBic(s_group *);
/*Funcion recursiva para Limpiar la estructura*/
void LimpiaBic(s_group *);
/*Imprime lista*/
void ImprimirBicluster(s_group *, int **);
/*Conecta los genes con las mismas condiciones*/
s_group *connectivity(s_edge *, s_group *);
/*Llena estructura*/
s_nodo *LlenarGenes(s_nodo *, s_nodo *);
/*Imprime los datos de la estructura*/
void ImprimirBiclusterFile(s_group *, FILE *, int **);

#endif  //HTMGENE_H
