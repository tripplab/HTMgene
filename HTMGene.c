/* HTMGene Angélica Alejandra Serrano Rubio, Mauricio Carrillo-Tripp */
/* Biomolecular Diversity Lab, http://tripplab.com */

/*Declaración de bibliotecas*/
#include "HTMGene.h"
#include "memory.h"
#include "math.h"

/*Llena estructura*/
s_edge *LlenarEstructura(s_edge *ptrLista, s_edge *ptr)
{
	s_edge *lp = ptrLista;
	/*Si la lista es distinto de null agregar un nuevo nodo*/
	if(ptrLista!=NULL)
    {
        while(ptrLista->siguiente!=NULL)
            ptrLista=ptrLista->siguiente;
        /*Reserva memoria dinamica*/
        ptrLista->siguiente=(struct edge*)malloc(sizeof(s_edge));
        /*Asignacion de parametros*/
        ptrLista=ptrLista->siguiente;
        ptrLista->siguiente=NULL;
        ptrLista->id=ptr->id;
        ptrLista->origen=ptr->origen;
        ptrLista->destino=ptr->destino;
        ptrLista->win=ptr->win; 
        ptrLista->probabilidad=ptr->probabilidad;
        ptrLista->score=ptr->score;        
        ptrLista->cols=(int*)malloc(ptrLista->score*sizeof(int));
        return lp;
    }
    /*Genera el nodo raiz de la estructura*/
    else
    {
        /*Reserva memoria dinamica*/
        ptrLista = (struct edge*)malloc(sizeof(s_edge));
        /*Asignacion de parámetros*/
        ptrLista->siguiente=NULL;
        ptrLista->id=ptr->id;
        ptrLista->origen=ptr->origen;
        ptrLista->destino=ptr->destino;
        ptrLista->win=ptr->win;  
        ptrLista->probabilidad=ptr->probabilidad;   
        ptrLista->score=ptr->score;             
        ptrLista->cols=(int*)malloc(ptrLista->score*sizeof(int));
        return ptrLista;
    }
}

/*Imprime los datos de la estructura*/
void ImprimirEstructura(s_edge *ptrLista)
{
    /*Apuntador al inicio de la lista*/
    s_edge *ptrInit=NULL;
    int *ptrInitc=NULL;
    int i=0;
    ptrInit=ptrLista;
    if(ptrLista==NULL)
        printf("Error -- No  data-- \n");
    else
        /*Mientras la lista no sea diferente de nulo*/
        while(ptrLista != NULL)
        {
        	if(ptrLista->win)
        	{
	            printf("%d: origen %d destino %d score %d %d\n", ptrLista->id, ptrLista->origen, ptrLista->destino, ptrLista->score, ptrLista->win);
	            ptrInitc=ptrLista->cols;
	            if(ptrLista->score>0 && ptrLista->cols!=NULL)
	            {
	            	printf("Columas\n");
	            	for(i=0; i<ptrLista->score; i++)
	            	{	
	            		printf("%d ",ptrLista->cols[i]);
	            	}
	            	printf("\n");
	            }
	            ptrLista->cols=ptrInitc;
        	}
            ptrLista=ptrLista->siguiente;
        }
    ptrLista=ptrInit;
}

/*Borra un nodo de la estructura*/
s_edge *BorrarDato(s_edge *ptrLista)
{
    s_edge *temp;
    temp=ptrLista->siguiente;
    free(ptrLista->cols);
    free(ptrLista);
    return temp;
}

/*Funcion recursiva para Limpiar la estructura*/
void LimpiaDatos(s_edge *ptrLista)
{
    /*Si la lista no es nula hace la llamada recursiva*/
    if(ptrLista != NULL)
        LimpiaDatos(BorrarDato(ptrLista));
}

/*Borra un nodo de la estructura*/
s_group *BorrarBic(s_group *ptrGroup)
{
    s_group *temp;
    temp=ptrGroup->siguiente;
    free(ptrGroup->rows);    
    free(ptrGroup->cols);
    free(ptrGroup);
    return temp;
}

/*Funcion recursiva para Limpiar la estructura*/
void LimpiaBic(s_group *ptrGroup)
{
    /*Si la lista no es nula hace la llamada recursiva*/
    if(ptrGroup != NULL)
        LimpiaBic(BorrarBic(ptrGroup));
}


/*Calcula la interseccion de dos filas*/
int IntersectCols(int **matrix, int cols, int row_a, int row_b)
{
	int i=0;
	int counter=0;
	for(i=0; i<cols; i++)
	{
		if(matrix[row_a][i]==matrix[row_b][i] && matrix[row_a][i]!=0)
		{
			counter++;
		}
	}
	return counter;
}

/*Intercambia los datos de la estructura*/
void Swap(s_edge **primero)
{
    struct edge *aux = NULL;
    struct edge *edge;
    /* Por cada nodo en la lista original apuntada por "*primero" */
    while (*primero != NULL)
    {
        /* Extraemos un nodo de la lista original */
        edge = *primero;
        *primero = (*primero)->siguiente;
        /*y lo insertamos (en orden) en la lista auxiliar */
        insertar(&aux, edge);
    }
    /* Ajustamos el puntero a la lista original */
    *primero = aux;
}

/*Funcion que inserta un nodo de acuerdo al ordenamiento del número aleatorio de menor a mayor*/
void insertar(struct edge **primero, struct edge *edge)
{
    /*Apuntador a estructura*/
    struct edge **p;
    p = primero;
    /*Mientras sea diferente de nulo y sea menor que el aleatorio siguiente*/
    while (*p != NULL && (*p)->score > edge->score)
        p = &((*p)->siguiente);
    edge->siguiente = *p;
    /*Apuntador a nodo*/
    *p = edge;
}

/*Calculo de la media aritmética*/
double media(s_edge *ptrLista)
{
	int suma=0;
	int counter=0;
	/*Apuntador al inicio de la lista*/
    s_edge *ptrInit=NULL;
    ptrInit=ptrLista;
    if(ptrLista==NULL)
        printf("Error -- No  data-- \n");
    else
        /*Mientras la lista no sea diferente de nulo*/
        while(ptrLista != NULL)
        {
            suma+=ptrLista->score;
            counter++;
            ptrLista=ptrLista->siguiente;
        }
    ptrLista=ptrInit;
    return(suma/counter);
}

/*Calculo de la desviacion estandar poblacional*/
double desviacion(s_edge *ptrLista)
{
	int suma=0;
	int counter=0;
	/*Apuntador al inicio de la lista*/
    s_edge *ptrInit=NULL;
    ptrInit=ptrLista;
    if(ptrLista==NULL)
        printf("Error -- No  data-- \n");
    else
    {
    	/*Calculo de la media aritmética*/
    	double mu=media(ptrLista);
        /*Mientras la lista no sea diferente de nulo*/
        while(ptrLista != NULL)
        {
            suma+=(ptrLista->score-mu)*(ptrLista->score-mu);
            counter++;
            ptrLista=ptrLista->siguiente;
        }
    }
    ptrLista=ptrInit;
    double desv=sqrt(suma/counter);
    return(desv);
}

/*InferenceSP*/
s_edge *InferenceSP(s_edge *ptrLista)
{
	double suma=0;
	/*Apuntador al inicio de la lista*/
    s_edge *ptrInit=NULL;
    ptrInit=ptrLista;
    if(ptrLista==NULL)
        printf("Error -- No  data-- \n");
    else
    {
    	/*Calculo de la desviacion estandar*/
    	double sigma=desviacion(ptrLista);
        /*Mientras la lista no sea diferente de nulo*/
        while(ptrLista != NULL)
        {
        	double tmp1=1/(sqrt(2*pi*sigma*sigma));
        	double tmp2=exp(2*sigma*sigma);
        	double tmp3=ptrLista->score/tmp2;
            ptrLista->probabilidad=tmp1*tmp3;
            suma+=ptrLista->probabilidad;
            ptrLista=ptrLista->siguiente;
        }
        ptrLista=ptrInit;
        /*Normalizacion*/
        while(ptrLista != NULL)
        {
            ptrLista->probabilidad/=suma;
            ptrLista=ptrLista->siguiente;
        }
    }
    ptrLista=ptrInit;
    return ptrLista;
}

/*Learning*/
s_edge *Preprocessing(s_edge *ptrLista, double umbral, int type)
{
	/*Apuntador al inicio de la lista*/
    s_edge *ptrInit=NULL;
    ptrInit=ptrLista;
    if(ptrLista==NULL)
        printf("Error -- No  data-- \n");
    else
    	/*Mientras la lista no sea diferente de nulo*/
        while(ptrLista != NULL)
        {
        	if(ptrLista->score >= umbral && type==1)
        		ptrLista->win=1;	
        	else if(ptrLista->probabilidad >= umbral && type==2)
	        	ptrLista->win=1;	
        	ptrLista=ptrLista->siguiente;
        }
    ptrLista=ptrInit;
    return ptrLista;
}

/*Lee las columnas con valores activos*/
s_edge *LearningSP(int **matrix, int cols, s_edge *ptrLista)
{
	/*Apuntador al inicio de la lista*/
    s_edge *ptrInit=NULL;
    int i=0;
    ptrInit=ptrLista;
    if(ptrLista==NULL)
        printf("Error -- No  data-- \n");
    else
    	/*Mientras la lista no sea diferente de nulo*/
        while(ptrLista != NULL)
        {
           int k=0;
	       for(i=0; i<cols; i++)
	       {
	          if((matrix[ptrLista->origen][i]==1 && matrix[ptrLista->destino][i]==1) ||
	           (matrix[ptrLista->origen][i]==-1 && matrix[ptrLista->destino][i]==-1) ||
	           (matrix[ptrLista->origen][i]==1 && matrix[ptrLista->destino][i]==-1) ||
	           (matrix[ptrLista->origen][i]==-1 && matrix[ptrLista->destino][i]==1) 
	           )
	          {
	        		ptrLista->cols[k]=i;
	        		k++;
	          }
	          if(k==ptrLista->score)
	        	break;
	        }
        	ptrLista=ptrLista->siguiente;
        }
    ptrLista=ptrInit;
    return ptrLista;
}

/*Overlap SP*/
s_edge *OverlapSP(int **matrix, s_edge *ptrLista, int rows, int cols)
{
	int i=0; 
	int j=0;
	int counter=0;
	for(i=0; i<rows; i++)
	{
		for(j=i+1;j<rows; j++)
		{
			s_edge *New=(struct edge*)malloc(sizeof(s_edge));
			New->siguiente=NULL;
			New->id=counter++;
			New->origen=i;
			New->destino=j;
			New->score=IntersectCols(matrix, cols, i, j);
			New->probabilidad=0.0;
			New->win=0;	
			New->cols=NULL;
			ptrLista=LlenarEstructura(ptrLista, New);
			free(New);
		}
	}  
	return ptrLista;
}

/*Llena estructura*/
s_group *LlenarBicluster(int **matrix, s_group *ptrGroup, s_group *ptr, int noRows, int noCols)
{
	s_group *lp = ptrGroup;
	/*Si la lista es distinto de null agregar un nuevo nodo*/
	if(ptrGroup!=NULL)
    {
        while(ptrGroup->siguiente!=NULL)
            ptrGroup=ptrGroup->siguiente;
        ptrGroup->siguiente=(struct group*)malloc(sizeof(s_group));
        ptrGroup=ptrGroup->siguiente;
        ptrGroup->siguiente=NULL;
        ptrGroup->id=ptr->id;
	   ptrGroup->condiciones=ptr->condiciones;
	   ptrGroup->genes=ptr->genes;	        
        ptrGroup->cols=ptr->cols;
        ptrGroup->pvalue=ptr->pvalue;
        return lp;
    }
    /*Genera el nodo raiz de la estructura*/
    else
    {
        /*Reserva memoria dinamica*/
        ptrGroup=(struct group*)malloc(sizeof(s_group));
        ptrGroup->siguiente=NULL;
        ptrGroup->id=ptr->id;
	    ptrGroup->condiciones=ptr->condiciones;
	    ptrGroup->genes=ptr->genes;        
        ptrGroup->cols=ptr->cols;
        ptrGroup->pvalue=ptr->pvalue;
        return ptrGroup;
    }
}

/*Llena estructura*/
s_nodo *LlenarGenes(s_nodo *ptrNodo, s_nodo *ptr)
{
	s_nodo *lp = ptrNodo;
	/*Si la lista es distinto de null agregar un nuevo nodo*/
	if(ptrNodo!=NULL)
    {
        while(ptrNodo->siguiente!=NULL)
        {
        	 if(ptrNodo->pos==ptr->pos)
        	 	return lp;
            ptrNodo=ptrNodo->siguiente;
        }
        if(ptrNodo->pos==ptr->pos)
        	 	return lp;
        ptrNodo->siguiente=(struct  nodo*)malloc(sizeof(s_nodo));
        ptrNodo=ptrNodo->siguiente;
        ptrNodo->siguiente=NULL;
        ptrNodo->pos=ptr->pos;
	     return lp;
    }
    /*Genera el nodo raiz de la estructura*/
    else
    {
    	/*Reserva memoria dinamica*/
        ptrNodo=(struct nodo*)malloc(sizeof(s_nodo));
        ptrNodo->siguiente=NULL;
        ptrNodo->pos=ptr->pos;
        return ptrNodo;
    }
}

/*Compara genes con las mismas condiciones*/
int compare(int **matrix, int *cond, int score, int origen, int destino, int Rows, int Cols)
{
	int i=0;
	int j=0;
	int counter=2;
	int band=0;
	for(i=0; i<Rows; i++)
	{
		band=0;
		if(i!=origen && i!=destino)
		{
			for(j=0; j<score; j++)
			{
				if(!matrix[i][cond[j]])
				{
					band=1;
				}
			}
			if(!band)
				counter++;
		}
	}
	return counter;
}

/*Conecta los genes con las mismas condiciones*/
s_group *connectivity(s_edge *ptrLista, s_group *ptrGroup)
{
	/*Apuntador al inicio de la lista*/
    s_edge *ptrInit=NULL;
    s_group *ptrInitg=NULL;
    ptrInitg=ptrGroup;
    int i=0;
    int band=0;
    if(ptrGroup==NULL)
        printf("Error -- No  data-- \n");
    else
    	/*Mientras la lista no sea diferente de nulo*/
        while(ptrGroup != NULL)
        {
	        	ptrInit=ptrLista;
	        	ptrGroup->rows=NULL;
	        	while(ptrLista != NULL)
	        	{
	        		if(ptrLista->win)
        			{
		        		band=0;
		        		if(ptrGroup->condiciones==ptrLista->score)
		        		{
		        			for(i=0; i<ptrLista->score; i++)
				       		if(ptrGroup->cols[i]==ptrLista->cols[i])
				       			band++;
					       if(band==ptrLista->score)
				       		band=1;
					       else
					        	band=0;
					       if(band)
				        	{
				        		/*insertar origen*/
				        		s_nodo *New=(struct nodo*)malloc(sizeof(s_nodo));
								New->siguiente=NULL;
								New->pos=ptrLista->origen;
								ptrGroup->rows=LlenarGenes(ptrGroup->rows, New);
								free(New);
								/*insertar destino*/
								New=(struct nodo*)malloc(sizeof(s_nodo));
								New->siguiente=NULL;
								New->pos=ptrLista->destino;
								ptrGroup->rows=LlenarGenes(ptrGroup->rows, New);
								free(New);
				        	}
		        		}
			       }
			       ptrLista=ptrLista->siguiente;
	        	}
	        	ptrLista=ptrInit;
	        	ptrGroup=ptrGroup->siguiente;
        }
    ptrLista=ptrInit;
    ptrGroup=ptrInitg;
    return ptrGroup;
}


/*Agrupar por condiciones*/
s_group *GroupingTP(int **matrix, s_group *ptrGroup, s_edge *ptrLista, int noRows, int noCols)
{
	/*Apuntador al inicio de la lista*/
    s_edge *ptrInit_e=ptrLista;
    int counter=0;
    int i=0;
    s_group *ptr=NULL;
    
    if(ptrLista==NULL)
        printf("Error -- No  data-- \n");
    else
    	/*Mientras la lista no sea diferente de nulo*/
        while(ptrLista != NULL)
        {
        	if(ptrLista->win)
        	{
	        	/*Primer elemento se almacena directamente*/
	        	if(!counter)
	        	{
	        		s_group *New=(struct group*)malloc(sizeof(s_group));
					New->siguiente=NULL;
					New->id=counter++;
					New->condiciones=ptrLista->score;
					/*Almacenar las condiciones*/	  
					New->cols=(int *)malloc(ptrLista->score*sizeof(int));
					for(i=0; i<ptrLista->score; i++)
						New->cols[i]=ptrLista->cols[i];
					New->genes=0;						  
					New->rows=NULL;
					New->pvalue=0.0;
					ptrGroup=LlenarBicluster(matrix, ptrGroup, New, noRows, noCols);
					free(New);
		          ptr=ptrGroup;
	        	}
	        	else
	        	{
	        		/*Compara elementos*/
	        		s_group *tmp=ptr;
	        		int band=0;
	        		while(ptr!=NULL)
	        		{
	        			band=0;
		        		for(i=0; i<ptrLista->score; i++)
		        		{
		        			if(ptr->cols[i]==ptrLista->cols[i])
								band++;
		        		}
		        		if(band==ptrLista->score && ptr->condiciones==ptrLista->score)
		        		{
		        			band=1;
		        			break;
		        		}
		        		else
		        			band=0;
		        		ptr=ptr->siguiente;
	        		}
	        		ptr=tmp;
	        		if(!band)
	        		{
	        			s_group *New=(struct group*)malloc(sizeof(s_group));
						New->siguiente=NULL;
						New->id=counter++;
						New->condiciones=ptrLista->score;
						/*Almacenar las condiciones*/	  
						New->cols=(int *)malloc(ptrLista->score*sizeof(int));
						for(i=0; i<ptrLista->score; i++)
							New->cols[i]=ptrLista->cols[i];
						New->genes=0;						  
						New->rows=NULL;
						New->pvalue=0.0;
						ptrGroup=LlenarBicluster(matrix, ptrGroup, New, noRows, noCols);
						free(New);
	        		}
	        	}
        	}
        	ptrLista=ptrLista->siguiente;
        }
    ptrLista=ptrInit_e;  
    return ptrGroup;
}

/*Imprime los datos de la estructura*/
void ImprimirBicluster(s_group *ptrGroup, int ** matrix)
{
    /*Apuntador al inicio de la lista*/
    s_group *ptrInit=NULL;
    int i=0;
    int j=0;
    int counter=0;
    s_nodo *ptr;
    ptrInit=ptrGroup;
    if(ptrGroup==NULL)
        printf("Error -- No  data-- \n");
    else
        /*Mientras la lista no sea diferente de nulo*/
        while(ptrGroup != NULL)
        {
        	  counter=0;
        	  printf("Bicluster [%d]\n", ptrGroup->id);
        	  printf("Conditions: ");
        	  for(i=0; i<ptrGroup->condiciones; i++)
        	  	printf("%d ", ptrGroup->cols[i]);
        	  printf("\n"); 	
        	  printf("Genes: ");
        	  ptr=ptrGroup->rows;
        	  while(ptrGroup->rows != NULL)
        	  {
        	  	  printf("%d ", ptrGroup->rows->pos);
        	  	  ptrGroup->rows=ptrGroup->rows->siguiente;
        	  	  counter++;
        	  }
        	  ptrGroup->rows=ptr;
        	  ptrGroup->genes=counter;
        	  printf("\n"); 	
        	  printf("Matrix: \n  \t");
        	  for(i=0; i<ptrGroup->condiciones; i++)
					printf("c_%d\t", ptrGroup->cols[i]);
			  printf("\n");
        	  for(i=0; i<ptrGroup->genes; i++)
        	  {
        	  	 printf("g_%d\t", ptrGroup->rows->pos);
        	  	 for(j=0; j<ptrGroup->condiciones; j++)
        	  	 {
        	  	 	 printf("%d\t", matrix[ptrGroup->rows->pos][ptrGroup->cols[j]]);
        	  	 }
        	  	 printf("\n");
        	  	 ptrGroup->rows=ptrGroup->rows->siguiente;
        	  }
        	  printf("\n");
           ptrGroup=ptrGroup->siguiente;
        }
    ptrGroup=ptrInit;
}

/*Imprime los datos de la estructura*/
void ImprimirBiclusterFile(s_group *ptrGroup, FILE *output, int **matrix)
{
        /*Apuntador al inicio de la lista*/
    s_group *ptrInit=NULL;
    int i=0;
    int j=0;
    int counter=0;
    s_nodo *ptr;
    ptrInit=ptrGroup;
    if(ptrGroup==NULL)
        fprintf(output, "Error -- No  data-- \n");
    else
        /*Mientras la lista no sea diferente de nulo*/
        while(ptrGroup != NULL)
        {
        	  counter=0;
        	  fprintf(output, "Bicluster [%d]\n", ptrGroup->id);
        	  fprintf(output, "Conditions: ");
        	  for(i=0; i<ptrGroup->condiciones; i++)
        	  	fprintf(output, "%d ", ptrGroup->cols[i]);
        	  fprintf(output, "\n"); 	
        	  fprintf(output, "Genes: ");
        	  ptr=ptrGroup->rows;
        	  while(ptrGroup->rows != NULL)
        	  {
        	  	  fprintf(output, "%d ", ptrGroup->rows->pos);
        	  	  ptrGroup->rows=ptrGroup->rows->siguiente;
        	  	  counter++;
        	  }
        	  ptrGroup->rows=ptr;
        	  ptrGroup->genes=counter;
        	  fprintf(output, "\n"); 	
        	  fprintf(output, "Matrix: \n  \t");
        	  for(i=0; i<ptrGroup->condiciones; i++)
					fprintf(output, "c_%d\t", ptrGroup->cols[i]);
			  fprintf(output, "\n");
        	  for(i=0; i<ptrGroup->genes; i++)
        	  {
        	  	 fprintf(output, "g_%d\t", ptrGroup->rows->pos);
        	  	 for(j=0; j<ptrGroup->condiciones; j++)
        	  	 {
        	  	 	 fprintf(output, "%d\t", matrix[ptrGroup->rows->pos][ptrGroup->cols[j]]);
        	  	 }
        	  	 fprintf(output, "\n");
        	  	 ptrGroup->rows=ptrGroup->rows->siguiente;
        	  }
        	  fprintf(output, "\n");
           ptrGroup=ptrGroup->siguiente;
        }
    ptrGroup=ptrInit;
}

/*Algoritmo HTMGene*/
int HTMGene(int **matrix, int noCols, int noRows, double minCols, int minRows, char *FileName, int type, int Umbral)
{
	/*Apuntador a agrupacion SP*/
	s_edge *ptrLista = NULL;
	/*Spatial Pooler*/
	/*Overlap SP*/
	ptrLista=OverlapSP(matrix, ptrLista, noRows, noCols);
	/*Intercambia las posiciones de la lista*/
    Swap(&ptrLista);
	/*Inference SP*/
	ptrLista=InferenceSP(ptrLista);
	/*Preprocessing SP*/
	ptrLista=Preprocessing(ptrLista, minCols, Umbral);
    /*Lee las columnas con valores activos*/
    ptrLista=LearningSP(matrix, noCols, ptrLista);
    /*Apuntador a agrupaciones TP*/
    s_group *ptrGroup = NULL;
    /*Temporal Pooler*/
    /*Grouping TP*/
    s_edge *ptrInit=ptrLista;
    ptrGroup=GroupingTP(matrix, ptrGroup, ptrLista, noRows, noCols);
    ptrLista=ptrInit;
    /*Conecta los genes con las mismas condiciones*/
    ptrInit=ptrLista;
	ptrGroup=connectivity(ptrLista, ptrGroup);
    ptrLista=ptrInit;	 
    /*Imprime estructura*/
    if(type==1)
	    ImprimirBicluster(ptrGroup, matrix); 
	 else
	 {
	    /*Archivo que escribe la salida*/
	    FILE *output=NULL;
	    output=fopen(FileName, "w");
	    /*Imprime estructura*/
	    ImprimirBiclusterFile(ptrGroup, output, matrix);     
	    /*Cierra el archivo*/
	    fclose(output);  
	 }
	 /*Libera memoria*/	
    LimpiaDatos(ptrLista);
    return 0;
}
