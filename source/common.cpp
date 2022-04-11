#include ".././include/common.h"
#include <math.h>

/*Print the table containing the tree given in parameter*/
void PrintGHCutTree_Table(Tree *ghct)
{
	long i;

	for (i = 0; i < ghct->n_nodes; i++)
	{
		printf("%d  ", ghct->Tab[i].id);

		if (ghct->Tab[i].lv != NULL)
			printf("%d  ", ghct->Tab[i].lv->id);
		else
			printf("-  ");

		if (ghct->Tab[i].lh != NULL)
			printf("%d  ", ghct->Tab[i].lh->id);
		else
			printf("-  ");

		printf("%f\n", ghct->Tab[i].cap);
	}
}

void PrintGHCutTree(Tree *ghct)
{
	PrintTree(&(ghct->Tab[0]));
}

void PrintTree(Tree_Node *t)
{
	Tree_Node *t_son;

	if (t != NIL_TN)
	{
		printf("%ld ", t->id);

		t_son = t->lv;
		while (t_son != NIL_TN)
		{
			PrintTree(t_son);

			t_son = t_son->lh;
		}
	}
}

/*Delete the memory allocated for a b_Node set*/
void Delete_b_Node_Set(b_Node **f_w)
{
	b_Node *b_cour;

	while (*f_w != NIL_BN)
	{
		b_cour = (*f_w);
		*f_w = (*f_w)->next;
		free(b_cour);
		b_cour = NIL_BN;
	}
}

/*Delete the memory allocated for a b_Node set*/
void Delete_b_Edge_Set(b_Edge **f_w)
{
	b_Edge *b_cour;

	while (*f_w != NIL_BE)
	{
		b_cour = (*f_w);
		*f_w = (*f_w)->next;
		free(b_cour);
		b_cour = NIL_BE;
	}
}

/*Delete the memory allocated to the tree*/
void Delete_Tree(Tree **tr)
{
	if ((*tr) != NULL)
	{
		free((*tr)->Tab);
		(*tr)->Tab = NULL;

		Delete_b_Edge_Set(&((*tr)->b_List));

		free(*tr);
	}
}

/*delete the memory allocated to a global mincut*/
void Delete_Mincut(Mincut *gmcu)
{
	Delete_b_Node_Set(&gmcu->f_w);
	Delete_b_Node_Set(&gmcu->f_w_b);
	free(gmcu);
}

/*Supprime une LCC de la m�moire*/
void Delete_LCC_Node(b_Node **f_w)
{
	b_Node *b_cour;

	while ((*f_w)->next != (*f_w))
	{
		b_cour = (*f_w)->next;
		(*f_w)->next = b_cour->next;
		free(b_cour);
	}

	free((*f_w));
	(*f_w) = NIL_BN;
}

/*Print the set of b_Node W*/
void Print_b_Node_Set(b_Node *first_w, FILE *out_file)
{
	b_Node *n_cour;

	n_cour = first_w;
	while (n_cour != NIL_BN)
	{
		std::cout << n_cour->id << "  ";

		n_cour = n_cour->next;
	}
}

void Print_LCC_b_Node(b_Node *first_w, FILE *out_file)
{
	b_Node *n_cour;

	n_cour = first_w;
	do
	{
		std::cout << n_cour->id <<"   ";

		n_cour = n_cour->next;
	} while (n_cour->id != first_w->id); 
}

/*Print the b_Edge set delta(W)*/
void Print_b_Edge_Set(b_Edge *f_e, FILE *out_file)
{
	b_Edge *e_cour;

	e_cour = f_e;
	while (e_cour != NIL_BE)
	{
		std::cout << e_cour->node1 << "  " << e_cour->node2 << "   " <<  e_cour->cap <<  "   " << e_cour->num << "  "<< std::endl;

		e_cour = e_cour->next;
	}
}

/*This function renumbers the nodes of a given set of nodes W, in order to apply*/
/*the ghct function (for instance) on the subgraph induced by W (ie G[W]). You must*/
/*notice that the ghct function and the NOI function work on graph which nodes are*/
/*consecutives (ie numerotate from 1 to n).*/
/*f_w is the head of the linked list of nodes representing the set W.*/
/*Tab1_2 gives the the correspondant of a node of W in the new numerotation.*/
/*n_12 is the number of nodes of the graph G from which the set W is comming from.*/
/*Tab_2_1 allows to make the correspondance in the back way.*/
/*n_21 is the number of nodes in the set W.*/
void RenumberTables(b_Node *f_w, simpleNode **Tab_12, long n_12, simpleNode **Tab_21, long *n_21)
{
	b_Node *b_cour;
	simpleNode *aux;
	long k;

	/*Allocation of Tab_12 and Tab_21:the default value will be 0*/
	(*Tab_12) = (simpleNode *)calloc(n_12, sizeof(simpleNode));
	(*Tab_21) = (simpleNode *)malloc(n_12 * sizeof(simpleNode));

	b_cour = f_w;
	*n_21 = 0;
	while (b_cour != NIL_BN)
	{
		(*n_21)++;

		k = b_cour->id;
		(*Tab_12)[k - 1] = (*n_21);

		(*Tab_21)[(*n_21) - 1] = k;

		b_cour = b_cour->next;
	}

	aux = *Tab_21;
	*Tab_21 = (simpleNode *)malloc((*n_21) * sizeof(simpleNode));

	for (k = 0; k < (*n_21); k++)
	{
		(*Tab_21)[k] = aux[k];
	}

	free(aux);
}

/*Create a new list of b_Node from the list f_w. corresp_tab says what*/
/*is the correspondant of each node.If you are doing a old-to-new correspondance,*/
/*you should give the table Tab_12 obtained with the function RenumberTables().*/
/*If you want to make a new-to-old correspondance, you should give Tab_21 as parameter.*/
b_Node *renumber_nodes(b_Node *f_w, simpleNode *corresp_tab, long n_corresp)
{
	b_Node *b_cour, *f_new, **new_cour;
	long k;

	f_new = NIL_BN;
	new_cour = &f_new;

	b_cour = f_w;
	while (b_cour != NIL_BN)
	{
		k = b_cour->id;

		*new_cour = (b_Node *)malloc(sizeof(b_Node));
		(*new_cour)->id = corresp_tab[k - 1];
		(*new_cour)->next = NIL_BN;

		b_cour = b_cour->next;
		new_cour = &((*new_cour)->next);
	}

	return f_new;
}

/*Create a new list of simpleEdge from the list Tab. corresp_tab says what*/
/*is the correspondant of each node.If you are doing a old-to-new correspondance,*/
/*you should give the table Tab_12 obtained with the function RenumberTables().*/
/*If you want to make a new-to-old correspondance, you should give Tab_21 as parameter.*/
simpleEdge *renumber_edges(simpleEdge *Tab, long n_tab, simpleNode *corresp_tab, long n_corresp)
{
	long i, k;
	simpleEdge *tab_new;

	tab_new = (simpleEdge *)malloc(n_tab * sizeof(simpleEdge));

	for (i = 0; i < n_tab; i++)
	{
		k = Tab[i].node1;
		tab_new[i].node1 = corresp_tab[k - 1];

		k = Tab[i].node2;
		tab_new[i].node2 = corresp_tab[k - 1];

		tab_new[i].cap = Tab[i].cap;
	}

	return tab_new;
}

/*This function similar to renumber_nodes() and renumber_edges().*/
/*The new function achieve in one the same time the job done by the*/
/*renumber_nodes(), renumber_edges() and RenumberTables().*/
/*It returns the new list of edges obtained by renumbering.*/
simpleEdge *renumber_edge_list(simpleEdge *Te, long m, simpleNode **T1, long n1, simpleNode **T2, long *n2)
{
	long i, node1, node2, new_number1, new_number2;
	simpleEdge *Tab;

	Tab = (simpleEdge *)malloc(m * sizeof(simpleEdge));

	(*T1) = (simpleNode *)malloc(n1 * sizeof(simpleNode)); /*The table must be initialized with 0*/
	(*T2) = (simpleNode *)malloc(n1 * sizeof(simpleNode));
	(*n2) = 0;

	for (i = 0; i < n1; i++)
	{
		(*T1)[i] = 0;
		(*T2)[i] = 0;
	}

	for (i = 0; i < m; i++)
	{
		node1 = Te[i].node1;
		node2 = Te[i].node2;

		/*If the node has already been relabel then we take its new label.*/
		/*Else we insert a new label for the node in the new labels table T1*/
		/*and we insert the node in the old labels table T2.*/
		if ((*T1)[node1 - 1] != 0)
		{
			new_number1 = (*T1)[node1 - 1];
		}
		else
		{
			(*n2)++;
			(*T1)[node1 - 1] = (*n2);
			(*T2)[(*n2) - 1] = node1;

			new_number1 = (*n2);
		}

		/*Same as node1*/
		if ((*T1)[node2 - 1] != 0)
		{
			new_number2 = (*T1)[node2 - 1];
		}
		else
		{
			(*n2)++;
			(*T1)[node2 - 1] = (*n2);
			(*T2)[(*n2) - 1] = node2;

			new_number2 = (*n2);
		}

		/*We insert a new edge in the new edge table*/
		Tab[i].node1 = new_number1;
		Tab[i].node2 = new_number2;
		Tab[i].cap = Te[i].cap;
	}

	return Tab;
}

/*Inserts a b_Node in the given b_Node linked list.*/
/*The function inserts the b_Node such as the first*/
/*node of the list is the one with the lower id.*/
void insert_b_Node(b_Node **f_list, b_Node **b_cour)
{
	if (*f_list != NIL_BN)
	{
		if ((*f_list)->id > (*b_cour)->id)
		{
			(*b_cour)->next = (*f_list);
			(*f_list) = (*b_cour);
		}
		else
		{
			(*b_cour)->next = (*f_list)->next;
			(*f_list)->next = (*b_cour);
		}
	}
	else
	{
		(*b_cour)->next = (*f_list);
		(*f_list) = (*b_cour);
	}
}

double max(double a, double b)
{
	if (a > b)
		return a;
	else
		return b;
}

double min(double a, double b)
{
	if (a < b)
		return a;
	else
		return b;
}

/*Renvoie l'indice dans la liste de l'arete e = (node1,node2)*/
long search_edge(simpleEdge *eList, long n_List, long node1, long node2)
{
	long i;
	Bool trouver;

	trouver = kECSP_False;
	i = -1;
	do
	{
		i++;

		if (((eList[i].node1 == node1) && (eList[i].node2 == node2)) ||
			((eList[i].node1 == node2) && (eList[i].node2 == node1)))
		{
			trouver = kECSP_True;
		}
	} while ((i < n_List - 1) && (!trouver));

	if (trouver == kECSP_True)
	{
		return i;
	}
	else
		return -1;
}

/*Calcul l'intersection des cercles (x1,y1,R1) et (x2,y2,R2)*/
/*Renvoie:						    */
/* +   -1 si le calcul est impossible.*/
/* +   0 s'il n'y a pas d'intersection.*/
/* +   1 s'il y a un seul point.Dans ce cas le point est stocke dans (xa,ya).*/
/* +   2 s'il y a deux points d'intersection.*/
int intersection_cercle(double x1, double y1, double R1, double x2, double y2, double R2, double *xa, double *ya, double *xb, double *yb)
{
	double cN, cA, cB, cC, delta;

	/*Impossible de calculer avec cette formule*/
	if (y1 == y2)
		return -1;

	cN = (R2 * R2 - R1 * R1 - x2 * x2 + x1 * x1 - y2 * y2 + y1 * y1) / (2 * (y1 - y2));
	cA = pow((x1 - x2) / (y1 - y2), 2) + 1.0;
	cB = 2 * y1 * (x1 - x2) / (y1 - y2) - 2 * cN * (x1 - x2) / (y1 - y2) - 2 * x1;
	cC = x1 * x1 + y1 * y1 + cN * cN - R1 * R1 - 2 * y1 * cN;

	delta = cB * cB - 4 * cA * cC;

	/*Pas d'intersection pour ces deux cercles*/
	if (delta < 0)
		return 0;

	/*Dans ce cas il y a un seul point d'intersection*/
	if (delta == 0)
	{
		*xa = (-cB + sqrt(delta)) / (2 * cA);
		*ya = cN - (*xa) * (x1 - x2) / (y1 - y2);

		return 1;
	}

	/*Dans ce cas il y a deux points d'intersection*/
	if (delta > 0)
	{
		*xa = (-cB + sqrt(delta)) / (2 * cA);
		*ya = cN - (*xa) * (x1 - x2) / (y1 - y2);

		*xb = (-cB - sqrt(delta)) / (2 * cA);
		*yb = cN - (*xb) * (x1 - x2) / (y1 - y2);

		return 2;
	}

	return 0;
}




Bool pile_init(Pile **p, long sz)
{
	if (*p != NULL)
		return kECSP_False;

	//(*p) = (Pile *)malloc(sizeof(Pile));
	(*p) = new Pile;

	if (*p == NULL)
		return kECSP_False;

	//(*p)->Tab = (data_type *)malloc(sz*sizeof(data_type));
	(*p)->Tab = new data_type[sz];

	if ((*p)->Tab == NULL)
		free(p);

	(*p)->size = sz;
	(*p)->index = -1;

	return kECSP_True;
}

Bool pile_empiler(Pile *p, data_type w)
{
	if (p->index < p->size - 1)
	{
		p->index++;

		p->Tab[p->index] = w;

		return kECSP_True;
	}
	else
		return kECSP_False;
}

data_type pile_depiler(Pile *p)
{
	data_type val;

	if (p->index > -1)
	{
		val = p->Tab[p->index];
		p->index--;

		return val;
	}
	else
		return data_type_undefined_value;
}

Bool pile_vide(Pile *p)
{
	if (p->index > -1)
		return kECSP_False;
	else
		return kECSP_True;
}

long pile_get_size(Pile *p)
{
	return p->size;
}

Bool pile_delete_pile(Pile **p)
{
	if ((*p) != NULL)
	{
		// free((*p)->Tab);
		delete[](*p)->Tab;
		(*p)->Tab = NULL;

		free(*p);

		return kECSP_True;
	}

	return kECSP_False;
}

/*Functions for using the Fifo*/

Bool fifo_init(Fifo **fifo, long f_size)
{
	//(*fifo) = (Fifo *)malloc(sizeof(Fifo));
	(*fifo) = new Fifo;

	if ((*fifo) == NULL)
	{
		fprintf(stdout, "%s", "Unable to allocate Fifo");
		return kECSP_False;
	}

	(*fifo)->first = NIL_BN;
	(*fifo)->last = NIL_BN;
	(*fifo)->size = f_size;
	(*fifo)->nb_elem = 0L;

	return kECSP_True;
}

long fifo_insert(Fifo *fifo, long u)
{

	if (fifo->nb_elem < fifo->size)
	{
		if (fifo->nb_elem == 0L)
		{
			// fifo->first = (b_Node *)malloc(sizeof(b_Node));
			fifo->first = new b_Node;

			fifo->first->id = u;
			fifo->first->next = NIL_BN;

			fifo->last = fifo->first;
		}
		else
		{
			// fifo->last->next = (b_Node *)malloc(sizeof(b_Node));
			fifo->last->next = new b_Node;

			fifo->last = fifo->last->next;

			fifo->last->id = u;
			fifo->last->next = NIL_BN;
		}

		fifo->nb_elem++;

		return fifo->nb_elem;
	}
	else
		return -1;
}

long fifo_get_first(Fifo *fifo) /*Deletes the first element from the list*/
{
	long u;
	b_Node *cour;

	if (fifo->nb_elem > 0L)
	{
		cour = fifo->first;

		fifo->first = fifo->first->next;

		fifo->nb_elem--;

		u = cour->id;

		// free(cour);
		delete cour;

		if (fifo->nb_elem == 0)
			fifo->last = NIL_BN;
	}
	else
		u = -1L;

	return u;
}

void fifo_delete_fifo(Fifo **fifo)
{
	b_Node *cour;

	while ((*fifo)->first != NIL_BN)
	{
		cour = (*fifo)->first;
		(*fifo)->first = (*fifo)->first->next;

		// free(cour);
		delete cour;
	}

	(*fifo)->last = NIL_BN;

	(*fifo)->nb_elem = 0;
	(*fifo)->size = 0;

	delete (*fifo);
}

void fifo_print_fifo(Fifo *fifo)
{
	b_Node *cour;

	cour = fifo->first;
	while (cour != NIL_BN)
	{
		// fprintf(stdout,"%d ",cour->id);

		cour = cour->next;
	}
}
