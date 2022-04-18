#ifndef COMMON_H
#define COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#define kECSP_True 1
#define kECSP_False 0

#define EPS 1.0E-10

using namespace std;

typedef enum
{
	FALSE = 0,
	TRUE = 1
} BOOL;

typedef short Bool;

/*This representation allows to create linked list of nodes*/
typedef struct b_Node
{
	long id;
	struct b_Node *next; /*Next in the linked list*/
} b_Node;

typedef long simpleNode; /*For contiguous list of nodes*/

/*Simple representation of an edge*/
typedef struct simpleEdge
{
	long node1;
	long node2;
	double cap;
} simpleEdge;

/*Representation of edges used to create*/
/*linked list of edges.*/
typedef struct b_Edge
{
	long node1;
	long node2;
	double cap;
	long num;
	struct b_Edge *next;
} b_Edge;

/*A Tree representation with horizontal link(for brothers) and*/
/*vertical link(for sons) Tree_Node represents a node of the tree*/
typedef struct Tree_Node
{
	long id;			  /*Id of the node*/
	struct Tree_Node *lv; /*Vertical link of this node*/
	struct Tree_Node *lh; /*Horizontal link of this node*/
	double cap;			  /*capacity of the edge between the*/
						  /*node	and its unique parent*/
} Tree_Node;

typedef struct Tree
{
	long n_nodes;	/*Number of nodes of the tree*/
	Tree_Node *Tab; /*Tree nodes table */
	b_Edge *b_List; /*Linked list of the edges of the tree that*/
					/*don't statisfied some criterion*/
} Tree;

/*Data structure which represents a global mincut of a graph*/
typedef struct Mincut
{
	b_Node *f_w;
	b_Node *f_w_b;
	double mincap;
} Mincut;

#define NIL_BN (b_Node *)0
#define NIL_BE (b_Edge *)0
#define NIL_SE (simpleEdge *)0
#define NIL_TN (Tree_Node *)0

#define W_ONLY 1
#define W_B_ONLY 2
#define W_AND_W_B 3

std::string GetFileNameFromPath(const char* _buffer);
/*Print the table containing the tree given in parameter*/
void PrintGHCutTree_Table(Tree *ghct);

void PrintTree(Tree_Node *t);

/*Print the tree data structure given in parameter*/
void PrintGHCutTree(Tree *ghct);

/*Free the memory allocated to the tree*/
void Delete_Tree(Tree **tr);

/*delete the memory allocated to a global mincut*/
void Delete_Mincut(Mincut *gmcu);

/*Free the memory allocated for a b_Node set*/
void Delete_b_Node_Set(b_Node **f_w);

/*Free the memory allocated for a b_Edge set*/
void Delete_b_Edge_Set(b_Edge **f_w);

/*Print the set of b_Node W*/
void Print_b_Node_Set(b_Node *first_w, FILE *out_file);

/*Print the b_Edge set delta(W)*/
void Print_b_Edge_Set(b_Edge *f_e, FILE *out_file);

/*This function renumbers the nodes of a given set of nodes W, in order to apply*/
/*the ghct function (for instance) on the subgraph induced by W (ie G[W]). You must*/
/*notice that the ghct function and the NOI function work on graph which nodes are*/
/*consecutives (ie numerotate from 1 to n).*/
/*f_w is the head of the linked list of nodes representing the set W.*/
/*Tab_12 gives the correspondant of a node of W in the new numerotation(old -> new).*/
/*n_12 is the number of nodes of the graph G from which the set W is coming from.*/
/*Tab_21 allows to make the correspondance in the back way(new -> old).*/
/*n_21 is the number of nodes in the set W.*/
void RenumberTables(b_Node *f_w, simpleNode **Tab_12, long n_12, simpleNode **Tab_21, long *n_21);

/*Create a new list of b_Node from the list f_w. corresp_tab says what*/
/*is the correspondant of each node.If you are doing a old-to-new correspondance,*/
/*you should give the table Tab_12 obtained with the function RenumberTables().*/
/*If you want to make a new-to-old correspondance, you should give Tab_21 as parameter.*/
b_Node *renumber_nodes(b_Node *f_w, simpleNode *corresp_tab, long n_corresp);

/*Create a new list of simpleEdge from the list Tab. corresp_tab says what*/
/*is the correspondant of each node.If you are doing a old-to-new correspondance,*/
/*you should give the table Tab_12 obtained with the function RenumberTables().*/
/*If you want to make a new-to-old correspondance, you should give Tab_21 as parameter.*/
simpleEdge *renumber_edges(simpleEdge *Tab, long n_tab, simpleNode *corresp_tab, long n_corresp);

/*This function similar to renumber_nodes() and renumber_edges().*/
/*The new function achieve in one the same time the job done by the*/
/*renumber_nodes(), renumber_edges() and RenumberTables().*/
/*It returns the new list of edges obtained by renumbering.*/
simpleEdge *renumber_edge_list(simpleEdge *Te, long m, simpleNode **T1, long n1, simpleNode **T2, long *n2);

/*Inserts a b_Node in the given b_Node linked list.*/
/*The function inserts the b_Node such as the first*/
/*node of the list is the one with the lower id.*/
void insert_b_Node(b_Node **f_list, b_Node **b_cour);

double min(double a, double b);

double max(double a, double b);

/*Renvoie l'indice dans la liste de l'arete e = (node1,node2)*/
long search_edge(simpleEdge *eList, long n_List, long node1, long node2);

/*Calcul l'intersection des cercles (x1,y1,R1) et (x2,y2,R2)*/
/*Renvoie:						    */
/* +   -1 si le calcul est impossible.*/
/* +   0 s'il n'y a pas d'intersection.*/
/* +   1 s'il y a un seul point.Dans ce cas le point est stocke dans (xa,ya).*/
/* +   2 s'il y a deux points d'intersection.*/
int intersection_cercle(double x1, double y1, double R1, double x2, double y2, double R2, double *xa, double *ya, double *xb, double *yb);

/*Supprime une LCC*/
void Delete_LCC_Node(b_Node **list);

void Print_LCC_b_Node(b_Node *first_w, FILE *out_file);




/*Type of the data stored in the stack*/
/*If you want to store int for example*/
/*you can specify it in the following typedef.*/
/*You should also redefine the data_type_undefined_value*/
typedef b_Node *data_type;

#define data_type_undefined_value (b_Node *)0

typedef struct Pile
{
  data_type *Tab;
  long index; /*index of the head of the stack*/
  long size;  /*Size of the stack*/
} Pile;

Bool pile_init(Pile **p, long sz);

Bool pile_empiler(Pile *p, data_type w);

data_type pile_depiler(Pile *p);

long pile_get_size(Pile *p);

Bool pile_delete_pile(Pile **p);

Bool pile_vide(Pile *p); /*Return True if the stack is empty*/

/*Data structure for representing a queue*/
/*The data that will be stored in the queue*/
/*will be pointers to Node.*/
typedef struct Fifo
{
  b_Node *first; /*Pointer to the head of the list*/
  b_Node *last;  /*Pointer to the end of the list*/
  long size;     /*Size of the Fifo*/
  long nb_elem;

} Fifo;

/*Functions for using the Fifo*/

long fifo_get_first(Fifo *fifo); /*Deletes the last element from the table*/

long fifo_insert(Fifo *fifo, long u);

Bool fifo_init(Fifo **fifo, long f_size);

void fifo_delete_fifo(Fifo **fifo);

void fifo_printf_fifo(Fifo *fifo);
#endif
