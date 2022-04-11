#ifndef _GRAPH_H
#define _GRAPH_H

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <list>
#include <stack>
#include <cmath>
#include <iterator>
#include <cstdio>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <float.h>
#include "common.h"

#define INFINI 1000000000

//#define EPS 1E-04
#define MON_INFINITY 10E+20
#define NIL_E (Edge *)0
#define NIL_N (Node *)0
#define EPSILON 1e-4

#define CAP_VALUE_CAP 1
#define CAP_VALUE_X 2
#define CAP_VALUE_1 3
#define CAP_VALUE_0 4
#define CAP_VALUE_0_FOR_FRAC_EDGE 5
#define CAP_VALUE_M_FOR_FRAC_EDGE 6

#define NO_0_EDGES 1
#define ALL_EDGES 2
#define NO_CAP_GREATER_THAN 3
#define N_SUR_2 210
#define N_TAB 420
#define M_TAB 88000
#define N_SUR_2 210

#define NODE_NOT_MARKED 0
#define NODE_MARKED 1
#define NODE_LIST_TERMINATED 2

#define MARKED 10

using namespace ::std;

class C_edge;
/*************************************************************************************************/
//					***	Classe des Noeuds	***
/*************************************************************************************************/
class C_node
{
public:
  int num;                 // node number
  char *name;              // node name
  double coord_x, coord_y; // coordonnees du sommet
  // int connexite;	//type de connexite

  // list<C_node*> nd_incid;
  list<C_edge *> ed_incid;

  C_node();
  ~C_node();
};
/*************************************************************************************************/
//					***	Classe des Arêtes	***
/*************************************************************************************************/
class C_edge
{
public:
  int num;       // edge number
  double length; // length of the edge
  double dist_euclid(C_node *v1, C_node *v2);
  C_node *end1; // first end of the edge
  C_node *end2; // second end of the edge

  C_edge();
  C_edge(C_node *v1, C_node *v2);
  ~C_edge();

  C_node *get_end1();
  C_node *get_end2();
};

/*************************************************************************************************/
//					***	Classe des Graphes	***
/*************************************************************************************************/
class C_graph
{
public:
  int nb_nodes; // number of nodes
  int nb_edges; // number of edges 
  std::vector<int> con; // connectivity types 
  vector<C_node *> Nodes; // pointers of nodes
  vector<C_edge *> Edges; // pointers of edges  

  C_graph();
  ~C_graph();

  void read_instance(const char *nom, bool sndlib_tsplib);
  void affiche();

  int **coupe2(int &i, int &j, int &k, int **tab);
  void coupe3(int &i, int &j, int &k, int **tab);
  void coupe4(int &i, int &j, int &k, int **tab); 
  void coupez(int itab, int max_jtab, int **tab, int **tabz);
};

/*************************************************************************************************/

typedef struct Node
{
  long id;
  struct Edge *first_edge;    /*Pointer to the head of the adjacency list*/
  struct Node *next_gr;       /*Suivant dans le graphe réduit*/
  long id_W;                  /*Id de l'ensemble Wi auquel appartient le sommet*/
  struct Node *next_W;        /*the following node of this node in the set of nodes to which they belongs to*/
  struct Node *n_sh;          /*pointeur sur le représentant du sommet dans le graphe réduit*/
  struct Node *next_sh;       /*suivant dans la liste des sommets qui sont contractés ensemble*/
  struct Edge *first_edge_gr; /*premier dans la liste d'incidence du sommet dans le graphe réduit*/

  /*Element annexe pour effectuer*/
  /*d'autres contractions*/
  struct Node *a_next_gr; /*Suivant dans le graphe réduit*/
  struct Node *a_n_sh;    /*pointeur sur le représentant du sommet dans le graphe réduit*/
  struct Node *a_next_sh; /*suivant dans la liste des sommets qui sont contractés ensemble*/

  /*Pour un deuxi�me niveau de sauvegarde*/
  struct Node *b_next_gr; /*Suivant dans le graphe réduit*/
  struct Node *b_n_sh;    /*pointeur sur le représentant du sommet dans le graphe réduit*/
  struct Node *b_next_sh; /*suivant dans la liste des sommets qui sont contractés ensemble*/
} Node;

typedef struct Edge
{
  struct Node *adjac;    /*The destination node to which this edge is adjacent*/
  struct Edge *next;     /*the following edge in the adjacency list to wich this edge belongs to*/
  struct Edge *back;     /*reverse edge*/
  double cap;            /*Capacity of the edge*/
  struct Node *adjac_gr; /*the destination node of this edge in the shrunk graph*/
  struct Edge *next_gr;  /*the following edge in an adjacency list in the shrunk graph*/

  float X;  /*Value of X(e)*/
  long num; /*number of this edge in the graph edge list*/

} Edge;

typedef struct Graph
{
  long n_Nodes;
  long m_Edges;
  Node *Nodes;
  Edge *Edges;
  std::vector<int> con; // connectivity types 
} Graph;

#define NO_MEM fprintf(stderr, "Unable to allocate memory\n");

/*typedef enum
  { FALSE = 0,
    TRUE  = 1
  } BOOL;*/

typedef struct HO_Node
{
  long id;
  struct HO_Edge *first_edge;
  /* in cyclic list of incident edges */
  struct HO_Edge *scan_ptr;
  /* next edge to be scanned when node
     will be visited again */
  /* subsequent entries for use by gmincut */
  long dist;
  double excess;
  struct HO_Node *bfs_link;   /* for one way BFS working queue */
  struct HO_Node *stack_link; /* for stack of active node */
  struct HO_Node *left_link;  /* for doubly linked lists */
  struct HO_Node *right_link; /* of dormant and W nodes  */
  BOOL unmarked;              /* while BFS in progress */
  BOOL in_S;                  /* in set of source nodes  */
  BOOL in_W;                  /* in set of W-valid nodes */
  BOOL shore1;                /* final mark for placement of
                                 node in one of the cut shores
       */
} HO_node;

typedef struct HO_Edge
{
  HO_node *adjac;       /* pointer to adjacent node */
  struct HO_Edge *next; /* in incidence list of node
      from which edge is emanating */
  struct HO_Edge *back; /* pointer to reverse edge */
  double cap;
  double rcap; /* residual capacity used in gmincut */
} HO_edge;

typedef struct HO_Graph
{
  long n_nodes;
  HO_node *nodes;
  long n_edges;
  HO_edge *edges;
} HO_graph;

#define HO_NILN (HO_node *)0
#define HO_NILE (HO_edge *)0
#define EPS 1.0E-10

BOOL initialize(HO_graph *gr);

long bfs0(HO_node *t);

void bfs1(HO_node *t);

BOOL gmincut_ho(HO_graph *gr, double *mincap, long *n_shore);

/*Initialize a graph from the given list of edges*/
BOOL InitHO(simpleEdge *TSrc, HO_graph **GDest, long n_Nodes, long m_Edges);

/*Calculate the global Mincut of the graph gr*/
/*and returns pointer to a structure which*/
/*represents the set W and its complementary W_b.*/
/*Its also contains the value of the global mincut.*/
Mincut *global_mincut_ho(simpleEdge *TSrc, long n, long m);

BOOL gmincut_ho(HO_graph *, double *, long *);

BOOL alloc_graph_ho(long n, long m, HO_graph **gr);

void dealloc_graph_ho(HO_graph *gr);
static Bool *getMark_Tab();

Bool AllocateGraph(long n, long m, Graph *g);                       /*Allocates memory for the internal graph representation*/
Bool DeleteGraph(Graph *g);                                         /*Deletes graph from memory*/
Bool InitGraph_from_file(char *fich, Graph *g, int *k);             /*Reads graph from file input*/
Bool InitGraph_from_list(simpleEdge *Te, long n, long m, Graph *g); /*Reads graph from input list*/
void PrintGraph(Graph *g);                                          /*Prints the graph*/
void PrintReducedGraph(Graph *g, FILE *fout);                       /*Print the shrunk graph*/
void PrintGraph_Table(Graph *g);

/*Affiche les ar�tes avec leurs sommets d'origine*/
void Print_b_Edge_Set_Gr(Graph *gr, b_Edge *f_e, FILE *out_file);

/*Renvoie la liste des ar�tes du graphe r�duit sans les ar�tes nulles*/
simpleEdge *GetReducedGraph_eList(Graph *g, long *m_edge);

Bool copy_graph(Graph *src, Graph *dest);

/*Give the list of edges of the graph*/
simpleEdge *Graph2SimpleEdge(Graph *g);

/*Return the list of the edges of the subgraph G[W]*/
/*n_tab will be the length of the returned list.*/
/*edge_flag says if you don't want the edges that have X(e)=0*/
/*cap_flag says what kind of value you want as capacity.*/
void get_edge_list(Graph *g, b_Node *w, simpleEdge *edge_tab, long *n_tab, int edge_flag, int cap_flag);

void Print_eList(simpleEdge *Tab, long n, long m);

/*Return a subset of nodes which constitutes a violated cut of the graph.*/
/*ghct represents the the Gomory-Hu cut tree given by the GHCutTree_Tree function.*/
/*b_cour represents an edge of the ghct that is violated*/
b_Node *get_cut_set(Tree *ghct, b_Edge *b_cour);

/*Recursive function for searching a violated cut in the graph.*/
/*It needs first to calculate the ghct associated to the graph.*/
/*b_cour is the ghct edge that doesn't satisfied the criterion*/
/*we have fixed and that we shouldn't explore.*/
/*f_list is first b_Node of the final list that will represent*/
/*the set W.*/
void create_cut_set(Tree_Node *t, b_Node **f_list, b_Edge *b_cour, Bool *found);

/*Create the complement of the set W that we have first built with create_cut_set()*/
/*If b_cour=uv is a bad edge of the ghct, with get_cut_set we create the set W above the*/
/*node above u.The complement of the set W will be the set under the node v.*/
void create_complement_cut_set(Tree_Node *t, b_Node **f_list, b_Edge *b_cour);

/*Return the set W_b of b_Node that should represent the complementary*/
/*of the set obtained with get_cut_set().*/
b_Node *get_complement_cut_set(Tree *ghct, b_Edge *b_cour);

/*Separate a cut using the ghct and create the two sets W and its complementary W_b*/
/*option says if you want to have only W, only W_b or both W and W_b.*/
/*f_w is the list of W and f_w_b is the list of w_b. If you choose W_ONLY then*/
/*the parameter f_w_b will not be used. You could give NULL as value of the parameter f_w_b*/
void separate_cut(Tree *ghct, b_Edge *b_cour, b_Node **f_w, b_Node **f_w_b, int option);

/*Contract the nodes which are in the same id*/
Bool contract_set_id(Graph *gr, long id_W);

/*Contract the set of nodes contained in the linked list f_w*/
Bool contract_set_w(Graph *gr, b_Node *f_w);

/*Decontract the node pointed by n_node*/
Bool decontract_node(Graph *gr, Node *n_node);

/*Merges two cyclic linked lists*/
/*The lists must have at least 1 node.*/
Bool cyclic_set_merge(Node *f1, Node *f2);

/*Checks if the cardinal of the set W is greater or equal to a.*/
/*You can use this function to check if |W| < b by checking if*/
/*the result of check_w_card(W,b) is false.*/
Bool check_w_card(b_Node *w, long a);

/*Checks if the cardinal of the set W is greater or equal to a.*/
/*You can use this function to check if |W| < b by checking if*/
/*the result of check_w_card(W,b) is false.*/
Bool check_delta_w_card(b_Edge *d_w, long a);

/*This function check if all the edges given by Tab*/
/*statisfied X(e)=a. It will return -1 if all the edges*/
/*statisfied the condition X(e)=1 and i if the i-th edge*/
/*doesn't satisfied this constraint.*/
long find_edge(simpleEdge *Tab, long m_tab, float a);

/*Says if the set of edges contains a fractionnary edge*/
/*It returns the index of the first fractionnary edge.*/
/*It will return -1 if there isn't any frac edge in the table.*/
long frac_edge_index(simpleEdge *Tab, long m_tab);

/*Gives the list of the edges of delta(v).*/
/*The method is similar to the method used in get_delta_w().*/
b_Edge *get_delta_v(Graph *g, long v, int edge_flag, int cap_flag);

/*Gives the list of the edges of delta(W).*/
/*The method is similar to the method used in get_edge_list().*/
b_Edge *get_delta_w(Graph *g, b_Node *w, int edge_flag, int cap_flag);

/*Calculate the connected components of a graph*/
b_Node **connected_components(simpleEdge *Te, long n_Te, long m_Te, long *n_comp, int edge_flag, double cap_value);

/*Add a new node in the connected components which number is n_comp*/
void add_to_connected_comp(b_Node **tab_comp, long *n_comp, long u);

/*Calculate the connected components of a graph without relabelling the nodes*/
/*You should use this function when your working on the original graph.*/
b_Node **graph_connected_components(Graph *gr, long *n_comp, int edge_flag);

/*Serch a fractionnary cycle on the graph given*/
/*by frac_edge_list.p is the cardinality of the cycle.*/
b_Node *find_frac_cycle(Graph *gr, long n, long m, long *p, Bool *marquage, Bool *cycle_type, short *status);

/*Retourne Vrai si l'entier z est pair et Faux sinon*/
Bool is_pair(long z);

long sup_part(long num, long denom);

Bool InitGraph_from_file(char *fich, Graph *g, int *k);             /*Reads graph from file input*/
Bool InitGraph_from_list(simpleEdge *Te, long n, long m, Graph *g); /*Reads graph from input list*/
Bool AllocateGraph(long n, long m, Graph *g);                       /*Allocates memory for the internal graph representation*/
Bool DeleteGraph(Graph *g);                                         /*Deletes graph from memory*/

/*This function returns a set of nodes which constitute a k-clique*/
void get_clique(Graph *gr, int k, int k_card_1, long k_card_2, b_Node *clique[], long *n_clique);

/*Permet de sauvegarder la contraction courante*/
/*dans les champs annexes*/
void save_contraction_info(Graph *gr, int niveau);

/*Contract the set of nodes contained in the linked list f_w*/
Bool contract_set_w(Graph *gr, b_Node *f_w);

/*Returns a list of subsets of V that should*/
/*be contracted by the operation Theta_4.*/
/*Notice that we work on the shrunk graph.*/
/*Tab is a list of the subset of V that*/
/*should be contracted by theta_4.*/
void operation_theta_4(Graph *gr, b_Node *w, int k, b_Node **Tab, long *n_subset);
/*REcherche les ar�tes multiples  � 1 et les contracte*/
void operation_theta_4_4(Graph *gr, b_Node *w, int k, b_Node **Tab, long *n_subset);
/*This function contracts the ((k/2)+1)-clique whose degree is equal to k+1*/
/*This function should work only for complete graphs such as TSPLIB graphs.*/
void operation_theta_4_bis(Graph *gr, b_Node *w, int k, b_Node **Tab, long *n_subset);

void operation_theta_3_bis(Graph *gr, b_Node *w, int k, b_Node **Tab, long *n_subset);

typedef struct GUS_Node
{
  long id;
  struct GUS_Edge *first_edge;
  /* in list of incident edges */
  struct GUS_Edge *scan_ptr;
  /* next edge to be scanned when node
     will be visited again */
  struct GUS_Node *parent; /* pointer of Gomory-Hu cut tree */
  double mincap;           /* capacity of minimum cut between
                  node and parent in GH cut tree */
  /* subsequent entries for use by maxflow */
  long dist;
  double excess;
  struct GUS_Node *bfs_link;   /* for one way BFS working queue */
  struct GUS_Node *stack_link; /* for stack of active nodes */
  BOOL alive;
  BOOL unmarked; /* while BFS in progress */
} GUS_node;

typedef struct GUS_Edge
{
  GUS_node *adjac;       /* pointer to adjacent node */
  struct GUS_Edge *next; /* in incidence list of node
      from which edge is emanating */
  struct GUS_Edge *back; /* pointer to reverse edge */
  double cap;
  double rcap; /* residual capacity */
} GUS_edge;

typedef struct GUS_Graph
{
  long n_nodes;
  GUS_node *nodes;
  long n_edges;  /* number of edges with non-zero capacity */
  long n_edges0; /* number of all edges */
  GUS_edge *edges;
} GUS_graph;

#define GUS_NILN (GUS_node *)0
#define GUS_NILE (GUS_edge *)0
/*#define  EPS  1.0E-10*/

/*Initialize a graph from the given list of edges*/
BOOL InitGusfield(simpleEdge *TSrc, GUS_graph **GDest, long n_Nodes, long m_Edges);

/*Calculate the Gomory-Hu cut tree*/
/*and return the list of the edges of the tree*/
simpleEdge *get_GHCT_eList(simpleEdge *TSrc, long n, long m, long *n_ghct);

/*Calculate the GH cut tree and return the tree under*/
/*the Tree data structure representation(cf simpleEdge.h)*/
Tree *GHCutTree(simpleEdge *TSrc, long n, long m, double min, double max);

/*Print the tree data structure given in parameter*/
void PrintGHCutTree(Tree *ghct);

/*Compute the st-maxflow with the Goldberg-Tarjan*/
/*algorithm implemented by Gusfield.The graph is*/
/*the simpleEdge list Te. m is the length of Te.*/
Mincut *st_maxflow(simpleEdge *Te, long n, long m, long s, long t);

float user_time();
BOOL ghc_tree(GUS_graph *);

BOOL alloc_graph(long n, long m, GUS_graph **gr);

void dealloc_graph(GUS_graph *gr);

char skip_comment(FILE *fd);

GUS_graph *get_graph(FILE *fd);

void print_tree(GUS_graph *gr);

BOOL init_maxflow(long n);

void global_relabel(GUS_graph *gr, GUS_node *tptr);

double maxflow(GUS_graph *gr, GUS_node *s_ptr, GUS_node *t_ptr);

BOOL ghc_tree(GUS_graph *gr);

void fini_maxflow();

/*Permet de marquer un sommet et tous les sommets*/
/*qui sont contract�s avec lui. Le marquage se fait*/
/*dans le tableau mark_t*/
void marquage_sommet(Graph *gr, long v, Bool *mark_t, Bool mark_val);


/*Calcule le plus court chemin entre*/
/*tous les couples de sommets dans le graphe G*/
/*On utilise l'algorithme de Floyd*/
long **Floyd(b_Edge *sh_eList, long n);

Bool separation_sp_partition_2(Graph *gr, b_Node *chaine[], long *chaine_sz, long nb_chaine, long ***sp_part_list, long *sp_p, long *nb_sp_part, long *rhs);

#endif
