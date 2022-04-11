
#ifndef Cuts_H
#define Cuts_H

#include <ilcplex/ilocplex.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <list>
#include <vector>
#include <algorithm>
#include <iterator>
#include <string.h>

#include ".././include/graph.h"
#include ".././include/Graph_Plot.h"
#include ".././include/Functions.h"

using namespace std;
 
struct cutinfo
{
   CPXLPptr model;
   int numcols;
   C_graph Gr;
   Graph *G_aux; // Utiliser pour les séparation provenant du kECSP.
   Graph_Plot *G_plot; 
   int K; // nbre d'arete
   string instanceName; 

   simpleEdge *bestSolution;
   double *bestValue;
};
typedef struct cutinfo CUTINFO, *CUTINFOptr;

struct lazycutinfo
{
   CPXLPptr model;
   int numcols;
   C_graph Gr;
   Graph *G_aux; // Utiliser pour les séparation provenant du kECSP.
   Graph_Plot *G_plot; 
   int K;
   string instanceName; 

   simpleEdge *bestSolution;
   double *bestValue;
};
typedef struct lazycutinfo LAZYCUTINFO, *CUTINFOptrLazy;

class C_cut
{
public:
   static vector<C_cut *> GeneratedEdgeCuts; 
   static vector<C_cut *> GeneratedSPPartition;
   static vector<C_cut *> GeneratedFPartition; 
   static vector<C_cut *> GeneratedPartition;

public:
   vector<int> vectInd;         // indices of var that appears in the cut
   vector<double> vectCoeff;    // coefficients of var that appears in the cut
   vector<string> VarNameArray; // names of var that appears in the cut
   int rhs;
   string type;
   C_cut();
   ~C_cut();

   static bool CutAlreadyExists(C_cut *coupe, vector<C_cut *> PoolsOfCuts);
};

vector<C_edge *> coupeMin(CPXCENVptr env, CPXLPptr model, C_graph *G, double *sol_y, int numNodeI, int numNodeII, double *value); // OK !

void updateCutinfo(CUTINFO *cutinfo, CPXLPptr model, int cur_numcols, C_graph G, Graph *G_aux, Graph_Plot *g_plot, simpleEdge *theBestSolution, double *theBestValue, int K, string instance);

void updateLazyCutinfo(LAZYCUTINFO *lazycutinfo, CPXLPptr model, int cur_numcols, C_graph G, Graph *G_aux, Graph_Plot *g_plot, simpleEdge *theBestSolution, double *theBestValue, int K, string instance);

int mycutcallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p);

int add_edge_cut_inequality(CPXCENVptr env, CPXLPptr model, void *cbdata, int wherefrom, C_graph G, double *sol_y);
int add_lazy_edge_cut_inequality(CPXCENVptr env, CPXLPptr model, void *cbdata, int wherefrom, C_graph G, double *sol_y);

int separation_globale(CPXCENVptr env, CPXLPptr model, void *cbdata, int wherefrom, void *cbhandle, C_graph G, double *sol_y, int type);
 
int AjouteSPPartition(CPXCENVptr env, CPXLPptr model, void *cbdata, void *cbhandle, int wherefrom, Graph *gr, b_Node *frac_cycle, long cycle_sz, FILE *sortie);
int AjouteFPartition(CPXCENVptr env, CPXLPptr model, void *cbdata, void *cbhandle, int wherefrom, Graph *gr, int k, b_Node *frac_cycle, long cycle_sz, FILE *sortie);
int AjoutePartition(CPXCENVptr env, CPXLPptr model, void *cbdata, void *cbhandle, int wherefrom, Graph *gr, int k, b_Node *frac_cycle, long cycle_sz, FILE *sortie);

void ProceedToContractions(Graph *gr, int k);
void RemoveAllContractions(Graph *gr);

/*Permet de r�initialiser la table des contractions*/
/*en rappelant la derni�re contraction sauvegard�e*/
void recall_contraction_info(Graph *gr, int niveau);
/*S�paration des F-partition*/
/*Retourne Vrai si F-Partition trouv�e et Faux sinon*/
Bool separation_f_partition(Graph *gr, int k, b_Node *frac_cycle, long cycle_sz, long **f_partition, long **V0_F, long *V0_F_sz, long *rhs, FILE *sortie_frac_gk);
/*S�paration des Partition*/
/*Retourne Vrai si Partition trouv�e et Faux sinon*/
Bool separation_partition(Graph *gr, int k, b_Node *frac_cycle, long cycle_sz, long **partition/*, long *rhs*/, FILE *sortie_frac_gk);
/*S�paration des SP-Partition*/
/*Retourne Vrai si SP-Partition trouv�e et Faux sinon*/
Bool separation_sp_partition(Graph *gr, int k, b_Edge **partition_edge, long *nb_sp_part, long *rhs, FILE *sortie_frac_gk);
Bool separation_sp_partition_halin(Graph *gr, int k, b_Node *cycle, b_Edge *sp_partition_edge[], long rhs[], long *nb_sp_partition, FILE *sortie_frac_gk);
/*S�paration des chemin impair*/
Bool separation_chemin_impair(Graph *gr, int k, b_Node *frac_cycle, long cycle_sz, b_Edge **L1, b_Edge **L2, b_Edge **L3, long *rhs, FILE *sortie_frac_gk);
void get_sp_partition_chaine(Graph *gr, b_Node **chaine, long *chaine_sz, long *nb_chaine);

#endif
