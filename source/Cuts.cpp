#include "../include/Cuts.h"
#include "../include/graph.h"
#include "../include/Functions.h"

#include <stdio.h>
#include <stdlib.h>
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
#include <stdlib.h>
#include <string>
#include <boost/concept_check.hpp>

vector<C_cut *> generatedCuts;
vector<C_cut *> generatedLazyCuts;

double cplexTimeBeforeSepMSI, cplexTimeAfterSepMSI;
int numcuts_edge = 0; 
double timeMSI = 0.0;

vector<C_cut *> C_cut::GeneratedEdgeCuts; 
vector<C_cut *> C_cut::GeneratedSPPartition;
vector<C_cut *> C_cut::GeneratedFPartition; 
vector<C_cut *> C_cut::GeneratedPartition;
vector<C_cut *> C_cut::GeneratedNewCut;

C_cut::C_cut()
{
	rhs = 0;
	type = "";
}

C_cut::~C_cut()
{
} 
bool C_cut::C_cut::CutAlreadyExists(C_cut *coupe, vector<C_cut *> PoolsOfCuts)
{
	int tailleCuts = PoolsOfCuts.size();

	if (tailleCuts == 0)
		return false;

	for (int i = 0; i < tailleCuts; i++)
	{
		if ((coupe->vectInd.size() == PoolsOfCuts[i]->vectInd.size()) && (coupe->rhs == PoolsOfCuts[i]->rhs))
		{

			int tailleIndices = PoolsOfCuts[i]->vectInd.size();
			int compteur = 0;

			while (compteur < tailleIndices && coupe->vectInd[compteur] == PoolsOfCuts[i]->vectInd[compteur])
			{

				compteur++;
			}

			if (compteur == tailleIndices)
				return true;
		}
	}

	return false;
}

//*********************************************************************************************************************************
//***************************************	coupeMin		***********************************************************
//*********************************************************************************************************************************

vector<C_edge *> coupeMin(CPXCENVptr env, CPXLPptr model, C_graph *G, double *sol_y, int numNodeI, int numNodeII, double *value)
{
	using namespace ::std;
	using namespace lemon;

	typedef lemon::ListGraph Graph; // création d'un type Graph (graphe indirect)

	typedef Graph::EdgeIt EdgeIt;
	typedef Graph::NodeIt NodeIt;
	typedef Graph::Node Node; // création d'un type Node lié à Graph (noeud)
	typedef Graph::Edge Edge; // création d'un type Edge lié à Graph (arête)

	// une Map sert à donner des valeurs aux composants (noeuds, arêtes, ...) d'un graphe !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	typedef Graph::EdgeMap<int> CapacityMap;	// Pour donner une capacité pour chaque arête
	typedef Graph::EdgeMap<int> RealNumEdgeMap; // Pour numéroter les arêtes
	typedef Graph::NodeMap<int> RealNumNodeMap; // Pour numéroter les noeuds

	using lemon::INVALID;

	int nb_nodes = G->nb_nodes;
	int nb_edges = G->nb_edges;

	Graph GP;

	vector<Node> NodesP;
	vector<Edge> EdgesP;

	CapacityMap pds1(GP); // capacité pour après
	RealNumEdgeMap pds2(GP);
	RealNumNodeMap pds3(GP);

	int poids;

	//*********************************	Création des noeuds	****************************************************
	for (int i1 = 0; i1 < nb_nodes; i1++)
	{
		NodesP.push_back(GP.addNode());
		pds3[NodesP[i1]] = G->Nodes[i1]->num; 
	}
 
	char name[100];
	for (int i2 = 0; i2 < nb_edges; i2++)
	{
		string s = getCollumName(env, model, i2);
		int i, j, k;
		if (s.at(0) == 'x')
		{
			strncpy(name, s.c_str(), 100);
			sscanf(name, "x%d(%d,%d)", &k, &i, &j);

			Edge a = GP.addEdge(NodesP[i], NodesP[j]);
			EdgesP.push_back(a);
			pds2[EdgesP[i2]] = G->Edges[k]->num;
			pds1[EdgesP[i2]] = (int)(sol_y[i2] * 1000);
		}
	}

	//*******************	algorithme GomoryHu pour le calcul de la coupe min	**************************************

	GomoryHu<Graph, CapacityMap> gom(GP, pds1);
	gom.run();

	//
	for (GomoryHu<Graph>::MinCutEdgeIt e(gom, NodesP[numNodeI], NodesP[numNodeII]); e != INVALID; ++e)
	{
		(*value) += pds1[e]; 
	} 
	vector<C_edge *> MinCut;
	C_edge *ptr; 
	for (GomoryHu<Graph>::MinCutEdgeIt e(gom, NodesP[numNodeI], NodesP[numNodeII]); e != INVALID; ++e)
	{
		ptr = new (C_edge);
		ptr = G->Edges[pds2[e]];
		MinCut.push_back(ptr);
		*value += ((double)pds1[e]) / 1000.0;
	}

	return MinCut;
}
//*********************************************************************************************************************************
//***************************************	updateCutinfo		***********************************************************
//*********************************************************************************************************************************

void updateCutinfo(CUTINFO *cutinfo, CPXLPptr model, int cur_numcols, C_graph G, Graph *G_aux,  simpleEdge *theBestSolution, double *theBestValue, int K, string instance)
{
	// cout<<"updateCutinfo !!!!"<<endl;
	cutinfo->model = model;
	cutinfo->numcols = cur_numcols;
	cutinfo->Gr = G;
	cutinfo->G_aux = G_aux; 
	cutinfo->K = K;
	cutinfo->instanceName = instance; 

	cutinfo->bestSolution = theBestSolution;
	cutinfo->bestValue = theBestValue;
}
//*********************************************************************************************************************************
//***************************************	updateCutinfo\end	***********************************************************
//*********************************************************************************************************************************

void updateLazyCutinfo(LAZYCUTINFO *lazycutinfo, CPXLPptr model, int cur_numcols, C_graph G, Graph *G_aux, simpleEdge *theBestSolution, double *theBestValue, int K, string instance)
{
	// cout<<"updateCutinfo !!!!"<<endl;
	lazycutinfo->model = model;
	lazycutinfo->numcols = cur_numcols;
	lazycutinfo->Gr = G;
	lazycutinfo->G_aux = G_aux; 
	lazycutinfo->K = K;
	lazycutinfo->instanceName = instance; 

	lazycutinfo->bestSolution = theBestSolution;
	lazycutinfo->bestValue = theBestValue;
}

//*********************************************************************************************************************************
//***************************************	add_edge_cut_inequality			***********************************************************
//*********************************************************************************************************************************

int add_edge_cut_inequality(CPXCENVptr env, CPXLPptr model, void *cbdata, int wherefrom, C_graph G, double *sol_y)
{
 	C_cut *ptr_cut;
	vector<C_edge *> coupe;
	double RHS;
	int status;
	double value = 0.0;
	int nbCtrViolee = 0;

	for (int i1 = 0; i1 < G.nb_nodes - 1; i1++)
	{
		for (int i2 = i1 + 1; i2 < G.nb_nodes; i2++)
		{
			coupe = coupeMin(env, model, &G, sol_y, G.Nodes[i1]->num, G.Nodes[i2]->num, &value);
			int nbActiveVar = coupe.size(); 
			int *cutind = new int[nbActiveVar];
			double *cutval = new double[nbActiveVar];
			ptr_cut = new (C_cut); 

			// Compute the set W et V\W from the cut
			std::vector<int> coupe_index;
			for(auto cou : coupe)
				coupe_index.push_back((*cou).num);
			std::vector< int> W(G.nb_nodes, 0);
			std::vector< int> COM_W(G.nb_nodes, 0);
			std::vector<bool> checked_link(G.nb_edges, 0);
 
			W[0] = 1;
			int TOTAL_NODES = 1;
			do
			{ 
				bool oneNodeAssigned = false;
				for(int link = 0; link < G.nb_edges; link++)
					if(checked_link[link] == 0)
					{ 
						bool notBelongToCut = (std::find(coupe_index.begin(), coupe_index.end(), link) == coupe_index.end()); 
						if( (W[G.Edges[link]->end1->num] == 1 ^ W[G.Edges[link]->end2->num] == 1) && notBelongToCut)
						{ 
							W[G.Edges[link]->end1->num] = 1;
							W[G.Edges[link]->end2->num] = 1; 
							TOTAL_NODES++;
							oneNodeAssigned = true;
							checked_link[link] = 1;
							break;
						}
						if( (COM_W[G.Edges[link]->end1->num] == 1 ^ COM_W[G.Edges[link]->end2->num] == 1) && notBelongToCut)
						{ 
							COM_W[G.Edges[link]->end1->num] = 1;
							COM_W[G.Edges[link]->end2->num] = 1; 
							TOTAL_NODES ++;
							oneNodeAssigned = true;
							checked_link[link] = 1;
							break;
						} 
						if( (W[G.Edges[link]->end1->num] == 1 && COM_W[G.Edges[link]->end2->num] == 0) && !notBelongToCut)
						{ 
							COM_W[G.Edges[link]->end2->num] = 1; 
							TOTAL_NODES ++;
							oneNodeAssigned = true;
							checked_link[link] = 1;
							break;
						} 
						if( (W[G.Edges[link]->end2->num] == 1 && COM_W[G.Edges[link]->end1->num] == 0) && !notBelongToCut)
						{ 
							COM_W[G.Edges[link]->end1->num] = 1; 
							TOTAL_NODES++;
							oneNodeAssigned = true;
							checked_link[link] = 1;
							break;
						} 
						if( (COM_W[G.Edges[link]->end1->num] == 1 && W[G.Edges[link]->end2->num] == 0) && !notBelongToCut)
						{ 
							W[G.Edges[link]->end2->num] = 1; 
							TOTAL_NODES++;
							oneNodeAssigned = true;
							checked_link[link] = 1;
							break;
						} 
						if( (COM_W[G.Edges[link]->end2->num] == 1 && W[G.Edges[link]->end1->num] == 0) && !notBelongToCut)
						{ 
							W[G.Edges[link]->end1->num] = 1; 
							TOTAL_NODES++;
							oneNodeAssigned = true;
							checked_link[link] = 1;
							break;
						}
					} 
				if(oneNodeAssigned == false)
				{ 
					for (int v = 0; v < G.nb_nodes; v++)
						if(W[v] == 0)
						{ 
							W[v] = 1; 
							TOTAL_NODES++;
						}
				}  
			} while (TOTAL_NODES != G.nb_nodes);
			
			/*std::cout << "W : ";
			for(auto v : W)
				std::cout << v << ", ";
			std::cout << std::endl;
			std::cout << "V\W : ";
			for(auto v : COM_W)
				std::cout << v << ", ";
			std::cout << std::endl;*/

			int r_W = 0, r_comp_W = 0;
			
			for (int v = 0; v < G.nb_nodes; v++)
			{
				if(W[v] == 1)
					r_W = std::max(r_W, G.con[v]);
				if(COM_W[v] == 1)
					r_comp_W = std::max(r_comp_W, G.con[v]);
			}

			RHS = std::min(r_W, r_comp_W);  
			ptr_cut->rhs = RHS;

			if (value < RHS)
			{
				for (int i3 = 0; i3 < nbActiveVar; i3++)
				{
					char varName[100];
					int colindex;
					sprintf(varName, "x%d(%d,%d)", coupe[i3]->num, coupe[i3]->end1->num, coupe[i3]->end2->num);
					CPXgetcolindex(env, model, varName, &colindex);
					cutind[i3] = colindex;
					cutval[i3] = 1.0;
					ptr_cut->vectInd.push_back(colindex);
					ptr_cut->vectCoeff.push_back(1.0);
					ptr_cut->VarNameArray.push_back(varName);
				}
				// voir si la coupe générée existe ou pas :
				if (!C_cut::CutAlreadyExists(ptr_cut, C_cut::GeneratedEdgeCuts))
				{
					nbCtrViolee++;

					C_cut::GeneratedEdgeCuts.push_back(ptr_cut);
					status = CPXcutcallbackadd(env, cbdata, wherefrom, nbActiveVar, RHS, 'G', cutind, cutval, 1);
					if (status)
						cout << "Failed to add cut\n";
				}
 
				string s = getCollumName(env, model, cutind[0]);
				for (int i4 = 1; i4 < nbActiveVar; i4++)
				{
					s = s + " + " + getCollumName(env, model, cutind[i4]);
				} 
				ptr_cut = NULL; 

				delete[] cutind;
				delete[] cutval;
				break;
			}

			delete ptr_cut;
			value = 0.0;
		}
	}

	numcuts_edge += nbCtrViolee; 
	int taille = generatedCuts.size(); 

	return nbCtrViolee;
}

//*********************************************************************************************************************************
//***************************************	add_lazy_edge_cut_inequality		***********************************************************
//*********************************************************************************************************************************

int add_lazy_edge_cut_inequality(CPXCENVptr env, CPXLPptr model, void *cbdata, int wherefrom, C_graph G, double *sol_y)
{
	C_cut *ptr_cut;
	vector<C_edge *> coupe;
	double RHS;
	int status;
	double value = 0.0;
	bool connectivity;
	int nbCtrViolee = 0;

	for (int i1 = 0; i1 < G.nb_nodes - 1; i1++)
	{
		for (int i2 = i1 + 1; i2 < G.nb_nodes; i2++)
		{
			coupe = coupeMin(env, model, &G, sol_y, G.Nodes[i1]->num, G.Nodes[i2]->num, &value); 
			int nbActiveVar = coupe.size(); 
			int *cutind = new int[nbActiveVar];
			double *cutval = new double[nbActiveVar];
			ptr_cut = new (C_cut); 
			RHS = std::min(G.con[i1],G.con[i2]); 
			ptr_cut->rhs = RHS;

			if (value < RHS)
			{
				for (int i3 = 0; i3 < nbActiveVar; i3++)
				{
					char varName[100];
					int colindex;
					sprintf(varName, "x%d(%d,%d)", coupe[i3]->num, coupe[i3]->end1->num, coupe[i3]->end2->num);
					CPXgetcolindex(env, model, varName, &colindex);
					cutind[i3] = colindex;
					cutval[i3] = 1.0;
					ptr_cut->vectInd.push_back(colindex);
					ptr_cut->vectCoeff.push_back(1.0);
					ptr_cut->VarNameArray.push_back(varName);
				} 
				nbCtrViolee++;
				generatedLazyCuts.push_back(ptr_cut); 
				status = CPXcutcallbackadd(env, cbdata, wherefrom, nbActiveVar, RHS, 'G', cutind, cutval, 1);
				if (status)
					cout << "Failed to add cut\n";
				ptr_cut = NULL;

				delete[] cutind;
				delete[] cutval;
				break;
			}

			delete ptr_cut;
			value = 0.0;
		}
	}

	numcuts_edge += nbCtrViolee; 

	return nbCtrViolee;
}


//*********************************************************************************************************************************
//***************************************	mycutcallback		***********************************************************
//*********************************************************************************************************************************

int mycutcallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int *useraction_p)
{
	static pthread_mutex_t cs_mutex = PTHREAD_MUTEX_INITIALIZER; // méthode de calcul (serie/parallele) a ne pas toucher
	int status = 0;												 // function return(status) 0:on n'ajoute pas de contraintes | 1:ont ajoute des contraintes

	int ndepth; // profondeur du noeuds

	CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_DEPTH, &ndepth); // information sur le noeud

	pthread_mutex_lock(&cs_mutex); // methode de calcul

	CUTINFOptr cutinfo = (CUTINFOptr)cbhandle;

	CPXLPptr model = cutinfo->model;
	int numcols = cutinfo->numcols;
	C_graph Gr = cutinfo->Gr;

	int K = cutinfo->K;

	string instance = cutinfo->instanceName;
	/// bool activeCutsOnly = cutinfo->activeCutsOnly;

	*useraction_p = CPX_CALLBACK_DEFAULT; //=0 -> on n'ajoute pas de contraintes | =1 -> on ajoute des contraintes

	double EpsForIntegrality = 0.000001;

	bool flagFrac = false; // Flag that checks if there is at least one variable with fractional cost
						   // bool feasible = false; // Flag that checks if the current solution satisfies capacity constraint (feasible)

	double *sol_y = new double[K];
	for (int i = 0; i < K; i++)
	{
		sol_y[i] = 0;
	}

	double *x = new double[numcols];
	CPXgetcallbacknodex(env, cbdata, wherefrom, x, 0, numcols - 1);

	double objValue;
	storeLPSolution(env, model, numcols - 1, x, sol_y, Gr, cutinfo->G_aux, &objValue); 

	flagFrac = isFrac(K, sol_y, EpsForIntegrality);
	const char *I = instance.c_str();
	if (flagFrac)
	{ 
		CPXgettime(env, &cplexTimeBeforeSepMSI);

		int nbCoupe = 0;

		nbCoupe = add_edge_cut_inequality(env, model, cbdata, wherefrom, Gr, sol_y);

		//********************************************************
		// Separation des F-partition et SP-Partition et Partition
		if (nbCoupe == 0)
		{ 
			// Affichage de la solution
			ProceedToContractions(cutinfo->G_aux, Gr.con[0]);

			simpleEdge *reducedGraphEdgeList;
			long m_ReducedGraphEdge;
			reducedGraphEdgeList = GetReducedGraph_eList(cutinfo->G_aux, &m_ReducedGraphEdge); 

			// Separation des contraintes
			separation_globale(env, model, cbdata, wherefrom, cbhandle, Gr, sol_y, 1); 

			RemoveAllContractions(cutinfo->G_aux);
		}
		//********************************************************

		*useraction_p = CPX_CALLBACK_SET;

		CPXgettime(env, &cplexTimeAfterSepMSI);
		timeMSI += cplexTimeAfterSepMSI - cplexTimeBeforeSepMSI;
	}
	else
	{
		CPXgettime(env, &cplexTimeBeforeSepMSI);

		int nbCoupe = 0; 
		nbCoupe = add_lazy_edge_cut_inequality(env, model, cbdata, wherefrom, Gr, sol_y);
		*useraction_p = CPX_CALLBACK_SET;

		CPXgettime(env, &cplexTimeAfterSepMSI);
		timeMSI += cplexTimeAfterSepMSI - cplexTimeBeforeSepMSI;

		if (nbCoupe == 0)
		{ // Alors on a une solution entière et réalisable
		  // On la compare donc avec la meilleure solution connue

			if (objValue < (*cutinfo->bestValue))
			{ // Dans ce cas, on met à jour la meilleure solution connue
 
				(*cutinfo->bestValue) = objValue;

				for (int i = 0; i < cutinfo->G_aux->m_Edges; i++) 
					cutinfo->bestSolution[i].cap = cutinfo->G_aux->Edges[i].X;  
			}
		}
	}
	//cout << "Nb Node Cuts generated: " << numcuts_node << endl;
	//cout << "Nb Edge Cuts generated: " << numcuts_edge << endl;
	delete[] sol_y;
	delete[] x;
	generatedCuts.clear();
	generatedLazyCuts.clear(); 
	pthread_mutex_unlock(&cs_mutex);
	return (status);
}

//*********************************************************
// Il n'est pas nécessaire d'écrire les séparations annexes
// pour les LazyConstraints car la LAZYCONSTRAINT Callback
// est appelée lorsque la solution est entière. Si elle
// vérifie déjà les coupes, alors elle vérifie aussi les
// F-partition, SP-partition et Partition.
//*********************************************************

int separation_globale(CPXCENVptr env, CPXLPptr model, void *cbdata, int wherefrom, void *cbhandle, C_graph G, double *sol_y, int type)
{
	// Séparation des F-partition
	b_Edge *f_part_list = NIL_BE;
	b_Edge *part_list = NIL_BE;
	b_Edge **sp_part_list;
	long rhs_f_part;
	long rhs_part;
	long *rhs_sp_part;
	Bool f_part_result, part_result, sp_part_result;

	long rhs_chemin = 0;
	b_Edge *Chemin_L1 = NIL_BE, *Chemin_L2 = NIL_BE, *Chemin_L3 = NIL_BE;
	Bool chemin_result;
	long n, m;
	Graph *gr, gr_frac;
	int k;
	Bool *marquage;
	long u, v, n_e_frac, cycle_sz;
	Bool cycle_type;
	simpleEdge *frac_edge_list;
	b_Node *frac_cycle;

	bool chemin_ok = true, fp_ok = true, p_ok = true, sp_ok = true;
	Bool chemin_total;

	int res_p, res_fp, res_sp, res_chem;
	int n_fp, n_p, n_sp, n_chem;

	n_fp = n_p = n_chem = n_sp = 0;

	// Récupération des informations provenant de l'extérieur de la fonction
	CUTINFOptr cutinfo = (CUTINFOptr)cbhandle;
	gr = cutinfo->G_aux;
	k = cutinfo->Gr.con[0];////////////////////////////////////////////////////////////////////////////////////////////////// TO CHECK

	n = gr->n_Nodes;
	m = gr->m_Edges;

	AjouteNewCut(env, model, cbdata, cbhandle, wherefrom, gr);

	// Frac_cycle n'est pas utilisé dans AjouteSP. Donc Frac_cycle == NIL_BN
	res_sp = 0;
	res_sp = AjouteSPPartition(env, model, cbdata, cbhandle, wherefrom, gr, frac_cycle, cycle_sz, stdout);

	if (res_sp == -1)
	{
		return -1;
	}
	else
	{
		n_sp += res_sp;
	}

	if (n_sp >= 1)
		return n_sp;
 
	//***************************************************
	//        Separation des autres contraintes

	// Recherche d'un cycle d'arêtes fractionnaires dans le graphe réduit 
	frac_edge_list = new simpleEdge[m];
	n_e_frac = 0;

	for (int i = 0; i < m; i++)
	{
		v = gr->Edges[i].adjac->id;
		u = gr->Edges[i].back->adjac->id;

		if ((0.0 < gr->Edges[i].X) && (gr->Edges[i].X < 1.0) && (gr->Nodes[u - 1].n_sh->id != gr->Nodes[v - 1].n_sh->id))
		{
			frac_edge_list[n_e_frac].node1 = gr->Nodes[u - 1].n_sh->id;
			frac_edge_list[n_e_frac].node2 = gr->Nodes[v - 1].n_sh->id;
			frac_edge_list[n_e_frac].cap = gr->Edges[i].X;

			n_e_frac++;
		}
	}

	// Il faut éviter la situation où le graphe
	// réduit ne contient plus d'arête fractionnaire:
	// Cela est possible dans la procédure annexe de séparation
	if (n_e_frac == 0)
	{
		// free(frac_edge_list);
		delete[] frac_edge_list;

		return 0;
	}

	// Initialisation du graphe dans lequel va s'effectuer la recherche
	if (!InitGraph_from_list(frac_edge_list, n, n_e_frac, &gr_frac))
	{
		puts("Separation Globale: problème lors de l'initialisation du graphe");

		delete[] frac_edge_list;

		return kECSP_False;
	}

	// Tableau de marquage pour les arêtes fractionnaires
	// marquage = (Bool *)malloc(n_e_frac*sizeof(Bool));
	marquage = new Bool[n_e_frac];

	for (int i = 0; i < n_e_frac; i++)
		marquage[i] = 0; 

	short status = 0;
	res_fp = res_p = res_chem = res_sp = 0;

	do
	{
		res_fp = res_p = res_chem = res_sp = 0;

		// Recherche de l'ensemble de sommets p
		status = 0;
		frac_cycle = find_frac_cycle(&gr_frac, n, n_e_frac, &cycle_sz, marquage, &cycle_type, &status);  
		//Print_b_Node_Set(frac_cycle, stdout); 

		// Dans le cas où on a tout exploré
		if(frac_cycle == NIL_BN)
		{
			//cout << "Cycle Fractionnaire NULL" << endl;

			free(frac_edge_list);
			free(marquage);

			return false;
		}

		if (status == 1)
		{
			// Tests sur la cardinalité de cycle_sz
			// Pour un cycle impair on lance la séparation des F-partition
			if (is_pair(cycle_sz) == kECSP_False)
			{
				res_fp = AjouteFPartition(env, model, cbdata, cbhandle, wherefrom, gr, k, frac_cycle, cycle_sz, stdout);
 
				if (res_fp == -1)
				{
					// free(frac_edge_list);
					// free(marquage);
					delete[] frac_edge_list;
					delete[] marquage;
					DeleteGraph(&gr_frac);
					Delete_b_Node_Set(&frac_cycle);

					return -1;
				}
				else
					n_fp += res_fp;
 

				if (res_fp == 0)
				{
					res_p = AjoutePartition(env, model, cbdata, cbhandle, wherefrom, gr, k, frac_cycle, cycle_sz, stdout);

					if (res_p == -1)
					{
						// free(frac_edge_list);
						// free(marquage);
						delete[] frac_edge_list;
						delete[] marquage;
						DeleteGraph(&gr_frac);
						Delete_b_Node_Set(&frac_cycle);

						return -1;
					}
					else
					{
						n_p += res_p;
					} 
				}
			}
			else
			{
				res_p = AjoutePartition(env, model, cbdata, cbhandle, wherefrom, gr, k, frac_cycle, cycle_sz, stdout);

				if (res_p == -1)
				{
					// free(frac_edge_list);
					// free(marquage);
					delete[] frac_edge_list;
					delete[] marquage;
					DeleteGraph(&gr_frac);
					Delete_b_Node_Set(&frac_cycle);

					return -1;
				}
				else
					n_p += res_p; 
			}
		}
		else
		{
			if (status == -1)
			{
				res_fp = AjouteFPartition(env, model, cbdata, cbhandle, wherefrom, gr, k, frac_cycle, cycle_sz, stdout);

				if (res_fp == -1)
				{
					// free(frac_edge_list);
					// free(marquage);
					delete[] frac_edge_list;
					delete[] marquage;
					DeleteGraph(&gr_frac);
					Delete_b_Node_Set(&frac_cycle);

					return -1;
				}
				else
					n_fp += res_fp;

				if (res_fp == 0)
				{
					// On essaie la 2e separation des F-Partition
					// On regroupe dans V0 tous les sommets qui ne
					// sont pas sérrés et les sommets à 1
					b_Edge *delta_v, *be_cour;
					b_Node *b_cour, *nouv;
					b_Node *frac_cycle_2;
					long cycle_sz_2;

					frac_cycle_2 = NIL_BN;
					b_cour = frac_cycle;
					cycle_sz_2 = 0;
					while (b_cour != NIL_BN)
					{
						delta_v = get_delta_v(gr, b_cour->id, NO_0_EDGES, CAP_VALUE_X);

						double X = 0;
						be_cour = delta_v;
						while (be_cour != NIL_BE)
						{
							X = X + be_cour->cap;

							be_cour = be_cour->next;
						}

						Delete_b_Edge_Set(&delta_v); 

						if (X <= (double)(k + EPSILON))
						{
							nouv = new b_Node;
							nouv->id = b_cour->id;

							insert_b_Node(&frac_cycle_2, &nouv);
							cycle_sz_2++;
						}

						b_cour = b_cour->next;
					} 
					//Print_b_Node_Set(frac_cycle_2, stdout); 

					// create_branch_file(Te,m);

					// res_fp = 0;
					res_fp = AjouteFPartition(env, model, cbdata, cbhandle, wherefrom, gr, k, frac_cycle_2, cycle_sz_2, stdout);
					Delete_b_Node_Set(&frac_cycle_2);

					if (res_fp == -1)
					{
						// free(frac_edge_list);
						// free(marquage);
						delete[] frac_edge_list;
						delete[] marquage;
						DeleteGraph(&gr_frac);
						Delete_b_Node_Set(&frac_cycle);

						return -1;
					}
					else
						n_fp += res_fp;

					if (res_fp >= 1)
					{
						//cout << "F-Part2 violée" << endl;
						////create_branch_file(Te,m);
					}

					// Si là aussi on échoue alors, on sépare les partitions
					if (res_fp == 0)
					{
						res_p = AjoutePartition(env, model, cbdata, cbhandle, wherefrom, gr, k, frac_cycle, cycle_sz, stdout);

						if (res_p == -1)
						{
							// free(frac_edge_list);
							// free(marquage);
							delete[] frac_edge_list;
							delete[] marquage;
							DeleteGraph(&gr_frac);
							Delete_b_Node_Set(&frac_cycle);

							return -1;
						}
						else
							n_p += res_p;
					}
				}
			}
		}

		Delete_b_Node_Set(&frac_cycle);

	} while (status != -1);

	// On essaie de séparer les F-partition en
	// utilisant GHCT et en mettant 1-x comme capacité
	if ((n_fp + n_p + n_chem) == 0)
	{
		Tree *ghct_f;
		long sh_edge = 0;
		simpleEdge *sh_eList = GetReducedGraph_eList(gr, &sh_edge);
		b_Edge *be_cour;
		b_Node *W1, *W2, *cour;

		// Il ne faut pas oublier de rajouter EPSILON
		for (int i = 0; i < sh_edge; i++)
		{
			sh_eList[i].cap = 1 - sh_eList[i].cap + EPSILON;
		} 

		ghct_f = GHCutTree(sh_eList, n, sh_edge, 0.0, (double)k - EPSILON); 
		// free(sh_eList);
		delete[] sh_eList; 

		be_cour = ghct_f->b_List;
		while (be_cour != NIL_BE)
		{
			W1 = get_cut_set(ghct_f, be_cour);
			W2 = get_complement_cut_set(ghct_f, be_cour);
  
			if (check_w_card(W1, 2) && check_w_card(W2, 2))
			{
				frac_cycle = W1;
				cour = frac_cycle;
				cycle_sz = 0;
				while (cour != NIL_BN)
				{
					cycle_sz++;

					cour = cour->next;
				}

				res_fp = AjouteFPartition(env, model, cbdata, cbhandle, wherefrom, gr, k, frac_cycle, cycle_sz, stdout);

				if (res_fp == -1)
				{
					// free(frac_edge_list);
					// free(marquage);
					delete[] frac_edge_list;
					delete[] marquage;
					DeleteGraph(&gr_frac);
					Delete_b_Node_Set(&frac_cycle);
					Delete_Tree(&ghct_f);

					return -1;
				}
				else
					n_fp += res_fp;

				if (res_fp == 0)
				{
					// On essaie la 2e separation des F-Partition
					// On regroupe dans V0 tous les sommets qui ne
					// sont pas sérrés et les sommets à 1
					b_Edge *delta_v, *be_cour;
					b_Node *b_cour, *nouv;
					b_Node *frac_cycle_2;
					long cycle_sz_2;

					frac_cycle_2 = NIL_BN;
					b_cour = frac_cycle;
					cycle_sz_2 = 0;
					while (b_cour != NIL_BN)
					{
						delta_v = get_delta_v(gr, b_cour->id, NO_0_EDGES, CAP_VALUE_X);

						double X = 0;
						be_cour = delta_v;
						while (be_cour != NIL_BE)
						{
							X = X + be_cour->cap;

							be_cour = be_cour->next;
						}

						Delete_b_Edge_Set(&delta_v); 
						if (X <= (double)(k + EPSILON))
						{
							nouv = new b_Node;
							nouv->id = b_cour->id;

							insert_b_Node(&frac_cycle_2, &nouv);
							cycle_sz_2++;
						}

						b_cour = b_cour->next;
					}
 
					//Print_b_Node_Set(frac_cycle_2, stdout); 

					// create_branch_file(Te,m);

					// res_fp = 0;
					res_fp = AjouteFPartition(env, model, cbdata, cbhandle, wherefrom, gr, k, frac_cycle_2, cycle_sz_2, stdout);
					Delete_b_Node_Set(&frac_cycle_2);

					if (res_fp == -1)
					{
						// free(frac_edge_list);
						// free(marquage);
						delete[] frac_edge_list;
						delete[] marquage;
						DeleteGraph(&gr_frac);
						Delete_b_Node_Set(&frac_cycle);
						Delete_Tree(&ghct_f);

						return -1;
					}
					else
						n_fp += res_fp;

					if (res_fp >= 1)
					{
						//cout << "F-Part2 violée" << endl;
						////create_branch_file(Te,m);
					}

					// Si là aussi on échoue alors, on sépare les partitions
					if (res_fp == 0)
					{
						res_p = AjoutePartition(env, model, cbdata, cbhandle, wherefrom, gr, k, frac_cycle, cycle_sz, stdout);

						if (res_p == -1)
						{
							// free(frac_edge_list);
							// free(marquage);
							delete[] frac_edge_list;
							delete[] marquage;
							DeleteGraph(&gr_frac);
							Delete_b_Node_Set(&frac_cycle);
							Delete_Tree(&ghct_f);

							return -1;
						}
						else
							n_p += res_p;
					}
				}

				// On sépare sur W2
				frac_cycle = W2;
				cour = frac_cycle;
				cycle_sz = 0;
				while (cour != NIL_BN)
				{
					cycle_sz++;

					cour = cour->next;
				}

				res_fp = AjouteFPartition(env, model, cbdata, cbhandle, wherefrom, gr, k, frac_cycle, cycle_sz, stdout);

				if (res_fp == -1)
				{
					// free(frac_edge_list);
					// free(marquage);
					delete[] frac_edge_list;
					delete[] marquage;
					DeleteGraph(&gr_frac);
					Delete_b_Node_Set(&frac_cycle);
					Delete_Tree(&ghct_f);

					return -1;
				}
				else
					n_fp += res_fp;

				if (res_fp == 0)
				{
					// On essaie la 2e separation des F-Partition
					// On regroupe dans V0 tous les sommets qui ne
					// sont pas sérrés et les sommets à 1
					b_Edge *delta_v, *be_cour;
					b_Node *b_cour, *nouv;
					b_Node *frac_cycle_2;
					long cycle_sz_2;

					frac_cycle_2 = NIL_BN;
					b_cour = frac_cycle;
					cycle_sz_2 = 0;
					while (b_cour != NIL_BN)
					{
						delta_v = get_delta_v(gr, b_cour->id, NO_0_EDGES, CAP_VALUE_X);

						double X = 0;
						be_cour = delta_v;
						while (be_cour != NIL_BE)
						{
							X = X + be_cour->cap;

							be_cour = be_cour->next;
						}

						Delete_b_Edge_Set(&delta_v); 

						if (X <= (double)(k + EPSILON))
						{
							nouv = new b_Node;
							nouv->id = b_cour->id;

							insert_b_Node(&frac_cycle_2, &nouv);
							cycle_sz_2++;
						}

						b_cour = b_cour->next;
					}
 
					//Print_b_Node_Set(frac_cycle_2, stdout);
					
					// create_branch_file(Te,m);

					// res_fp = 0;
					res_fp = AjouteFPartition(env, model, cbdata, cbhandle, wherefrom, gr, k, frac_cycle_2, cycle_sz_2, stdout);
					Delete_b_Node_Set(&frac_cycle_2);

					if (res_fp == -1)
					{
						// free(frac_edge_list);
						// free(marquage);
						delete[] frac_edge_list;
						delete[] marquage;
						DeleteGraph(&gr_frac);
						Delete_b_Node_Set(&frac_cycle);
						Delete_Tree(&ghct_f);

						return -1;
					}
					else
						n_fp += res_fp;

					if (res_fp >= 1)
					{
						//cout << "F-Part2 violée" << endl;
						////create_branch_file(Te,m);
					}

					// Si là aussi on échoue alors, on sépare les partitions
					if (res_fp == 0)
					{
						res_p = AjoutePartition(env, model, cbdata, cbhandle, wherefrom, gr, k, frac_cycle, cycle_sz, stdout);

						if (res_p == -1)
						{
							// free(frac_edge_list);
							// free(marquage);
							delete[] frac_edge_list;
							delete[] marquage;
							DeleteGraph(&gr_frac);
							Delete_b_Node_Set(&frac_cycle);
							Delete_Tree(&ghct_f);

							return -1;
						}
						else
							n_p += res_p;
					}
				}
			}

			Delete_b_Node_Set(&W1);
			Delete_b_Node_Set(&W2);

			be_cour = be_cour->next;
		}

		Delete_Tree(&ghct_f);
 
		if ((n_fp + n_p + n_chem) >= 1)
		{
			//cout << "Separation GHCT réussie" << endl;

			// create_branch_file(Te,m);
		}
	}

	// Libération de la mémoire
	// free(frac_edge_list);
	// free(marquage);
	delete[] frac_edge_list;
	delete[] marquage;
	DeleteGraph(&gr_frac);
	// Delete_b_Node_Set(&frac_cycle);

	if((n_fp + n_p) == 0)
	{
		// Frac_cycle n'est pas utilisé dans AjouteSP. Donc Frac_cycle == NIL_BN
		//	res_sp = AjouteSPPartition(gr,k,frac_cycle,cycle_sz,sortie);

		if (res_sp == -1)
		{
			return -1;
		}
		else
		{
			n_sp += res_sp;
		}
	} 
	return (n_fp + n_p + n_sp);

	//***************************************************
}

int AjouteNewCut(CPXCENVptr env, CPXLPptr model, void *cbdata, void *cbhandle, int wherefrom, Graph *gr)
{
	/*
	lemon::ListGraph::ArcMap<double> weights((*gr).lemonGraph, 1.0);
	lemon::ListGraph::ArcMap<double> capacities((*gr).lemonGraph, 1.0);
	lemon::ListGraph::ArcMap<double> flows((*gr).lemonGraph);

	for (lemon::ListGraph::EdgeIt i((*gr).lemonGraph); i!=lemon::INVALID; ++i)
	{
		capacities[(*gr).lemonGraph.arcFromId((*gr).lemonGraph.id(i))] = 0.0;
		lemon::NetworkSimplex<lemon::ListGraph, double, double> ns((*gr).lemonGraph);
    	ns.costMap(weights).upperMap(capacities).stSupply((*gr).lemonGraph.u(i), (*gr).lemonGraph.v(i), 2);

		lemon::NetworkSimplex<lemon::ListGraph, double, double>::ProblemType status = ns.run(); 
		if(status == lemon::NetworkSimplex<lemon::ListGraph, double, double>::OPTIMAL)
		{
			ns.flowMap(flows);
			std::cout << "Solution found for link " <<  (*gr).lemonGraph.id(i) << "   : " << ns.totalCost() << std::endl;
			for (lemon::ListGraph::EdgeIt j((*gr).lemonGraph); j!=lemon::INVALID; ++j)
				//if(ns.flow((*gr).lemonGraph.arcFromId((*gr).lemonGraph.id(j))) > 0)
					std::cout << (*gr).lemonGraph.id(j) << " : " << ns.flow((*gr).lemonGraph.arcFromId((*gr).lemonGraph.id(j))) << std::endl; 
		} 
	}*/
	return 0;
}



int AjouteSPPartition(CPXCENVptr env, CPXLPptr model, void *cbdata, void *cbhandle, int wherefrom, Graph *gr, b_Node *frac_cycle, 
long cycle_sz, FILE *sortie)
{

	// return 0;

	bool sp_ok = false;
	Bool sp_part_result = kECSP_False;
	long nb_sp_part = 0;
	long n, m;
	n = gr->n_Nodes;
	m = gr->m_Edges;

	long nb_chaine = 0;
	b_Node **chaine;
	long *chaine_sz;
	long **sp_part_list;
	long *rhs;
	long *sp_p;

	int nbCtrViolee = 0;
	int status;

	chaine = new b_Node *[(n / 2)];
	chaine_sz = new long[(n / 2)];
	rhs = new long[(n / 2)];
	sp_p = new long[(n / 2)];
	sp_part_list = new long *[(n / 2)];

	for (int i = 0; i < (n / 2); i++)
	{
		chaine[i] = NIL_BN;
		chaine_sz[i] = 0;

		rhs[i] = 0;
		sp_part_list[i] = NULL;
		sp_p[i] = 0;
	}

	CUTINFOptr cutinfo = (CUTINFOptr)cbhandle;

	save_contraction_info(gr, 2);
	
	get_sp_partition_chaine(gr, chaine, chaine_sz, &nb_chaine);
	/*
	std::cout << "***Nmber Components " << nb_chaine << std::endl;
	for (int j = 0; j < nb_chaine; j++)
	{
		k = 0;
		b_Node* b_cour = chaine[j];
		while (b_cour != NIL_BN)
		{
			std::cout << "  " << b_cour->id << ", ";

			k++;
			b_cour = b_cour->next;
		}
		std::cout << std::endl;
	}*/

	simpleEdge *reducedGraphEdgeList;
	long m_ReducedGraphEdge;
	reducedGraphEdgeList = GetReducedGraph_eList(gr, &m_ReducedGraphEdge);

	sp_part_result = separation_sp_partition_2(gr, chaine, chaine_sz, nb_chaine, &sp_part_list, sp_p, &nb_sp_part, rhs);

/*
if(nb_sp_part>0)
{
std::cout << "Nmber Components " << nb_sp_part << std::endl;
	for (int j = 0; j < nb_sp_part; j++)
	{
		for (int l = 0; l < sp_p[j]; l++)
		{
			std::cout << "  " << sp_part_list[j][l] << ", ";
		}
		std::cout << std::endl;
	}
}*/
	// On rappelle le 2e niveau de contraction
	recall_contraction_info(gr, 2);

	sp_ok = true;

	if (sp_part_result == kECSP_True)
	{
		// Ajout de la contrainte dans le pool
		C_cut *new_sp;

		for (long i = 0; i < nb_sp_part; i++) // There is several parition 
		{
			/*
			std::cout << std::endl << " partition ";
			for (int j = 0; j < n; j++)
			{
				std::cout << "( " << sp_part_list[i][j] << ", " << (*gr).con[j] << ") , " ;
			}
			std::cout << std::endl;*/

			int max_set = 0;
			for (int j = 0; j < n; j++)
				if(sp_part_list[i][j]> max_set)
					max_set = sp_part_list[i][j];

			std::vector<int> r_V(max_set, 0); //Connectivty type of every set in the partition (max of r_v of all element in V)
			for (int j = 0; j < n; j++)
				r_V[sp_part_list[i][j]] = std::max(r_V[sp_part_list[i][j]], (*gr).con[j]);
			
			std::vector<int> con_V(max_set, 0); // con(W ) = min{r(W ), r(V \W )}
			for (int l = 0; l < max_set; l++)
			{
				int max_val = 0;
				for (int j = 0; j < max_set; j++)
					if(l!=j) 
						max_val = std::max(max_val, r_V[l]);
				con_V[l] = std::min(r_V[l], max_val); 
			}
			int p_1= 0, p_2 = 0, p_3 = 0, r_W = 0; // Number of sets in the partition with connectivity type 1, 2 and 3 respectivly.
			for (int j = 0; j < max_set; j++)
			{
				if(con_V[j] == 1)
					p_1++;
				if(con_V[j] == 2)
					p_2++;
				if(con_V[j] == 3)
					p_3++;
			}
			r_W = *std::max_element(con_V.begin(), con_V.end());
			/*std::cout <<std::endl<< "********************" <<std::endl;
			std::cout << p_1 << "  " << p_2 << "  " << p_3 << "  " << r_W << std::endl;
			std::cout <<std::endl<< "********************" <<std::endl;*/

			new_sp = new C_cut;
			int nbActiveVar = m;

			int *cutind_aux = new int[nbActiveVar];
			double *cutval_aux = new double[nbActiveVar];

			int nbCutNZ = 0;

			for (int i3 = 0; i3 < m; i3++)
			{ 
				// Calcul du coefficient de chaque variable
				double coeff;

				int t = gr->Edges[i3].back->adjac->id; // source of the link
				int h = gr->Edges[i3].adjac->id; // Destnation of the link

				if (sp_part_list[i][t - 1] != sp_part_list[i][h - 1])
				{
					coeff = 1.0;
				/*	if ((sp_part_list[i][t - 1] == sp_p[i]) || (sp_part_list[i][h - 1] == sp_p[i]))
					{
						coeff = 1.0;
					}
					else
					{
						if (fabs(sp_part_list[i][t - 1] - sp_part_list[i][h - 1]) <= 1)
						{
							coeff = 1.0;
						}
						else
						{
							coeff = 1.0;//2.0;
						}
					}*/
				}
				else
				{
					coeff = 0.0;
				}

				char varName[100];
				int colindex;
				sprintf(varName, "x%d(%d,%d)", i3, t - 1, h - 1);
				CPXgetcolindex(env, model, varName, &colindex);

				cutind_aux[i3] = colindex;
				cutval_aux[i3] = coeff;

				if (coeff != 0)
					nbCutNZ++;
			}

			double *cutval = new double[nbCutNZ];
			int *cutind = new int[nbCutNZ];

			if(r_W == 1)
				new_sp->rhs = p_1 - 1;
			if(r_W == 2)
				new_sp->rhs = p_1 + p_2;
			if(r_W == 3)
				new_sp->rhs = p_1 + p_2 + 2 * p_3 - 1;

			rhs[i] = new_sp->rhs;
			int indiceNZ = 0;
			double valueCut = 0.0;
			for (int i3 = 0; i3 < m; i3++) 
				if (cutval_aux[i3] > 0.5)
				{
					cutval[indiceNZ] = cutval_aux[i3];
					cutind[indiceNZ] = cutind_aux[i3];

					new_sp->vectInd.push_back(cutind_aux[i3]);
					//std::cout << "x_" << i3 << " + ";
					new_sp->vectCoeff.push_back(cutval_aux[i3]); 

					valueCut += (*gr).Edges[i3].X;

					indiceNZ++;
				} 
			//std::cout << " >= " << rhs[i] << std::endl;
			// voir si la coupe générée existe ou pas :
			if (!C_cut::CutAlreadyExists(new_sp, C_cut::GeneratedSPPartition) && valueCut < (rhs[i] - 1e-04))
			{
				nbCtrViolee++;

				C_cut::GeneratedSPPartition.push_back(new_sp);
				status = CPXcutcallbackadd(env, cbdata, wherefrom, indiceNZ, rhs[i], 'G', cutind, cutval, 1);
				if (status)
					cout << "Failed to add cut\n";
			}
			else
			{ 
				delete new_sp;
			}

			delete[] cutind;
			delete[] cutval;
			delete[] cutval_aux;
			delete[] cutind_aux;
		}
	}
	else
	{
		nb_sp_part = 0;
	}

	for (int i = 0; i < nb_chaine; i++)
	{
		Delete_b_Node_Set(&chaine[i]);
	} 
	delete[] chaine;
	delete[] chaine_sz;
	delete[] sp_part_list;
	delete[] sp_p;
	delete[] rhs; 

	return nbCtrViolee;
}

int AjouteFPartition(CPXCENVptr env, CPXLPptr model, void *cbdata, void *cbhandle, int wherefrom, Graph *gr, int k, b_Node *frac_cycle, long cycle_sz, FILE *sortie)
{

	long *f_part_list = NULL;
	long *V0_F, V0_F_sz;

	long rhs_f_part;
	Bool f_part_result;
	int nb_f_part = 0;
	bool fp_ok;

	fp_ok = true;

	int nbCtrViolee = 0;
	int status;

	int n = gr->n_Nodes;
	int m = gr->m_Edges;

	// Séparation des F-Partition
	f_part_result = kECSP_False;
	f_part_result = separation_f_partition(gr, k, frac_cycle, cycle_sz, &f_part_list, &V0_F, &V0_F_sz, &rhs_f_part, sortie);

	if (f_part_result == kECSP_True)
	{
		// Ajout de la contrainte dans le pool
		C_cut *new_fp;

		// cout << "F-partition RHS " << rhs_f_part << endl;

		new_fp = new C_cut;
		int nbActiveVar = m;

		int *cutind_aux = new int[nbActiveVar];
		double *cutval_aux = new double[nbActiveVar];

		int nbCutNZ = 0;
		for (int i3 = 0; i3 < m; i3++)
		{

			// Calcul du coefficient de chaque variable
			double coeff;
			bool dansF = false;

			int t = gr->Edges[i3].back->adjac->id;
			int h = gr->Edges[i3].adjac->id;

			if (f_part_list[h - 1] != f_part_list[t - 1])
			{
				if ((f_part_list[h - 1] != 0) && (f_part_list[t - 1] != 0))
				{
					coeff = 1.0;
				}
				else
				{
					dansF = false;
					int i = 0;

					while ((i < V0_F_sz) && (!dansF))
					{
						if (V0_F[i] == i3 + 1)
							dansF = true;

						i++;
					}

					if (dansF)
					{
						coeff = 0.0;
					}
					else
					{
						coeff = 1.0;
					}
				}
			}
			else
			{
				coeff = 0.0;
			}

			char varName[100];
			int colindex;
			sprintf(varName, "x%d(%d,%d)", i3, t - 1, h - 1);
			CPXgetcolindex(env, model, varName, &colindex);
			cutind_aux[i3] = colindex;
			cutval_aux[i3] = coeff;

			if (coeff != 0)
				nbCutNZ++;
		}

		double *cutval = new double[nbCutNZ];
		int *cutind = new int[nbCutNZ];

		new_fp->rhs = rhs_f_part;

		int indiceNZ = 0;
		for (int i3 = 0; i3 < m; i3++)
		{

			if (cutval_aux[i3] != 0)
			{
				cutval[indiceNZ] = cutval_aux[i3];
				cutind[indiceNZ] = cutind_aux[i3];
				//std::cout << cutval[indiceNZ] << " x_" << cutind[indiceNZ] << " + ";
				new_fp->vectInd.push_back(cutind_aux[i3]);
				new_fp->vectCoeff.push_back(cutval_aux[i3]);
				// new_fp->VarNameArray.push_back(varName);

				indiceNZ++;
			}
		}
		//std::cout << " >= " << new_fp->rhs << std::endl;
		// voir si la coupe générée existe ou pas :
		if (!C_cut::CutAlreadyExists(new_fp, C_cut::GeneratedFPartition))
		{
			nbCtrViolee++;

			C_cut::GeneratedFPartition.push_back(new_fp);
			status = CPXcutcallbackadd(env, cbdata, wherefrom, indiceNZ, rhs_f_part, 'G', cutind, cutval, 1);
			if (status)
				cout << "Failed to add cut\n";
		}
		else
		{

			//cout << "F-Partition - Contrainte déjà obtenue" << endl;

			delete new_fp;
		}

		delete[] cutind;
		delete[] cutval;
		delete[] cutind_aux;
		delete[] cutval_aux;
	}

	// cout << "NbCtrViolee == " << nbCtrViolee << endl;

	/*if (nbCtrViolee != 0)
		cout << "Nb F-Partition == " << nbCtrViolee << endl;*/

	return nbCtrViolee;
}


int AjoutePartition(CPXCENVptr env, CPXLPptr model, void *cbdata, void *cbhandle, int wherefrom, Graph *gr, int k, b_Node *frac_cycle, long cycle_sz, FILE *sortie)
{
	// return 0;

	long *part_list = NULL;
	long rhs_part;
	Bool part_result;
	int nb_part = 0;
	bool p_ok = true;

	int n = gr->n_Nodes;
	int m = gr->m_Edges;

	int nbCtrViolee = 0;
	int status;

	part_result = kECSP_False;
	part_result = separation_partition(gr, k, frac_cycle, cycle_sz, &part_list, sortie);

	if (part_result == kECSP_True)
	{
		// Ajout de la contrainte dans le pool
		C_cut *new_p;

		// cout << "*Partition Violee - RHS == " << rhs_part << endl;

		/****************************************************************************************************/
		int max_set = 0;
		for (int j = 0; j < n; j++)
			if(part_list[j]> max_set)
				max_set = part_list[j];

		std::vector<int> r_V(max_set, 0); //Connectivty type of every set in the partition (max of r_v of all element in V)
		for (int j = 0; j < n; j++)
			r_V[part_list[j]] = std::max(r_V[part_list[j]], (*gr).con[j]);
		
		std::vector<int> con_V(max_set, 0); // con(W ) = min{r(W ), r(V \W )}
		for (int l = 0; l < max_set; l++)
		{
			int max_val = 0;
			for (int j = 0; j < max_set; j++)
				if(l!=j) 
					max_val = std::max(max_val, r_V[l]);
			con_V[l] = std::min(r_V[l], max_val); 
		}
		int p_1= 0, p_2 = 0, p_3 = 0, r_W = 0; // Number of sets in the partition with connectivity type 1, 2 and 3 respectivly.
		for (int j = 0; j < max_set; j++)
		{
			if(con_V[j] == 1)
				p_1++;
			if(con_V[j] == 2)
				p_2++;
			if(con_V[j] == 3)
				p_3++;
		}
		r_W = *std::max_element(con_V.begin(), con_V.end());
		/*std::cout <<std::endl<< "********************" <<std::endl;
		std::cout << p_1 << "  " << p_2 << "  " << p_3 << "  " << r_W << std::endl;
		std::cout <<std::endl<< "********************" <<std::endl;*/
		/************************************************************************************************/
		new_p = new C_cut;
		int nbActiveVar = m;
		int *cutind_aux = new int[nbActiveVar];
		double *cutval_aux = new double[nbActiveVar];

		int nbCutNZ = 0;

		for (int i3 = 0; i3 < m; i3++)
		{

			// Calcul du coefficient de chaque variable
			double coeff;

			int t = gr->Edges[i3].back->adjac->id;
			int h = gr->Edges[i3].adjac->id;

			if (part_list[t - 1] != part_list[h - 1])
			{
				coeff = 1.0;
			}
			else
			{
				coeff = 0.0;
			}

			char varName[100];
			int colindex;
			sprintf(varName, "x%d(%d,%d)", i3, t - 1, h - 1);
			CPXgetcolindex(env, model, varName, &colindex);
			cutind_aux[i3] = colindex;
			cutval_aux[i3] = coeff;

			if (coeff != 0)
				nbCutNZ++;
		}

		double *cutval = new double[nbCutNZ];
		int *cutind = new int[nbCutNZ];

		if( p_2 + p_3 <= 0)
			new_p->rhs = p_1 -1;
		else
		{
			if(p_3 <= 0)
				new_p->rhs = p_1 + p_2;
			else
				new_p->rhs = p_1 + p_2 + p_3 + int(double(p_3)/2.0);
		}
		rhs_part = new_p->rhs;
		int indiceNZ = 0;
		for (int i3 = 0; i3 < m; i3++)
		{

			if (cutval_aux[i3] != 0)
			{
				cutval[indiceNZ] = cutval_aux[i3];
				cutind[indiceNZ] = cutind_aux[i3];

				//std::cout << cutval[indiceNZ] << " x_" << cutind[indiceNZ] << " + ";
				new_p->vectInd.push_back(cutind_aux[i3]);
				new_p->vectCoeff.push_back(cutval_aux[i3]);
				// new_ch->VarNameArray.push_back(varName);

				indiceNZ++;
			}
		}

		//std::cout << " >= " << rhs_part << std::endl;
		// voir si la coupe générée existe ou pas :
		if (!C_cut::CutAlreadyExists(new_p, C_cut::GeneratedPartition))
		{
			nbCtrViolee++;

			C_cut::GeneratedPartition.push_back(new_p);
			status = CPXcutcallbackadd(env, cbdata, wherefrom, indiceNZ, rhs_part, 'G', cutind, cutval, 1);
			if (status)
				cout << "Failed to add cut\n";
		}
		else
		{
			//cout << "Partition - Contrainte déjà obtenue" << endl;

			delete new_p;
		}

		delete[] cutind;
		delete[] cutval;
		delete[] cutval_aux;
		delete[] cutind_aux;
	}

	/*if (nbCtrViolee != 0)
		cout << "Nb Partition == " << nbCtrViolee << endl;
	*/
	return nbCtrViolee;
}

void ProceedToContractions(Graph *gr, int k)
{
	//*********************************************************************************
	//****  Si toutes les coupes sont vérifiées et que le point est fractionnaire  ****
	//****  alors on réduit le graphe et on sépare les nouvelles contraintes:      ****
	//****			partitions,F-Partitions,roues paires...                ****
	//*********************** Séparation nouvelles contraintes ************************
	//*********************************************************************************

	//*********  Réduction du graphe:opération Theta_2,3 et 4  *************
	b_Node *v_p, **bn_cour, *b_cour1;
	long nombre_contraction, n_contract;
	b_Node *contract_tab[N_SUR_2];
	Node *nptr;

	void (*operation[7])(Graph *, b_Node *, int, b_Node **, long *);
	long n_tab_contraction[6] = {0L, 0L, 0L, 0L, 0L, 0L}, nb_total_contraction;
	long n, m;

	n = gr->n_Nodes;
	m = gr->m_Edges;

	for (int i = 0; i < (n / 2) + 1; i++)
	{
		contract_tab[i] = NIL_BN;
	}

	// cout << "Réduction du graphe" << endl;

	operation[0] = operation_theta_3_bis; // Sans toucher aux arêtes fractionnaires
	operation[1] = operation_theta_3_bis;
	operation[2] = operation_theta_4_bis; // Recherche de clique
	operation[3] = operation_theta_4_4;	  // Contraction arêtes doubles
	operation[4] = operation_theta_4;	  // Theta_4 normal

	int nb_operation = 5;

	nb_total_contraction = 0;

	// PrintGraph_Table(gr);

	for (int i = 0; i < N_SUR_2; i++)
		contract_tab[i] = NIL_BN;

	for (int id_op = 0; id_op < nb_operation; id_op++)
	{
		nombre_contraction = 0;
		do
		{
			// Construction de l'ensemble V'
			if (id_op < 6)
			{
				v_p = NIL_BN;
				bn_cour = &v_p;
				nptr = &(gr->Nodes[0]);
				do
				{
					//(*bn_cour) = (b_Node *)malloc(sizeof(b_Node));
					(*bn_cour) = new b_Node;
					(*bn_cour)->id = nptr->id;
					(*bn_cour)->next = NIL_BN;

					bn_cour = &((*bn_cour)->next);
					nptr = nptr->next_gr;

				} while (nptr != &(gr->Nodes[0]));
			}

			// Démarrage de l'opération
			n_contract = 0;
 
			//Print_b_Node_Set(v_p, stdout);
			operation[id_op](gr, v_p, k, contract_tab, &n_contract); 

			n_tab_contraction[id_op] = n_tab_contraction[id_op] + n_contract;

			for (long i = 0; i < n_contract; i++)
			{ 
				//Print_b_Node_Set(contract_tab[i], stdout); 
				contract_set_w(gr, contract_tab[i]); 
				Delete_b_Node_Set(&contract_tab[i]);
			}

			nombre_contraction = n_contract;

			Delete_b_Node_Set(&v_p); 

		} while (nombre_contraction != 0);
	} 
}

void RemoveAllContractions(Graph *gr)
{

	int n = gr->n_Nodes;

	for (int i = 0; i < n - 1; i++)
	{

		gr->Nodes[i].next_gr = &(gr->Nodes[i + 1]);
		gr->Nodes[i].id_W = 1;
		gr->Nodes[i].next_W = &(gr->Nodes[i]);
		gr->Nodes[i].n_sh = &(gr->Nodes[i]);
		gr->Nodes[i].next_sh = &(gr->Nodes[i]);
	}

	gr->Nodes[n - 1].next_gr = &(gr->Nodes[0]);
	gr->Nodes[n - 1].id_W = 1;
	gr->Nodes[n - 1].next_W = &(gr->Nodes[n - 1]);
	gr->Nodes[n - 1].n_sh = &(gr->Nodes[n - 1]);
	gr->Nodes[n - 1].next_sh = &(gr->Nodes[n - 1]);

	save_contraction_info(gr, 1);
	save_contraction_info(gr, 2);
}



