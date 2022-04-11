#include <ilcplex/ilocplex.h>

#include <lemon/list_graph.h>

#include <stdio.h>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string.h>
#include <list>
#include <ctime>
#include <iostream>
#include <string>

#include ".././include/graph.h"
#include ".././include/Functions.h"
#include ".././include/Cuts.h"

#include "../include/Graph_Plot.h"

/*****************************************************/

#include <lemon/gomory_hu.h>

/******************************************************/

const double EpsForIntegrality = 0.000001;

using namespace lemon;
using namespace std;

int main(int argc, char **argv)
{

	time_t tbegin, tend;
	double texec = 0.;
	tbegin = time(NULL);

	char *inputFileName;
	C_graph G;
	Graph G_aux;
	Graph_Plot G_plot; 
	simpleEdge *theBestSolution;
	double theBestValue;
 
	if (argc != 5)
	{
		cerr << "usage: " << argv[0] << " : Number of arguments is not correct." << endl;
		return 1;
	}
	inputFileName = new char[strlen(argv[1]) + 1];
	sprintf(inputFileName, "%s", argv[1]);
	string instanceName = getInstanceName(argv);
	cout << "\n\nThe name of the instance is: " << instanceName;
	
	G.read_instance(inputFileName, stoi(argv[2])); // true: sndlib, false: tsplib
 
	for(int n = 0; n < int(stof(argv[3])*G.nb_nodes); n++)
		G.con.push_back(1);
	for(int n = 0; n < int(stof(argv[4])*G.nb_nodes); n++)
		G.con.push_back(2);
	int start = G.con.size();
	for(int n = start; n < G.nb_nodes; n++)
		G.con.push_back(3);
	std::random_shuffle(G.con.begin(), G.con.end());

	std::copy(G.con.begin(), G.con.end(), std::back_inserter(G_aux.con));

	G.affiche();

	/*************************************************************************/
	/*Création du graphe auxilliaire pour les séparations de F-Partition, etc.*/
	/*************************************************************************/

	/*C=On crée la liste des arêtes du graphe au bon format*/
	simpleEdge *listeAreteAuxilliaire;

	theBestSolution = new simpleEdge[G.nb_edges];
	theBestValue = 0;

	for (int i = 0; i < G.nb_edges; i++)
	{
		theBestSolution[i].node1 = G.Edges[i]->end1->num + 1;
		theBestSolution[i].node2 = G.Edges[i]->end2->num + 1;
		theBestSolution[i].cap = G.Edges[i]->length;

		theBestValue += theBestSolution[i].cap;
	}

	if (!InitGraph_from_list(theBestSolution, G.nb_nodes, G.nb_edges, &G_aux))
	{
		cout << "ERROR:kECMASTER().Erreur pendant la lecture du graphe." << endl;
		exit(-1);
	}

	G_plot.setCoordinatesAndGraph(&G, &G_aux);

	/***********************************    B & C                       **********************************************/

	int K = G.nb_edges;
	createLP(G);
 
	// Chargement du modèle
	CPXENVptr env = NULL;
	CPXLPptr model = NULL;
	int status = 0;

	env = CPXopenCPLEX(&status); // If the routine is successful, then status is 0
								 //    if(status==0)
								 //    cout<<"chargement de l'environnement OK !\n";

	model = CPXcreateprob(env, &status, "dev"); // If the routine is successful, then status is 0. The problem that is created is an LP minimization problem with zero constraints, zero variables, and an empty constraint matrix
												//    if(status==0)
												//    cout<<"chargement du modèle OK !\n";
	// the problem is read from a file
	status = CPXreadcopyprob(env, model, "dev.sav", NULL); 

	int *index = new int[CPXgetnumcols(env, model)];
	for (int i = 0; i < CPXgetnumcols(env, model); i++)
	{
		index[i] = i;
	}
	char *xctype = new char[CPXgetnumcols(env, model)];
	for (int i = 0; i < CPXgetnumcols(env, model); i++)
	{
		xctype[i] = 'I';
	}
	status = CPXchgctype(env, model, CPXgetnumcols(env, model), index, xctype);
	status = CPXchgprobtype(env, model, 1); // to change to MILP
	int cur_numcols = CPXgetnumcols(env, model);

	/****************************************************************************************************************************/

	CUTINFO cutinfo;
	updateCutinfo(&cutinfo, model, cur_numcols, G, &G_aux, &G_plot, theBestSolution, &theBestValue, K, inputFileName);

	LAZYCUTINFO lazycutinfo;
	updateLazyCutinfo(&lazycutinfo, model, cur_numcols, G, &G_aux, &G_plot, theBestSolution, &theBestValue, K, inputFileName);

	// CPLEX output on the screen
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);

	CPXsetintparam(env, CPX_PARAM_MIPINTERVAL, 100);

	// This is very important, the callbacks do not work without this
	// Ensure linear mappings between the presolved and original models
	CPXsetintparam(env, CPX_PARAM_PRELINEAR, 0);

	// Let MIP callbacks work on the original model
	CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
	CPXsetintparam(env, CPX_PARAM_REDUCE, CPX_PREREDUCE_PRIMALONLY);
	// CPXsetintparam (env, CPX_PARAM_REDUCE, CPX_PREREDUCE_NOPRIMALORDUAL);
	// In order to use more threads, it is necesseray to deal with some extra issues. Better use 1 by now
	// In this context, maximum number of threads means the number of CPUs or cores available.
	CPXsetintparam(env, CPX_PARAM_THREADS, 1);
	CPXsetintparam(env, CPX_PARAM_PARALLELMODE, 1);
	CPXsetintparam(env, CPX_PARAM_VARSEL, 3);

	// CPXsetintparam (env, CPX_PARAM_TILIM, 18000.0);
	CPXsetdblparam(env, CPX_PARAM_TILIM, 18000.0);

	// CPXsetintparam (env, CPX_PARAM_MIPSEARCH, 1);
	CPXsetintparam(env, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL);

	// turn off presolve and heuristics
	CPXsetintparam(env, CPX_PARAM_PREIND, false);
	CPXsetintparam(env, CPX_PARAM_HEURFREQ, -1);

	CPXsetintparam(env, CPX_PARAM_DATACHECK, CPX_ON);

	CPXsetintparam(env, CPX_PARAM_CLIQUES, -1);
	CPXsetintparam(env, CPX_PARAM_COVERS, -1);
	CPXsetintparam(env, CPX_PARAM_DISJCUTS, -1);
	CPXsetintparam(env, CPX_PARAM_FLOWCOVERS, -1);
	CPXsetintparam(env, CPX_PARAM_FLOWPATHS, -1);
	CPXsetintparam(env, CPX_PARAM_FRACCUTS, -1);
	CPXsetintparam(env, CPX_PARAM_GUBCOVERS, -1);
	CPXsetintparam(env, CPX_PARAM_IMPLBD, -1);
	CPXsetintparam(env, CPX_PARAM_MIRCUTS, -1);
	CPXsetintparam(env, CPX_PARAM_ZEROHALFCUTS, -1);

	CPXsetintparam(env, CPX_PARAM_MIPDISPLAY, 2);

	CPXsetdblparam(env, CPX_PARAM_TILIM, 7200.0);

	CPXsetusercutcallbackfunc(env, mycutcallback, &cutinfo);
	CPXsetlazyconstraintcallbackfunc(env, mycutcallback, &lazycutinfo);

	// cout<<"after brunch and cut"<<endl;

	double cplexTimeBefore, cplexTimeAfter;
	CPXgettime(env, &cplexTimeBefore);
	status = CPXmipopt(env, model); // if(status==CPXERR_NO_MEMORY){cout<<"no memory\b";}if(status==CPXERR_NO_PROBLEM){cout<<"no memory\b";}
	// cout<<status<<endl;

	CPXgetstat(env, model);

	CPXgettime(env, &cplexTimeAfter);
	double time1 = cplexTimeAfter - cplexTimeBefore;
 

	double UB;

	CPXgetobjval(env, model, &UB);
	cout << UB << endl;

	// printSolution(env, model,cur_numcols, K, M);
	printResultsToFile(env, model, instanceName, time1, cur_numcols - 1, G, K);
	printResults(env, model, instanceName, time1);

	CPXfreeprob(env, &model);
	CPXcloseCPLEX(&env);

	cout << "\n"
		 << endl;

	delete[] index;
	delete[] xctype;
	delete[] inputFileName;
	tend = time(NULL);
	texec = difftime(tend, tbegin);
	cout << "CPU Time: " << tend - tbegin << endl;
	cout << "Nb Edge Cut:" << C_cut::GeneratedEdgeCuts.size() << endl;
	cout << "Nb SP-Partition:" << C_cut::GeneratedSPPartition.size() << endl;
	cout << "Nb F-Partition:" << C_cut::GeneratedFPartition.size() << endl;
	cout << "Nb Partition:" << C_cut::GeneratedPartition.size() << endl;

	/****************************************************************************************************************************/

	return 0;
} // end int main
