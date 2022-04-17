#include <ilcplex/ilocplex.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <list>
#include <vector>
#include <algorithm>
#include <iterator>
#include <pthread.h>

#include ".././include/graph.h"
#include ".././include/Functions.h"
#include ".././include/Cuts.h"

using namespace std;

typedef IloArray<IloNumVarArray> NumVarMatrix;

vect mol;
// const double EpsForIntegrality=0.000001;

/***************************************************************************************************************************/

lpstrz tableau_20[] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "a"};
lpstrz tableau_25[] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "a"};
lpstrz tableau_30[] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "a"};
lpstrz tableau_40[] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "a"};
lpstrz tableau_50[] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "a"};

char **tableau = NULL;

const lpstrz *GetTableauLPSTRZ(int n)
{

	if (n == 20)
		return tableau_20;
	else if (n == 25)
		return tableau_25;
	else if (n == 30)
		return tableau_30;
	else if (n == 35)
		return tableau_30;
	else if (n == 40)
		return tableau_40;
	else if (n == 50)
		return tableau_50;
	else
		return tableau_50;
}

char **GetTableau(int n)
{

	if (tableau == NULL)
	{

		stringstream aSS;

		tableau = new char *[n + 1];

		for (int i = 0; i < n; i++)
		{

			aSS.clear();
			aSS.str("");
			// aSS << i;
			aSS << i;
			tableau[i] = new char[3];
			strcpy(tableau[i], aSS.str().c_str());
			// cout << "Tab[i]" << tableau[i] << endl;
		}

		tableau[n] = "a";
	}

	/*for(int i=0;i<n+1;i++){
		cout << tableau[i] << " ";
	}
	cout << endl;
*/
	return tableau;
}

/*****************************************	getInstanceName			********************************************/

string getInstanceName(char **argv)
{
	string arquivo(argv[1]);
	string::size_type loc = arquivo.find_last_of(".", arquivo.size());
	string::size_type loc2 = arquivo.find_last_of("/", arquivo.size());
	string nomeDaInstancia;
	if (loc != string::npos)
	{
		nomeDaInstancia.append(arquivo, loc2 + 1, loc - loc2 - 1);
		// cout << "\n1-" << nomeDaInstancia << endl;
	}
	else
	{
		nomeDaInstancia.append(arquivo, loc2 + 1, arquivo.size());
		// cout << "\n2-" << nomeDaInstancia << endl;
	}

	return nomeDaInstancia;
}
/***************************************************************************************************************************/

/*****************************************	getCollumName			********************************************/

string getCollumName(CPXCENVptr cpxEnv, CPXLPptr cpxModel, int collumIndex)
{
	int surplus;
	CPXgetcolname(cpxEnv, cpxModel, NULL, NULL, 0, &surplus, collumIndex, collumIndex);
	// Syntaxe int CPXgetcolname(CPXCENVptr env, CPXCLPptr lp, char ** name, char * namestore,
	//				int storespace, int * surplus_p, int begin, int end)
	//
	int storespace = -surplus;

	char *aux = new char[storespace];
	char **colname;
	colname = new char *[1];
	strcat(colname[0], "");

	int largestNameLength = 100;
	CPXgetcolname(cpxEnv, cpxModel, colname, aux, largestNameLength, &surplus, collumIndex, collumIndex);

	string result(colname[0]);

	delete[] colname;
	delete[] aux;
	return result;
}

/***************************************************************************************************************************/

/*****************************************	createLP	********************************************/

void createLP(C_graph G, lemon::ListGraph &leGraph)
{

	int i, j;
	int M = G.nb_nodes;
	int K = G.nb_edges;
	IloEnv env;
	IloModel model(env);
	IloRangeArray con(env);

	/*Creating variables*/

	double dist[K];
	IloNumVarArray x(env, K, 0, 1, ILOINT);

	char var[100];
	for (int i = 0; i < K; i++)
	{
		sprintf(var, "x%d(%d,%d)", G.Edges[i]->num, G.Edges[i]->end1->num, G.Edges[i]->end2->num);
		x[i].setName(var); 
	}
 
	/******************************************* Objective function ***********************************************/

	IloExpr obj(env);

	for (i = 0; i < K; i++)
	{
		// obj += dist[i]*x[i];
		obj += G.Edges[i]->length * x[i];
		// obj += x[i];
	}

	model.add(IloMinimize(env, obj));

	/************************************       creating constraints         **************************************/

	int nbCts1 = 0;

	for (int k = 0; k < M; k++)
	{

		IloExpr exp1(env);

		for (list<C_edge *>::iterator q = G.Nodes[k]->ed_incid.begin(); q != G.Nodes[k]->ed_incid.end(); ++q)
		{
			exp1 += x[(*q)->num]; 
		}

		IloRange r = (exp1 >= *std::min_element(G.con.begin(), G.con.end()));
		char a[100];
		sprintf(a, "a%d", nbCts1);
		r.setName(a);
		model.add(r);
		nbCts1++;
	}

	/*************************************************************************************************************/
	
	lemon::ListGraph::EdgeMap<double> weights(leGraph), capacities(leGraph), flows(leGraph), coeff(leGraph);
	for (lemon::ListGraph::EdgeIt j(leGraph); j!=lemon::INVALID; ++j)
	{
		weights[leGraph.edgeFromId(leGraph.id(j))] = 1.0;
		capacities[leGraph.edgeFromId(leGraph.id(j))] = 1.0; 
	}
	size_t nbCts2; 
	for (lemon::ListGraph::EdgeIt i(leGraph); i!=lemon::INVALID; ++i)
	{ 
		int linkID = leGraph.id(i);
		auto edge = leGraph.edgeFromId(linkID); 

		// ******************** Count number of links with link i ***************************
		int nb_links_with = 0;
		lemon::NetworkSimplex<lemon::ListGraph, double, double> ns(leGraph);
    	ns.costMap(weights).upperMap(capacities).stSupply(leGraph.u(edge), leGraph.v(edge), 2);
		auto status = ns.run(); 
		if(status == lemon::NetworkSimplex<lemon::ListGraph, double, double>::OPTIMAL) 
			nb_links_with = ns.totalCost() - 1;  

		// ******************** Count number of links without link i ***************************
		int nb_links_without = 0;
		capacities[edge] = 0.0; 
    	ns.costMap(weights).upperMap(capacities).stSupply(leGraph.u(edge), leGraph.v(edge), 2);
		status = ns.run();
		if(status == lemon::NetworkSimplex<lemon::ListGraph, double, double>::OPTIMAL) 
			nb_links_without = ns.totalCost() - 2;  

		if((nb_links_without - nb_links_with) < 2)
			continue;
			
		std::vector <size_t> incidentLinks;
		for(auto el : G.Nodes[leGraph.id(leGraph.u(edge))]->ed_incid)
			incidentLinks.push_back(el->num);
		for(auto el : G.Nodes[leGraph.id(leGraph.v(edge))]->ed_incid)
			incidentLinks.push_back(el->num); 

		IloExpr exp(env);
		for (lemon::ListGraph::EdgeIt j(leGraph); j!=lemon::INVALID; ++j)
			if(std::find(incidentLinks.begin(), incidentLinks.end(), leGraph.id(j)) == incidentLinks.end()) 
			{
				std::cout << "x_" << leGraph.id(j) << " + "; 
				exp += x[leGraph.id(j)]; 

			} 

		exp += (nb_links_without - nb_links_with) * x[linkID];
		std::cout <<  (nb_links_without - nb_links_with) << " * x_" << linkID << " >= " << nb_links_without << std::endl;

		IloRange r = (exp >= nb_links_without);
		char a[100];
		sprintf(a, "a%d", nbCts2);
		r.setName(a);
		model.add(r); 
		nbCts2++;		
	}

	IloCplex cplex(model);
	cplex.exportModel("dev.sav");
	cplex.exportModel("dev.lp");
	env.end();
}
/***************************************************************************************************************************/

/*****************************************	storeLPSolution			********************************************/

void storeLPSolution(CPXCENVptr env, CPXLPptr model, int numcols, double *x, double *sol_y, C_graph G, Graph *G_aux, double *objValue)
{
	char name[100];
	// cout<<"numcols: "<<numcols<<endl;

	(*objValue) = 0;

	for (int i1 = 0; i1 < numcols; i1++)
	{
		string s = getCollumName(env, model, i1);
		int i, j, k;
		if (s.at(0) == 'x')
		{
			strncpy(name, s.c_str(), 100);
			sscanf(name, "x%d(%d,%d)", &k, &i, &j);

			if (fabs(x[i1]) < EPSILON)
				sol_y[i1] = 0;
			else if (fabs(x[i1] - 1) < EPSILON)
				sol_y[i1] = 1;
			else
				sol_y[i1] = x[i1];

			G_aux->Edges[i1].X = G_aux->Edges[i1].back->X = sol_y[i1];

			(*objValue) += G_aux->Edges[i1].X * G_aux->Edges[i1].cap;

			// cout<<"y"<<k<<"("<<i+1<<","<<j+1<<") = "<<sol_y[i1]<<endl;
		}
	}
}

//     cout<<"before getx"<<endl;
//     CPXgetx (env, model, x, 0, numcols-1);
//     cout<<"After getx\n";
// //     for(int g=0; g<numcols;g++){cout<<x[g]<<endl;}
//     cout<<"StoreLPSolution starts\n";
// 	char name[100];
// 	for (int l = 0; l < numcols ; l++) {			//numcols : nbre de variables
// 		string s = getCollumName(env, model, l);
// 		int i,j,k;
// 		if (s.at(0) == 'y') {
// // 			cout<<s<<endl;
// 			strncpy(name,s.c_str(),100);
// 			// char * strncpy ( char * destination, const char * source, size_t num )
// 			// Copies the first num characters of source to destination
//
// 			sscanf(name,"y%d(%d,%d)",&k,&i,&j);
// 			sol_y[l] = x[l];
// 			cout<<"y"<<k<<"("<<i<<","<<j<<") = "<<sol_y[l]<<endl;
// // 			if(sol_y[l]!= 0)
// // 			{
// //  			 cout<<"sol_y["<<l<<"]		"<<name<<"	"<<sol_y[l]<<endl;
// 			  //cout<<"modf(sol_y[l], &intpart) "<<modf(sol_y[l], &intpart)<<endl;
// // 			}
// 		}
// 	}
//         cout<<"StoreLPSolution ends\n";

/***************************************************************************************************************************/

/*****************************************	printSolution		********************************************/

void printSolution(CPXCENVptr env, CPXLPptr model, int cur_numcols, double *sol_y)
{
	char varName[100];

	for (int i1 = 0; i1 < cur_numcols; i1++)
	{
		string s = getCollumName(env, model, i1);
		// cout<<s<<endl;
		int i, j, k;
		if (s.at(0) == 'x')
		{
			strncpy(varName, s.c_str(), 100);
			sscanf(varName, "x%d(%d,%d)", &k, &i, &j);
			// if (sol_y[l] > 0.01) {
			cout << "x" << k << "(" << i << "," << j << ") = " << sol_y[i1] << endl;
			//}
		}
	}
}
/***************************************************************************************************************************/

/*****************************************	printResults			********************************************/

void printResults(CPXCENVptr env, CPXLPptr model, string nomeDaInstancia, double time)
{

	int nodecount = CPXgetnodecnt(env, model);
	// Syntaxe int CPXgetnodecnt(CPXCENVptr env, CPXCLPptr lp)
	// The routine CPXgetnodecnt accesses the number of nodes used to solve a mixed integer problem

	double nodes_left = CPXgetnodeleftcnt(env, model);
	// Syntaxe int CPXgetnodeleftcnt(CPXCENVptr env, CPXCLPptr lp)
	// The routine CPXgetnodeleftcnt accesses the number of unexplored nodes left in the branch-and-cut tree

	bool optimum = false;
	int status = CPXgetstat(env, model);
	// Syntaxe int CPXgetstat(CPXCENVptr env, CPXCLPptr lp)
	// The routine CPXgetstat accesses the solution status of the problem after an LP, QP, QCP, or MIP optimization,
	//  after CPXfeasopt and its extensions, after CPXrefineconflict and its extensions

	if (status == 101 || status == 102)
	{
		optimum = true;
	}

	double UB, LB, Gap;

	CPXgetobjval(env, model, &UB);
	// Syntaxe int CPXgetobjval(CPXCENVptr env, CPXCLPptr lp, double * objval_p)
	// The routine CPXgetobjval accesses the solution objective value

	CPXgetbestobjval(env, model, &LB);
	// Syntaxe int CPXgetbestobjval(CPXCENVptr env, CPXCLPptr lp, double * objval_p)
	// The routine CPXgetbestobjval accesses the currently best known bound of all the
	//  remaining open nodes in a branch-and-cut tree.

	cout << "\n\nBranch-and-cut_Results: \n"
		 << endl;

	cout << "Instance: " << nomeDaInstancia << endl;
	cout << "Tree_Size: " << nodecount + nodes_left + 1 << endl;
	cout << "Total_Time: " << time << endl;

	cout << "LB: " << LB << endl;
	cout << "UB: " << UB << endl;
	cout << "Optimum: " << optimum << endl;
}
/***************************************************************************************************************************/

/*****************************************	printResultsToFile		********************************************/

void printResultsToFile(CPXCENVptr env, CPXLPptr model, string nameOfInstance, double time, int cur_numcols, C_graph G, int K)
{
	double *x = new double[cur_numcols];
	CPXgetx(env, model, x, 0, cur_numcols - 1);
	//        cout<<"After CPXgetx"<<endl;

	double *sol_y = new double[K];

	for (int i = 0; i < K; i++)
	{
		sol_y[i] = 0;
	}

	char varName[100];

	int nodecount = CPXgetnodecnt(env, model);
	double nodes_left = CPXgetnodeleftcnt(env, model);

	bool optimum = false;
	int status = CPXgetstat(env, model);
	if (status == 101 || status == 102)
	{
		optimum = true;
	}

	double UB, LB, gap;

	int numcovers;

	CPXgetnumcuts(env, model, CPX_CUT_USER, &numcovers);
	CPXgetobjval(env, model, &UB);
	CPXgetbestobjval(env, model, &LB);
	CPXgetmiprelgap(env, model, &gap);

	ofstream results("Results.txt", ios::app);
	results << "Instance: " << nameOfInstance << endl;
	results << " User Cuts Applied: " << numcovers << endl;
	results << " Tree_Size: " << nodecount + nodes_left + 1 << endl;
	results << " Gap: " << gap << endl;
	results << " Total_Time: " << time << " = " << ((int)time) / 60 << " min " << ((int)time) % 60 << " sec " << round((time - (int)time) * 100) << " ms" << endl;
	results << " LB: " << LB << endl;
	results << " UB: " << UB << endl;
	results << " Optimum: " << optimum << endl;
	results << "\n\nSolution:\n " << endl;

	for (int l = 0; l < cur_numcols; l++)
	{
		string s = getCollumName(env, model, l);

		int i, j, k;
		if (s.at(0) == 'x')
		{
			strncpy(varName, s.c_str(), 100);
			sscanf(varName, "x%d(%d,%d)", &k, &i, &j);
			sol_y[l] = x[l];
			if (sol_y[l] > 0.01)
			{
				results << "x(" << i << "," << j << ") " << sol_y[l] << endl;
			}
		}
	}
	results << endl;
	results.close();

	delete[] x;
	delete[] sol_y;
}

/*****************************************	isFrac			********************************************/

bool isFrac(int K, double *sol_y, double EpsForIntegrality)
{
	for (int j = 0; j < K; j++)
	{
		if (sol_y[j] != (floor(sol_y[j])))
		{
			return true;
		}
	}
	return false;
}

/******************************************		distEucl		*******************************************/

double distEucl(C_graph G, int i)
{
	return sqrt(pow(G.Edges[i]->end1->coord_x - G.Edges[i]->end2->coord_x, 2) + pow(G.Edges[i]->end1->coord_y - G.Edges[i]->end2->coord_y, 2));
}

/******************************************		belongs		*******************************************/

/*
bool belongs(char* V, int i, int k)
{

cout<<"entre dans belongs"<<endl;

for (int j = 0; j < k; j++) {

cout<<"j: "<<j<<endl;
char *pch;

char *a;
sprintf(a,"%d",i);

pch =strchr(V,j);
	printf ("pch: %c " , *pch);

pch = strpbrk (V, '1');

while (pch != NULL)
  {
	printf ("pch: %c " , *pch);
	pch = strpbrk (pch+1,1);
  }

if (pch!=NULL)
{
	return true;
}
}
return false;

}
*/
/*****************************************		fact		***************************************/

int fact(int n)
{
	int f = 1;
	for (int i = 1; i < n + 1; i++)
	{
		f = f * i;
	}
	return f;
}

/*****************************************		combinaison		***************************************/

double comb(int n, int k)
{

	double c;

	c = fact(n) / (fact(n - k) * fact(k));

	return c;
}

/*****************************************	vecteur combinaison		***************************************/

void combvec(int k, const lpstrz l[], const string &s, string &retour)
{

	if (k == 0)
	{

		if (retour == "")
		{
			retour = s;
		}

		else
		{

			retour = retour + "," + s;
		}
		return;
	}

	if (*l == "a")
		return;
	if (s.empty())
	{
		combvec(k - 1, l + 1, *l, retour);
	}
	else
	{

		combvec(k - 1, l + 1, s + "_" + *l, retour);
	};
	combvec(k, l + 1, s, retour);
}

/***************************************************************************************************************************************/
