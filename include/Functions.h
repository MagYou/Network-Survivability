#ifndef FUNCTIONS_H
#define FUNCTIONS_H

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

// using namespace std;

typedef char *lpstrz;

typedef vector<string> vect;

string getInstanceName(char **argv);

string getCollumName(CPXCENVptr cpxEnv, CPXLPptr cpxModel, int collumIndex);

void createLP(C_graph G);

void storeLPSolution(CPXCENVptr env, CPXLPptr model, int numcols, double *x, double *sol_y, C_graph G, Graph *G_aux, double *objValue);

void printSolution(CPXCENVptr env, CPXLPptr model, int cur_numcols, double *sol_y);

void printResults(CPXCENVptr env, CPXLPptr model, string nomeDaInstancia, double time);

void printResultsToFile(CPXCENVptr env, CPXLPptr model, string nameOfInstance, double time, int cur_numcols, C_graph G, int K);

bool isFrac(int K, double *sol_y, double EpsForIntegrality);

double distEucl(C_graph G, int i);

bool belongs(char *V, int i, int k);

int fact(int n);

double comb(int n, int k);

void combvec(int k, const lpstrz l[], const string &s, string &retour);

const lpstrz *GetTableauLPSTRZ(int n);
char **GetTableau(int n);

#endif // FUNCTIONS_H
