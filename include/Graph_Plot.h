/*
 * Graph_Plot.h
 *
 *  Created on: 6 déc. 2014
 *      Author: myriam
 */

#ifndef GRAPH_PLOT_H_
#define GRAPH_PLOT_H_

#include ".././include/graph.h"

class C_graph;

class Graph_Plot
{

protected:
	Graph *g;
	double *x_coord;
	double *y_coord;

public:
	void setCoordinatesAndGraph(C_graph *Gr, Graph *G_aux);

	/*Methode de sauvegarde d'un point dans un fichier XFIG*/
	/*La methode ne traite que les fichiers de la TSPLib.*/
	/*Le point est donné par le tableau sol*/
	void save_xfig(const char *xfig_file, bool avec_val, int epais);
	void save_xfig(const char *xfig_file, simpleEdge *sol, bool avec_val, int epais);

	void save_shrunk_xfig(const char *xfig_file, simpleEdge *sh_eList, long sh_edge, bool avec_val, int epais);

protected:
	/*La forme correspond au style de trait (plein, tirets,...)*/
	void save_xfig_node(ostream &fic, double node_x, double node_y, int couleur, int epais, int forme, double lx, double ly, double dx, double dy);

	void save_xfig_node_name(ostream &fic, double node_x, double node_y, int num, int couleur, int epais, double lx, double ly, double dx, double dy);

	void save_xfig_edge(ostream &fic, long node1, long node2, int couleur, int epais, int forme, double lx, double ly, double dx, double dy);

	/*Ajoute la valuation (numero ou capacite,...) sur l'arete*/
	void save_xfig_valuation(ostream &fic, long node1, long node2, double valuation, int couleur, int epais, int forme, double lx, double ly, double dx, double dy);

	void save_xfig_sh_node_name(ostream &fic, double node_x, double node_y, int num, int couleur, int epais, double lx, double ly, double dx, double dy);
};

#endif /* _GRAPH_PLOT_H_ */