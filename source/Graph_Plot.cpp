/*
 * Graph_Plot.cpp
 *
 *  Created on: 6 déc. 2014
 *      Author: myriam
 */

#include ".././include/Graph_Plot.h"
#include <iostream>
#include <iomanip>

void Graph_Plot::setCoordinatesAndGraph(C_graph *Gr, Graph *G_aux)
{
	g = G_aux;

	x_coord = new double[Gr->nb_nodes];
	y_coord = new double[Gr->nb_nodes];

	for (int i = 0; i < Gr->nb_nodes; i++)
	{
		x_coord[i] = Gr->Nodes[i]->coord_x;
		y_coord[i] = Gr->Nodes[i]->coord_y;
	}
}

void Graph_Plot::save_xfig(const char *xfig_file, bool avec_val, int epais)
{
	ofstream fich(xfig_file);

	/*En-tête du fichier*/
	fich << "#FIG 3.2" << endl;
	fich << "Landscape" << endl;
	fich << "Center" << endl;
	fich << "Metric" << endl;
	fich << "A4" << endl;
	fich << "100.00" << endl;
	fich << "Single" << endl;
	fich << "-2" << endl;
	fich << "1200 2" << endl;

	// calcul du max des coordonnees pour l'echelle
	double max_x = x_coord[0];
	double max_y = y_coord[0];
	double min_x = x_coord[0];
	double min_y = y_coord[0];

	for (int i = 0; i < g->n_Nodes; i++)
	{
		if (x_coord[i] > max_x)
			max_x = x_coord[i];

		if (y_coord[i] > max_y)
			max_y = y_coord[i];

		if (x_coord[i] < min_x)
			min_x = x_coord[i];

		if (y_coord[i] < min_y)
			min_y = y_coord[i];
	}

	// affichage des sommets
	for (int i = 0; i < g->n_Nodes; i++)
	{
		save_xfig_node(fich, x_coord[i], y_coord[i], 0, epais, 0, max_x, max_y, min_x, min_y);
		save_xfig_node_name(fich, x_coord[i], y_coord[i], i + 1, 0, epais, max_x, max_y, min_x, min_y);

		/*G.vect_sommets[i]->sauve_xfig(fich,0,epais,0,max_x,max_y,min_x,min_y);
		G.vect_sommets[i]->sauve_xfig_nomme_sommet(fich,0,epais,max_x,max_y,min_x,min_y);*/
	}

	// affichage des aretes
	for (int i = 0; i < g->m_Edges; i++)
	{
		if (g->Edges[i].X == 1.0)
			save_xfig_edge(fich, g->Edges[i].back->adjac->id, g->Edges[i].adjac->id, 0, epais, 0, max_x, max_y, min_x, min_y);
		else
		{
			if (g->Edges[i].X > 0) /*Les aretes fractionnaires sont en rouge et pointilles*/
			{
				save_xfig_edge(fich, g->Edges[i].back->adjac->id, g->Edges[i].adjac->id, 4, epais, 1, max_x, max_y, min_x, min_y);

				if (avec_val == true)
					save_xfig_valuation(fich, g->Edges[i].back->adjac->id, g->Edges[i].adjac->id, g->Edges[i].X, 4, epais, 0, max_x, max_y, min_x, min_y);
			}
		}
	}

	fich.close();
}

void Graph_Plot::save_xfig(const char *xfig_file, simpleEdge *solution, bool avec_val, int epais)
{
	ofstream fich(xfig_file);

	/*En-tête du fichier*/
	fich << "#FIG 3.2" << endl;
	fich << "Landscape" << endl;
	fich << "Center" << endl;
	fich << "Metric" << endl;
	fich << "A4" << endl;
	fich << "100.00" << endl;
	fich << "Single" << endl;
	fich << "-2" << endl;
	fich << "1200 2" << endl;

	// calcul du max des coordonnees pour l'echelle
	double max_x = x_coord[0];
	double max_y = y_coord[0];
	double min_x = x_coord[0];
	double min_y = y_coord[0];

	for (int i = 0; i < g->n_Nodes; i++)
	{
		if (x_coord[i] > max_x)
			max_x = x_coord[i];

		if (y_coord[i] > max_y)
			max_y = y_coord[i];

		if (x_coord[i] < min_x)
			min_x = x_coord[i];

		if (y_coord[i] < min_y)
			min_y = y_coord[i];
	}

	// affichage des sommets
	for (int i = 0; i < g->n_Nodes; i++)
	{
		save_xfig_node(fich, x_coord[i], y_coord[i], 0, epais, 0, max_x, max_y, min_x, min_y);
		save_xfig_node_name(fich, x_coord[i], y_coord[i], i + 1, 0, epais, max_x, max_y, min_x, min_y);

		/*G.vect_sommets[i]->sauve_xfig(fich,0,epais,0,max_x,max_y,min_x,min_y);
		G.vect_sommets[i]->sauve_xfig_nomme_sommet(fich,0,epais,max_x,max_y,min_x,min_y);*/
	}

	// affichage des aretes
	for (int i = 0; i < g->m_Edges; i++)
	{
		if (solution[i].cap == 1.0)
			save_xfig_edge(fich, solution[i].node1, solution[i].node2, 0, epais, 0, max_x, max_y, min_x, min_y);
		else
		{
			if (solution[i].cap > 0) /*Les aretes fractionnaires sont en rouge et pointilles*/
			{
				save_xfig_edge(fich, solution[i].node1, solution[i].node2, 4, epais, 1, max_x, max_y, min_x, min_y);

				if (avec_val == true)
					save_xfig_valuation(fich, solution[i].node1, solution[i].node2, solution[i].cap, 4, epais, 0, max_x, max_y, min_x, min_y);
			}
		}
	}

	fich.close();
}

void Graph_Plot::save_xfig_node(ostream &fic, double node_x, double node_y, int couleur, int epais, int forme, double lx, double ly, double dx, double dy)
{
	int pas = 40 * epais;
	int x2 = (int)(pas * 0.7071);

	dx = dx - 40;
	dy = dy - 40;

	double echelle = min((12500 / (lx - dx)), (8500 / (ly - dy)));

	if (echelle > 10)
		echelle = max((12500 / (lx - dx)), (8500 / (ly - dy)));

	int x = (int)(node_x * echelle);
	int y = (int)(node_y * echelle);

	int h = (int)(0.87 * pas);

	/*Pour un point*/
	if (forme == 0)
	{
		fic << "1 4 0 0 0 " << couleur << " 10 0 20 0.000 1 0.0000 " << x << " " << y << " " << pas << " " << pas << " 2047 967 2273 967" << endl;
	}

	/*Pour une ellipse*/
	if (forme == 1)
	{
		fic << "2 2 0 1 " << couleur << " " << couleur << " 10 0 20 0.000 0 0 -1 0 0 5";
		fic << endl;
		fic << "          " << x - h << " " << y - h << " " << x + h << " " << y - h << " ";
		fic << x + h << " " << y + h << " " << x - h << " " << y + h << " " << x - h << " " << y - h;
		fic << endl;
	}

	/*Pour une ligne*/
	if (forme == 2)
	{
		fic << "1 4 0 1 " << couleur << " 7 10 0 20 0.000 1 0.0000 " << x << " " << y << " " << pas << " " << pas << " 2047 967 2273 967" << endl;
		fic << "2 1 0 1 " << couleur << " 7 9 0 -1 0.000 0 0 -1 0 0 2" << endl;
		fic << " " << x + x2 << " " << y + x2 << " " << x - x2 << " " << y - x2 << endl;
		fic << "2 1 0 1 " << couleur << " 7 9 0 -1 0.000 0 0 -1 0 0 2" << endl;
		fic << " " << x + x2 << " " << y - x2 << " " << x - x2 << " " << y + x2 << endl;
	}
}

void Graph_Plot::save_xfig_node_name(ostream &fic, double node_x, double node_y, int num, int couleur, int epais, double lx, double ly, double dx, double dy)
{
	int pas = 60 * epais;

	dx = dx - 40;
	dy = dy - 40;

	double echelle = min((12500 / (lx - dx)), (8500 / (ly - dy)));

	if (echelle > 10)
		echelle = max((12500 / (lx - dx)), (8500 / (ly - dy)));

	int x = (int)(node_x * echelle);
	int y = (int)(node_y * echelle);

	fic << "4 0 0 100 " << 0 << " 14 " << 9 + epais << " 0.0000 4 90 90 " << (int)(x - (int)pas / 2) << " " << (int)(y - pas) << " " << num << "\\001" << endl;
}

void Graph_Plot::save_xfig_edge(ostream &fic, long node1, long node2, int couleur, int epais, int forme, double lx, double ly, double dx, double dy)
{
	// gestion de l'echelle
	dx = dx - 40;
	dy = dy - 40;

	double echelle = min((12500 / (lx - dx)), (8500 / (ly - dy)));

	if (echelle > 10)
		echelle = max((12500 / (lx - dx)), (8500 / (ly - dy)));

	int x1 = (int)((x_coord[node1 - 1]) * echelle);
	int y1 = (int)((y_coord[node1 - 1]) * echelle);
	int x2 = (int)((x_coord[node2 - 1]) * echelle);
	int y2 = (int)((y_coord[node2 - 1]) * echelle);

	if (forme == 0)
	{
		fic << "2 1 0 " << epais << " " << couleur << " 7 100 0 -1 0.000 0 0 -1 0 0 2";
		fic << endl;
		fic << " " << x1 << " " << y1 << " " << x2 << " " << y2 << endl;
	}

	if (forme == 1)
	{
		fic << "2 1 1 " << epais << " " << couleur << " 7 100 0 -1 6.000 0 0 -1 0 0 2";
		fic << endl;
		fic << " " << x1 << " " << y1 << " " << x2 << " " << y2 << endl;
	}

	if (forme == 2)
	{
		fic << "2 1 2 " << epais << " " << couleur << " 7 100 0 -1 1.500 0 0 -1 0 0 2";
		fic << endl;
		fic << " " << x1 << " " << y1 << " " << x2 << " " << y2 << endl;
	}
}

void Graph_Plot::save_xfig_valuation(ostream &fic, long node1, long node2, double valuation, int couleur, int epais, int forme, double lx, double ly, double dx, double dy)
{
	int pas = 60 * epais;

	// gestion de l'echelle
	dx = dx - 40;
	dy = dy - 40;

	double echelle = min((12500 / (lx - dx)), (8500 / (ly - dy)));

	if (echelle > 10)
		echelle = max((12500 / (lx - dx)), (8500 / (ly - dy)));

	int x1 = (int)((x_coord[node1 - 1]) * echelle);
	int y1 = (int)((y_coord[node1 - 1]) * echelle);
	int x2 = (int)((x_coord[node2 - 1]) * echelle);
	int y2 = (int)((y_coord[node2 - 1]) * echelle);

	int x = (int)((x1 + x2) / 2);
	int y = (int)((y1 + y2) / 2);

	fic << "4 0 0 100 " << 0 << " 14 " << 7 + epais << " 0.0000 4 90 90 " << (int)(x - (int)pas / 2) << " " << (int)(y - pas) << " " << setprecision(2) << valuation << "\\001" << endl;
}

void Graph_Plot::save_shrunk_xfig(const char *xfig_file, simpleEdge *sh_List, long sh_edge, bool avec_val, int epais)
{
	Graph *gr = g;

	ofstream fich(xfig_file);

	/*En-tête du fichier*/
	fich << "#FIG 3.2" << endl;
	fich << "Landscape" << endl;
	fich << "Center" << endl;
	fich << "Metric" << endl;
	fich << "A4" << endl;
	fich << "100.00" << endl;
	fich << "Single" << endl;
	fich << "-2" << endl;
	fich << "1200 2" << endl;

	/*Calcul de la liste des arêtes du graphe réduit*/
	// long sh_edge;
	// simpleEdge *sh_List = GetReducedGraph_eList(&gr,&sh_edge);

	/*calcul du min et max des coordonnees pour l'echelle*/
	double max_x = x_coord[0];
	double max_y = y_coord[0];
	double min_x = x_coord[0];
	double min_y = y_coord[0];

	/*Le sommet 1 est toujours dans le graphe réduit*/
	long first = 1L;
	long i = 1L;
	do
	{
		if (x_coord[i - 1] > max_x)
			max_x = x_coord[i - 1];

		if (y_coord[i - 1] > max_y)
			max_y = y_coord[i - 1];

		if (x_coord[i - 1] < min_x)
			min_x = x_coord[i - 1];

		if (y_coord[i - 1] < min_y)
			min_y = y_coord[i - 1];

		i = gr->Nodes[i - 1].next_gr->id;
	} while (i != first);

	// affichage des sommets
	first = 1L;
	i = 1L;
	do
	{
		save_xfig_node(fich, x_coord[i - 1], y_coord[i - 1], 0, epais, 0, max_x, max_y, min_x, min_y);
		save_xfig_sh_node_name(fich, x_coord[i - 1], y_coord[i - 1], i, 0, epais, max_x, max_y, min_x, min_y);

		i = gr->Nodes[i - 1].next_gr->id;
	} while (i != first);

	// affichage des aretes
	for (int i = 0; i < sh_edge; i++)
	{
		if (sh_List[i].cap == 1.0)
			save_xfig_edge(fich, sh_List[i].node1, sh_List[i].node2, 0, epais, 0, max_x, max_y, min_x, min_y);
		else
		{
			if (sh_List[i].cap > 0) /*Les aretes fractionnaires sont en rouge et pointilles*/
			{
				save_xfig_edge(fich, sh_List[i].node1, sh_List[i].node2, 4, epais, 1, max_x, max_y, min_x, min_y);

				if (avec_val == true)
					save_xfig_valuation(fich, sh_List[i].node1, sh_List[i].node2, sh_List[i].cap, 4, epais, 0, max_x, max_y, min_x, min_y);
			}
		}
	}

	fich.close();
}

void Graph_Plot::save_xfig_sh_node_name(ostream &fic, double node_x, double node_y, int num, int couleur, int epais, double lx, double ly, double dx, double dy)
{
	Graph *gr = g;

	int pas = 60 * epais;

	dx = dx - 40;
	dy = dy - 40;

	double echelle = min((12500 / (lx - dx)), (8500 / (ly - dy)));

	if (echelle > 10)
		echelle = max((12500 / (lx - dx)), (8500 / (ly - dy)));

	int x = (int)(node_x * echelle);
	int y = (int)(node_y * echelle);

	/*On affiche (1,2,5,7,13) où 2,5,7,13 sont contractés sur le sommet 1*/
	string sh_node_label = "";

	/*Parcours de la liste de contraction du sommet num*/
	long n_contr = num;
	char s[256] = "";

	do
	{
		sprintf(s, ",%ld\0", n_contr);
		sh_node_label = s + sh_node_label;

		n_contr = gr->Nodes[n_contr - 1].next_sh->id;
	} while (n_contr != num);

	sh_node_label = "(" + sh_node_label.substr(1) + ")";

	fic << "4 0 0 100 " << 0 << " 14 " << 9 + epais << " 0.0000 4 90 90 " << (int)(x - (int)pas / 2) << " " << (int)(y - pas) << " " << sh_node_label << "\\001" << endl;
}
