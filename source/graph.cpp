//***********************//
//*** FILE graph.cpp ***//
//*********************//

/*This file contains the classes required to describe the graph : Node, Arc, Demand, Graph, ...*/
/*It also contains the functions that read the instances*/

#include ".././include/graph.h"

#include <math.h>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>

Bool *getMark_Tab()
{

	static Bool mark_tab[N_TAB];
	return mark_tab;
}

static Bool mark_tab[N_TAB];

/*************************************************************************************************/
//					***	Classe des Noeuds	***
/*************************************************************************************************/
C_node::C_node()
{
	name = new char[60];
}

C_node::~C_node()
{
	delete[] name;
}
/*************************************************************************************************/
//					***	Classe des Arêtes	***
/*************************************************************************************************/
C_edge::C_edge()
{
}

C_edge::C_edge(C_node *v1, C_node *v2)
{
	end1 = v1;
	end2 = v2;
}

C_edge::~C_edge()
{
}

C_node *C_edge::get_end1()
{
	return end1;
}

C_node *C_edge::get_end2()
{
	return end2;
}

double C_edge::dist_euclid(C_node *v1, C_node *v2)
{
	double dist = 0;

	double d1, d2;

	d1 = pow(v1->coord_x - v2->coord_x, 2);
	d2 = pow(v1->coord_y - v2->coord_y, 2);
	dist = ceil(sqrt(d1 + d2));

	return dist;
}

/*************************************************************************************************/
//					***	Classe des Demandes	***
/*************************************************************************************************/
/*
C_demand::C_demand()
{}

C_demand::~C_demand()
{
}
*/
/*************************************************************************************************/
//					***	Classe des Graphes	***
/*************************************************************************************************/

C_graph::C_graph()
{
	nb_nodes = 0;
	nb_edges = 0;
}

C_graph::~C_graph()
{
}

//***********************************	Fonction read_instance	**********************************/
void C_graph::read_instance(const char *name, bool sndlib_tsplib)
{

	char *nom_ext;
	char *buffer;

	C_node *ptr_node;
	C_edge *ptr_edge;
	// C_demand *ptr_demand;

	bool open = false;
	bool dimensions_found = false;
	bool nodes_found = false;
	// bool edges_found=false;
	// bool demands_found=false;

	int i, j, k;

	nom_ext = new char[100];
	sprintf(nom_ext, "%s", name);

	ifstream fic(nom_ext);
	buffer = new char[1050];
	if (!(fic))
	{
		cout << "Fail to read the file of random instance" << endl;
		cout << "The file of random instance " << nom_ext << " could not be open" << endl;
	}
	else // file open
	{
		open = true;

		/*Reading the dimensions*/

		fic >> buffer;
		while (strcmp(buffer, "Dimensions") != 0)
			fic >> buffer;
		if (strcmp(buffer, "Dimensions") == 0)
		{
			dimensions_found = true; // dimensions were found
			fic >> buffer;
			fic >> nb_nodes;
			fic >> buffer;
		} // end if 293
		else
			dimensions_found = false; // dimensions were not found

		if (!dimensions_found)
			fic >> buffer;

		/************************************************ Reading nodes***********************************************/
		while (strcmp(buffer, "Nodes") != 0)
			fic >> buffer;
		if (strcmp(buffer, "Nodes") == 0)
		{
			nodes_found = true; // nodes section was found
		}						// end if 326
		else
			nodes_found = false;

		if (!nodes_found)
			cout << "nodes section was not found" << endl;
		cout << endl;
		if (open && dimensions_found && nodes_found)
		{
			Nodes.resize(nb_nodes);
			for (i = 0; i < nb_nodes; i++)
			{
				// cout<<"node: "<<i<<endl;
				ptr_node = new (C_node);

				if (sndlib_tsplib)
				{
					/*****sndlib*****/

					fic >> ptr_node->num;
					fic >> ptr_node->name;
					/******************/
				}
				else
				{
					/*******tsplib*********/
					ptr_node->num = i;
					ptr_node->name = "node";
					/************************/
				}
				fic >> ptr_node->coord_x;
				fic >> ptr_node->coord_y;
				// fic>>ptr_node->connexite;
				Nodes[i] = ptr_node;
			} // end for
		}	  // end if 337

		/******************************************Building edges************************************************/

		nb_edges = nb_nodes * (nb_nodes - 1) / 2;

		for (i = 0; i < nb_nodes; i++)
		{

			Nodes[i]->ed_incid.erase(Nodes[i]->ed_incid.begin(), Nodes[i]->ed_incid.end());
		}

		// cout<<"nb_edges: "<<nb_edges<<endl;
		Edges.resize(nb_edges);
		int m = 0;
		for (int i = 0; i < nb_nodes - 1; i++)
		{
			for (int j = i + 1; j < nb_nodes; j++)
			{
				// cout<<"edge: "<<m<<endl;
				ptr_edge = new (C_edge);
				ptr_edge->num = m;
				ptr_edge->end1 = Nodes[i];
				ptr_edge->end2 = Nodes[j];
				ptr_edge->length = ptr_edge->dist_euclid(Nodes[i], Nodes[j]);
				Edges[m] = ptr_edge;

				Nodes[i]->ed_incid.push_back(ptr_edge);
				Nodes[j]->ed_incid.push_back(ptr_edge);

				// Nodes[i]->nd_incid.push_back(Nodes[j]);
				// Nodes[j]->nd_incid.push_back(Nodes[i]);

				m = m + 1;
			}
		}

		/********************************************** Reading edges ************************************************/
		/*
			  fic>>buffer;
			  while (strcmp(buffer,"Edges")!=0)
			{
			  //Edges section not found
			  fic>>buffer;
			}//end while 356
			  if(strcmp(buffer,"Edges")==0)
			{
			  edges_found=true;
			} //end if 361
			  else
			edges_found=false;
			  if(!edges_found)
			cout<<"Edges were not found"<<endl;
			  if(open && dimensions_found && edges_found)
			{
			  Edges.resize(nb_edges);
			  for(i=0;i<nb_edges;i++)
				{
					  ptr_edge=new(C_edge);
				  ptr_edge->num=i;
				  fic>>buffer;
				  ptr_edge->end1=Nodes[atoi(buffer)];// atoi=ascii to integer;
				  fic>>buffer;
				  ptr_edge->end2=Nodes[atoi(buffer)];
				  fic>>buffer;
				  ptr_edge->length=atof(buffer);
					  Edges[i]=ptr_edge;
				  }// end if 375
				} //end if 371
		*/
		/************************************** Reading the demands **************************************************/
		/*
		 fic>>buffer;
			  while (strcmp(buffer,"Demands")!=0)
			{
			  fic>>buffer;
				}// end while 392
			  if (strcmp(buffer,"Demands")==0)
			{
			  demands_found=true;
			} //end if 396
			  else
			  demands_found=false;

			  if (!demands_found)
			cout<<"Demands were not found"<<endl;

			  if (open && dimensions_found && demands_found)
			{
			  Demands.resize(nb_demands);

			  for(i=0;i<nb_demands;i++)
				{

					  ptr_demand=new(C_demand);
				  ptr_demand->num=i;
				  fic>>buffer;
				  ptr_demand->source=Nodes[atoi(buffer)];
				  fic>>buffer;
				  ptr_demand->destination=Nodes[atoi(buffer)];
				  fic>>buffer;
					  ptr_demand->traffic=atof(buffer);
					  Demands[i]=ptr_demand;


				}// end for 411


			 } //end if 406
		*/
		fic.close();

		delete[] nom_ext;
		delete[] buffer;
	}
} // end read_instance

//**************************	Fonction affiche    ***************************//
void C_graph::affiche()
{

	cout << "The graph of the instance is defined by:" << endl;
	// cout<<endl;
	cout << "	" << nb_nodes << " nodes" << endl;
	cout << "	" << nb_edges << " edges" << endl; 
}

Bool InitGraph_from_file(char *fich, Graph *gr, int *k)
{
	long n, m, i, j;
	char c[2];

	Edge *eptr1, *eptr2;
	long node1, node2;
	double cap;
	Node *nptr1, *nptr2;

	FILE *f = fopen(fich, "r");

	if (f == NULL)
	{
		fprintf(stdout, "%s %s\n", "ERROR:InitGraph_from_file().Unable to open file == ", fich);

		return kECSP_False;
	}

	/*Reading of the order of the problem*/
	fscanf(f, "%s %d", c, k);
	/*printf("k = %d\n",(*k));*/

	/*Reading nodes and edges numbers*/
	fscanf(f, "%s", c);

	if (c[0] != 'p')
		return kECSP_False;

	fscanf(f, "%ld", &n);
	fscanf(f, "%ld", &m);

	/*fscanf(f,"%d",k);
	fscanf(f,"%ld",&n);
	fscanf(f,"%ld",&m);*/

	if (!AllocateGraph(n, m, gr))
	{
		return kECSP_False;
	}

	/*Reading of graph edges*/
	eptr1 = &(gr->Edges[0L]);
	eptr2 = &(gr->Edges[m]);

	for (j = 0L; j < m; j++)
	{
		if (fscanf(f, "%s %ld %ld %lf", c, &node1, &node2, &cap) == EOF)
		{
			fprintf(stderr, "EOF reached in input file, %ld edges read\n", j);
			DeleteGraph(gr);

			return kECSP_False;
		}

		/*/printf("%s %ld %ld %f\n",c,node1,node2,cap);*/

		if ((c[0] != 'e') || (node1 < 1) || (node1 > n) || (node2 < 1) || (node2 > n))
		{
			puts("Error in the input file");
			DeleteGraph(gr);
			return kECSP_False;
		}

		--node1;
		--node2;
		nptr1 = &(gr->Nodes[node1]);
		nptr2 = &(gr->Nodes[node2]);
		eptr1->adjac = nptr2;
		eptr2->adjac = nptr1;
		eptr1->cap = cap;
		eptr1->X = cap;
		/*eptr1->X = 0.0;*/
		eptr2->cap = cap;
		eptr2->X = cap;
		/*eptr2->X = 0.0;*/
		eptr1->back = eptr2;
		eptr2->back = eptr1;

		eptr1->num = j + 1;
		eptr2->num = j + 1;

		if (nptr1->first_edge == NULL)
		{
			nptr1->first_edge = eptr1;
			eptr1->next = NIL_E;
		}
		else
		{
			eptr1->next = nptr1->first_edge;
			nptr1->first_edge = eptr1;
		}

		if (nptr2->first_edge == NULL)
		{
			nptr2->first_edge = eptr2;
			eptr2->next = NIL_E;
		}
		else
		{
			eptr2->next = nptr2->first_edge;
			nptr2->first_edge = eptr2;
		}

		++eptr1;
		++eptr2;
	}

	fclose(f);

	/*Initialisation des autres champs dans la table des sommets*/
	for (i = 0L; i < n - 1; i++)
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

	return kECSP_True;
}

Bool InitGraph_from_list(simpleEdge *Te, long n, long m, Graph *gr)
{
	int i, j;
	Edge *eptr1 = NULL, *eptr2 = NULL;
	int node1, node2;
	double cap;
	Node *nptr1 = NULL, *nptr2 = NULL;

	if (!AllocateGraph(n, m, gr))
	{
		return kECSP_False;
	}

	/*Reading of graph edges*/
	eptr1 = &(gr->Edges[0L]);
	eptr2 = &(gr->Edges[m]);

	for (j = 0; j < m; j++)
	{
		node1 = Te[j].node1;
		node2 = Te[j].node2;
		cap = Te[j].cap;

		node1--;
		node2--;
		nptr1 = &(gr->Nodes[node1]);
		nptr2 = &(gr->Nodes[node2]);
		eptr1->adjac = nptr2;
		eptr2->adjac = nptr1;
		eptr1->cap = cap;
		eptr1->X = cap;
		/*eptr1->X = 0.0;*/
		eptr2->cap = cap;
		eptr2->X = cap;
		/*eptr2->X = 0.0;*/
		eptr1->back = eptr2;
		eptr2->back = eptr1;

		eptr1->num = j + 1;
		eptr2->num = j + 1;

		if (nptr1->first_edge == NULL)
		{
			nptr1->first_edge = eptr1;
			eptr1->next = NULL;
		}
		else
		{
			eptr1->next = nptr1->first_edge;
			nptr1->first_edge = eptr1;
		}

		if (nptr2->first_edge == NULL)
		{
			nptr2->first_edge = eptr2;
			eptr2->next = NULL;
		}
		else
		{
			eptr2->next = nptr2->first_edge;
			nptr2->first_edge = eptr2;
		}

		++eptr1;
		++eptr2;
	}

	/*Initialisation des autres champs dans la table des sommets*/
	for (i = 0L; i < n - 1; i++)
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

	return kECSP_True;
}

Bool AllocateGraph(long n, long m, Graph *g)
{
	long i;

	g->n_Nodes = n;
	g->m_Edges = m;

	/*Allocation of nodes table*/
	// g->Nodes = (Node *)malloc(n*sizeof(Node));
	g->Nodes = new Node[n];

	if (g->Nodes == NULL)
	{
		printf("Unable to allocate nodes memory");
		return kECSP_False;
	}

	/*Allocation of edges table*/
	// g->Edges = (Edge *)malloc(2L*m*sizeof(Edge));
	g->Edges = new Edge[2 * m];
	if (g->Edges == NULL)
	{
		printf("Unable to allocate edges memory");

		// free(g->Nodes);
		delete[] g->Nodes;
		return kECSP_False;
	}

	for (i = 0; i < n; i++)
	{
		(g->Nodes[i]).id = i + 1;
		(g->Nodes[i]).first_edge = NIL_E;
	}

	return kECSP_True;
}

Bool DeleteGraph(Graph *g)
{
	/*free(g->Nodes);
	free(g->Edges);*/
	delete[] g->Nodes;
	delete[] g->Edges;
	/*free(g);*/

	return kECSP_True;
}

/*Permet de sauvegarder la contraction courante*/
/*dans les champs annexes*/
void save_contraction_info(Graph *gr, int niveau)
{
	long n, i;

	n = gr->n_Nodes;

	if (niveau == 1)
	{
		for (i = 0; i < n; i++)
		{
			gr->Nodes[i].a_next_gr = gr->Nodes[i].next_gr;
			gr->Nodes[i].a_n_sh = gr->Nodes[i].n_sh;
			gr->Nodes[i].a_next_sh = gr->Nodes[i].next_sh;
		}
	}
	else if (niveau == 2)
	{
		for (i = 0; i < n; i++)
		{
			gr->Nodes[i].b_next_gr = gr->Nodes[i].next_gr;
			gr->Nodes[i].b_n_sh = gr->Nodes[i].n_sh;
			gr->Nodes[i].b_next_sh = gr->Nodes[i].next_sh;
		}
	}
}

/*Contract the set of nodes given by f_w*/
Bool contract_set_w(Graph *gr, b_Node *f_w)
{
	b_Node *b_cour;
	Node *n_cour, *n_first, *n_aux;
	long k, i;

	/*We don't contract a set contains less than 2 nodes*/
	if ((f_w == NULL) || ((f_w != NULL) && (f_w->next == NULL)))
		return kECSP_False;

	/*We start the contraction by creating*/
	/*the list of the node with the field next_sh*/
	/*and by updating the representant of each node*/
	/*of the set. This representant will be the first*/
	/*node of the set.*/

	n_first = &(gr->Nodes[f_w->id - 1]);
	b_cour = f_w->next;
	while (b_cour != NIL_BN)
	{
		k = b_cour->id - 1;
		n_cour = &(gr->Nodes[k]);

		/*Update of the representant of the nodes*/
		/*that have been already contracted in the node n_cour.*/
		/*There is no matter if the node n_cour hasn't been contracted yet.*/
		n_aux = n_cour;
		do
		{
			gr->Nodes[n_aux->id - 1].n_sh = n_first;

			/*n_aux->n_sh = n_first;*/

			n_aux = n_aux->next_sh;
		} while (n_aux->id != n_cour->id);

		/*Merging the next_sh list of n_cour to the list of n_first*/
		cyclic_set_merge(n_first, n_cour);

		/*Update of the shrunk graph list*/
		/*We search the previous node of b_cour in the shunk list*/
		i = 0; /*We suppose that the node 1 is the first of the shrunk graph*/
		while (gr->Nodes[i].next_gr->id != b_cour->id)
		{
			i = gr->Nodes[i].next_gr->id - 1;
		}

		gr->Nodes[i].next_gr = gr->Nodes[k].next_gr;

		/*************/
		b_cour = b_cour->next;
	}

	return kECSP_True;
}


/*Returns a list of subsets of V that should*/
/*be contracted by the operation Theta_4.*/
/*Notice that we work on the shrunk graph.*/
/*Tab is a list of the subset of V that*/
/*should be contracted by theta_4.*/
void operation_theta_4(Graph *gr, b_Node *w, int k, b_Node **Tab, long *n_subset)
{
	long n, m;
	simpleEdge *Te;
	Mincut *gmcu;
	b_Edge *b_cour;
	b_Node *w1, *w2;
	Bool finished, w1_contracted, w2_contracted;
	Tree *ghct;
	float k_sur_2;
	ldiv_t k2;

	/*puts("Operation Theta_4");*/

	if (check_w_card(w, 4)) /*We check if |W|>= 4*/
	{
		Te = (simpleEdge *)malloc(gr->m_Edges * sizeof(simpleEdge));

		/*k_sur_2 = ceil( ((double)k/2.0) );*/
		k2 = ldiv(k, 2);
		k_sur_2 = k2.quot + k2.rem;

		n = gr->n_Nodes;

		/*CAV_VALUE_1: We search (k+1)-Edges cut*/
		m = 0;
		get_edge_list(gr, w, Te, &m, NO_0_EDGES, CAP_VALUE_1);

		/*puts("Set W");
		Print_b_Node_Set(w,stdout);*/
		/*printf("M == %ld\n",m);
		Print_eList(Te,0,m);*/

		/*Calculation of the Gomory-Hu Cut Tree of G[W] to find the (k+1)-Edges cuts*/
		ghct = (Tree *)GHCutTree(Te, n, m, (double)(k + 1 - EPSILON), (double)(k + 1 + EPSILON));

		if(ghct == NULL)
			return;
		/*puts("Tree");
		PrintGHCutTree_Table(ghct);
		Print_b_Edge_Set(ghct->b_List,stdout);*/

		/*We should search a (k+1)-edges cut with |W1|>= 2 and |W2|>= 2*/
		b_cour = ghct->b_List;
		w1 = NIL_BN;
		w2 = NIL_BN;
		finished = kECSP_False;
		while ((b_cour != NIL_BE) && (!finished))
		{
			w1 = get_cut_set(ghct, b_cour);
			w2 = get_complement_cut_set(ghct, b_cour);

			/*puts("Set W1");
			Print_b_Node_Set(w1,stdout);
			puts("Set W2");
			Print_b_Node_Set(w2,stdout);*/

			if ((check_w_card(w1, 2)) && (check_w_card(w2, 2)))
			{
				w1_contracted = kECSP_False;
				w2_contracted = kECSP_False;

				/*We check if W1 can be contracted*/
				get_edge_list(gr, w1, Te, &m, NO_0_EDGES, CAP_VALUE_X);

				if (frac_edge_index(Te, m) == -1L)
				{
					gmcu = (Mincut *)global_mincut_ho(Te, n, m);
					if (gmcu->mincap >= k_sur_2)
					{
						/*W1 will be contracted*/
						Tab[(*n_subset)] = w1;
						(*n_subset)++;

						finished = kECSP_True;
						w1_contracted = kECSP_True;

						Delete_b_Node_Set(&w2);
					}
					else
						w1_contracted = kECSP_False;

					Delete_Mincut(gmcu);
					/*free(Te);*/
				}
				else
					w1_contracted = kECSP_False;

				if (!w1_contracted)
					Delete_b_Node_Set(&w1);

				/*free(Te);*/

				if (!w1_contracted)
				{
					/*We check if W2 can be contracted*/
					get_edge_list(gr, w2, Te, &m, NO_0_EDGES, CAP_VALUE_X);

					if (frac_edge_index(Te, m) == -1L)
					{
						gmcu = (Mincut *)global_mincut_ho(Te, n, m);

						if (gmcu->mincap >= k_sur_2)
						{
							/*W2 will be contracted*/
							Tab[(*n_subset)] = w2;
							(*n_subset)++;

							finished = kECSP_True;
							w2_contracted = kECSP_True;
						}
						else
							w2_contracted = kECSP_False;

						Delete_Mincut(gmcu);
						/*free(Te);*/
					}
					else
						w2_contracted = kECSP_False;

					if (!w2_contracted)
						Delete_b_Node_Set(&w2);
				}

				if ((!w1_contracted) && (!w2_contracted))
				{
					finished = kECSP_False;

					/*Delete_b_Node_Set(&w1);
					Delete_b_Node_Set(&w2);*/
				}
			}
			else
			{
				Delete_b_Node_Set(&w1);
				Delete_b_Node_Set(&w2);
			}

			b_cour = b_cour->next;
		}

		Delete_Tree(&ghct);
		free(Te);
	}
}

/*REcherche les aretes multiples  e 1 et les contracte*/
void operation_theta_4_4(Graph *gr, b_Node *w, int k, b_Node **Tab, long *n_subset)
{
	long i;
	long u_cour = 0, u = 0, v = 0;
	b_Edge *delta_v, *delta_w, *be_cour;
	Node *nptr, *new_node;
	long n, m, nb_voisin, w_card;
	Bool nouveau_sommet = kECSP_True;
	short *marquage;
	Bool Trouver;

	n = gr->n_Nodes;
	m = gr->m_Edges;

	marquage = new short[n];

	for (i = 0; i < n; i++)
	{
		mark_tab[i] = 0;
	}

	Trouver = kECSP_False;
	nptr = &(gr->Nodes[0]);
	do
	{
		mark_tab[nptr->id - 1] = 1;

		for (i = 0; i < n; i++)
		{
			marquage[i] = 0;
		}

		delta_v = get_delta_v(gr, nptr->id, NO_0_EDGES, CAP_VALUE_X);
		be_cour = delta_v;
		while (be_cour != NIL_BE)
		{
			if (be_cour->cap == 1.0)
			{
				marquage[be_cour->node2 - 1]++;
			}

			be_cour = be_cour->next;
		}

		/*Comptage du nombre de voisin par arete double*/
		nb_voisin = 0;
		for (i = 0; i < n; i++)
		{
			if (marquage[i] >= 2)
			{
				nb_voisin++;
				u_cour = i + 1;
			}
		}

		if (nb_voisin >= 1)
		{
			/*On verifie qu'il n'y a pas d'arete*/
			/*fractionnaire entre les deux*/
			be_cour = delta_v;
			while ((be_cour != NIL_BE) && ((be_cour->cap == 1.0) || ((be_cour->cap < 1.0) && (be_cour->node2 !=
																							  u_cour))))
			{
				be_cour = be_cour->next;
			}

			/*Dans ce cas il n'y a pas d'arete factionnaire*/
			if (be_cour == NIL_BE)
			{
				w = (b_Node *)malloc(sizeof(b_Node));
				w->id = nptr->id;

				w->next = (b_Node *)malloc(sizeof(b_Node));
				w->next->id = u_cour;

				w->next->next = NIL_BN;

				// mark_tab[nptr->id-1] = 1;
				// mark_tab[u_cour-1] = 1;

				/*On verifie en la cardinalite*/
				delta_w = get_delta_w(gr, w, NO_0_EDGES, CAP_VALUE_X);
				w_card = 0;
				be_cour = delta_w;
				while (be_cour != NIL_BE)
				{
					w_card++;
					be_cour = be_cour->next;
				}

				Delete_b_Edge_Set(&delta_w);

				if ((k <= w_card) && (w_card <= k + 1))
				{
					Tab[*n_subset] = w;
					(*n_subset)++;

					Trouver = kECSP_True;
				}
				else
				{
					Delete_b_Node_Set(&w);
					Trouver = kECSP_False;
				}
			}
			else
			{
				Trouver = kECSP_False;
			}
		}

		Delete_b_Edge_Set(&delta_v);

		if (Trouver == kECSP_False)
		{
			/*Recherche d'un nouveau sommet*/
			new_node = nptr;
			do
			{
				new_node = new_node->next_gr;
			} while ((new_node != &(gr->Nodes[0])) && (mark_tab[new_node->id - 1] != 0));

			/*nptr = nptr->next_gr;*/
			if (new_node != &(gr->Nodes[0]))
			{
				nouveau_sommet = kECSP_True;

				nptr = new_node;
			}
			else
			{
				nouveau_sommet = kECSP_False;
			}
		}
		else
		{
			nouveau_sommet = kECSP_False;
		}

	} while (nouveau_sommet == kECSP_True);

	delete[] marquage;

	/*puts("CH == ");
	for(i=0;i<(*nb_chaine);i++)
	{
		Print_b_Node_Set(chaine[i],stdout);
	}*/
}

/*This function contracts the ((k/2)+1)-clique whose degree is equal to k+1*/
/*This function should work only for complete graphs such as TSPLIB graphs.*/
void operation_theta_4_bis(Graph *gr, b_Node *w, int k, b_Node **Tab, long *n_subset)
{
	long n, m;
	long k_sur_2;
	b_Node **clique;
	long nb_clique;
	long i;

	/*n = gr->n_Nodes;
	m = gr->m_Edges;*/

	n = 100;
	m = gr->m_Edges;

	/*Calcul de k/2*/
	k_sur_2 = sup_part(k, 2);

	clique = (b_Node **)malloc(n * sizeof(b_Node *));
	nb_clique = 0;

	*n_subset = 0;

	for (i = 0; i < n; i++)
	{
		clique[i] = NIL_BN;
	}

	/*Calcul des cliques e (k/2)+1 sommets*/
	get_clique(gr, k_sur_2 + 1, k, k + 1, clique, &nb_clique);

	if (nb_clique > 0)
	{
		/*Comptage des degres de |delta(W)| pour chaque clique*/
		/*Si le degre est egal e k+1 alors on contracte la clique*/
		for (i = 0; i < nb_clique; i++)
		{
			Tab[*n_subset] = clique[i];
			(*n_subset)++;
		}
	}
	else
	{
		for (i = 0; i < n; i++)
			Delete_b_Node_Set(&clique[i]);
	}

	free(clique);
}

/*Returns a list of subsets of V that should*/
/*be contracted by the operation Theta_3.*/
/*Notice that we work on the shrunk graph.*/
/*Tab is a list of the subset of V that*/
/*should be contracted by theta_3.*/
void operation_theta_3_bis(Graph *gr, b_Node *w, int k, b_Node **Tab, long *n_subset)
{
	long n, m;
	simpleEdge *Te;
	b_Edge *b_cour;
	b_Node *w1, *w2;
	Bool finished;
	Tree *ghct;

	/*We check if |W|>=4*/
	if (check_w_card(w, 4))
	{
		Te = (simpleEdge *)malloc(gr->m_Edges * sizeof(simpleEdge));

		/*puts("Set W");
		Print_b_Node_Set(w,stdout);*/

		n = gr->n_Nodes;
		m = 0;
		/*Search of a k-edges cut*/
		/*Relabelling of the edges of G[W]*/
		get_edge_list(gr, w, Te, &m, NO_0_EDGES, CAP_VALUE_1); /*CAV_VALUE_1: We search k-Edges cut*/

		/*printf("M == %ld\n",m);
		Print_eList(Te,0,m);*/

		/*Calculation of the Gomory-Hu Cut Tree of G[W] to find the (k)-Edges cuts*/
		ghct = (Tree *)GHCutTree(Te, n, m, (double)(k - EPSILON), (double)(k + EPSILON));

		if(ghct == NULL)
			return;
		/*puts("GH Cut Tree");
		PrintGHCutTree_Table(ghct);
		puts("B_List");
		Print_b_Edge_Set(ghct->b_List,stdout);*/

		/*free(Te);*/

		b_cour = ghct->b_List;
		finished = kECSP_False;
		while ((b_cour != NIL_BE) && (!finished))
		{
			w1 = get_cut_set(ghct, b_cour);
			w2 = get_complement_cut_set(ghct, b_cour);

			/*We check if |W1|>=2 and |W2|>=2*/
			/*and then we check if G[W1] contains*/
			/*a fractionnary edge. In this case we*/
			/*contract the set W2. Else we check for W2.*/
			if ((check_w_card(w1, 2)) && (check_w_card(w2, 2)))
			{
				m = 0;
				get_edge_list(gr, w1, Te, &m, NO_0_EDGES, CAP_VALUE_X);

				/*We check if G[W1] doesn't contains any fractionnary edge*/
				/*In this case we contract W1*/
				if (frac_edge_index(Te, m) == -1L)
				{
					Tab[*n_subset] = w1;
					(*n_subset)++;

					finished = kECSP_True;

					Delete_b_Node_Set(&w2);
				}
				else
				{
					m = 0;
					get_edge_list(gr, w2, Te, &m, NO_0_EDGES, CAP_VALUE_X);

					if (frac_edge_index(Te, m) == -1L)
					{
						Tab[*n_subset] = w2;
						(*n_subset)++;

						finished = kECSP_True;

						Delete_b_Node_Set(&w1);
					}
					else
					{
						Delete_b_Node_Set(&w1);
						Delete_b_Node_Set(&w2);
					}
				}
			}
			else
			{
				Delete_b_Node_Set(&w1);
				Delete_b_Node_Set(&w2);
			}

			b_cour = b_cour->next;
		}

		Delete_Tree(&ghct);
		free(Te);
	}
}


/*Returns the list of the edges of the subgraph G[W]*/
/*n_tab will be the length of the returned table.*/
/*edge_flag says if you don't want the edges that have X(e)=0 or*/
/*if you want all the edges in G[W].*/
/*cap_flag says what kind of value you want as capacity.*/
/*See graph.h to have the different values of edge_flag or cap_flag.*/
void get_edge_list(Graph *g, b_Node *w, simpleEdge *edge_tab, long *n_tab, int edge_flag, int cap_flag)
{
	b_Node *b_cour;
	long k, nod1, nod2;
	double cap;
	Edge *eptr;
	long n_u, n_v, i;

	if (w == NIL_BN)
	{
		fprintf(stdout, "%s", "ERROR.Get_edge_List. Set W is empty.");
		/*return NULL;*/
		Print_b_Node_Set(w, stdout);
		return;
	}

	/*Marking up all the nodes of the set W*/
	if (g->n_Nodes > N_TAB)
	{
		puts("N nodes superieur e N_TAB");

		exit(-1);
	}

	/*if(mark_tab == NULL)
	{
		fprintf(stdout,"%s","ERROR.Get_edge_List. Unable to allocate memory for node marking.");
		return NULL;
	}*/

	for (i = 0; i < g->n_Nodes; i++)
	{
		mark_tab[i] = NODE_NOT_MARKED;
	}

	b_cour = w;
	while (b_cour != NIL_BN)
	{
		/*printf("%ld  ",b_cour->id);*/

		i = b_cour->id - 1;

		mark_tab[i] = NODE_MARKED;

		b_cour = b_cour->next;
	}

	/*Construction of the list of edges*/
	if (g->m_Edges > M_TAB)
	{
		puts("M edges superieur e M_TAB");

		exit(-1);
	}

	/*We mark 0 on a node to say that this node is not*/
	/*in W. If we mark a non zero value on the node(1 or 2),*/
	/*it means that the node is in W. To build the list E(W),*/
	/*we follow the incidence list of each node of W.*/
	/*When we finish traversing the incidence list of*/
	/*a node u, we mark 2 in mark_tab to say that we*/
	/*have already traverse the incidence list of u.*/
	/*Then when we are traversing the incidence list*/
	/*of a node v(let vk be an edge of this list),we*/
	/*check if k is marked with 2 or not.If k is marked,*/
	/*it means that we have already traverse the incidence*/
	/*list of k, and we have already taken the edge kv.*/
	/*Then we don't take the edge vk.*/
	/*If the node k is not marked with 2, then we can*/
	/*take the edge vk.*/

	/*puts("xxxx");
	for(i=0;i<2*g->m_Edges;i++)
	{
		printf("e == %x    %ld  %ld  %lx  %lx  %f\n",&(g->Edges[i]),g->Edges[i].adjac->id,g->Edges[i].back->adjac->id,g->Edges[i].next,g->Edges[i].back,g->Edges[i].cap);
		//printf("e == %ld %ld    %x\n",Gr.Edges[i].adjac->id,Gr.Edges[i].back->adjac->id,Gr.Edges[i].next);

	}*/

	k = 0;

	switch (edge_flag)
	{
	case NO_0_EDGES: /*Suppression des arètes nulles: opération theta_0*/
		b_cour = w;

		while (b_cour != NIL_BN)
		{
			n_u = b_cour->id;
			do
			{
				/*printf("u == %ld  u' == %ld\n",b_cour->id,n_u);*/

				eptr = g->Nodes[n_u - 1].first_edge;

				while (eptr != NIL_E)
				{
					/*puts("xxxx");

					printf("mark_Tab == %ld  mark_Tab+m == %ld  Nodes == %ld  Nodes+n == %ld\n",mark_tab,mark_tab+g->m_Edges-1,g->Nodes,g->Nodes+g->n_Nodes-1);

					for(i=0;i<g->n_Nodes;i++)
					{
						printf("Nodes[%ld] == %ld \n",i+1,g->Nodes[i].id);
					}*/

					/*for(i=0;i<2*g->m_Edges;i++)
					{
						printf("e == %x    %ld  %ld  %lx  %lx  %f\n",&(g->Edges[i]),g->Edges[i].adjac->id,g->Edges[i].back->adjac->id,g->Edges[i].next,g->Edges[i].back,g->Edges[i].cap);
						//printf("e == %ld %ld    %x\n",Gr.Edges[i].adjac->id,Gr.Edges[i].back->adjac->id,Gr.Edges[i].next);
					}*/

					/*printf("eptr->id == %ld
%x\n",eptr->adjac->id,eptr);*/

					n_v = eptr->adjac->id;
					/*printf("u == %ld  u' == %ld  v == %ld\n",b_cour->id,n_u,n_v);*/

					if ((g->Nodes[n_v - 1].n_sh->id != g->Nodes[n_u - 1].n_sh->id) /*Si u et v ne sont pas contractes ensembles*/
						&&
						(mark_tab[g->Nodes[n_v - 1].n_sh->id - 1] != NODE_NOT_MARKED) /*Si le sommet representant de v est marque*/
						&&
						(mark_tab[g->Nodes[n_v - 1].n_sh->id - 1] != NODE_LIST_TERMINATED) /*Si le sommet v n'a pas deje ete parcouru*/
						&&
						(eptr->X != 0.0)) /*Si l'arete est non nulle*/

					{
						/*printf("u == %ld  u' == %ld  v == %ld\n",b_cour->id,n_u,n_v);*/

						nod1 = g->Nodes[n_u - 1].n_sh->id;
						nod2 = g->Nodes[n_v - 1].n_sh->id;

						switch (cap_flag)
						{
						case CAP_VALUE_CAP:
							cap = eptr->cap;
							break;

						case CAP_VALUE_X:
							cap = eptr->X;
							break;

						case CAP_VALUE_1:
							cap = 1.0;
							break;

						case CAP_VALUE_0:
							cap = 0.0;
							break;

						case CAP_VALUE_0_FOR_FRAC_EDGE:
							if ((1.0 - eptr->X) > 0.0) /*Si arète fractionnaire*/
							{
								cap = 0.0;
							}
							else
							{
								cap = eptr->X;
							}
							break;

						case CAP_VALUE_M_FOR_FRAC_EDGE:
							if ((1.0 - eptr->X) > 0.0) /*Si arète fractionnaire*/
							{
								cap = (double)g->m_Edges;
							}
							else
							{
								cap = eptr->X;
							}
							break;

						default: /*Par defaut on prendra la valeur X*/
							puts("Erreur dans les FLAGS de GET_EDGE_LIST");
							exit(-1);
							cap = eptr->X;
							break;
						}

						edge_tab[k].node1 = nod1;
						edge_tab[k].node2 = nod2;
						edge_tab[k].cap = cap;
						k++;

						/*printf("K == %ld  Tab.node1 == %ld Tab.node2 == %ld Node1 == %ld  Node2 ==
						%ld\n",k,edge_tab[k].node1,edge_tab[k].node2,nod1,nod2);*/
					}

					eptr = eptr->next;
				}

				n_u = g->Nodes[n_u - 1].next_sh->id;

			} while (n_u != b_cour->id);

			/*When we finished the contraction list of b_cour*/
			/*we must set that nodes of the contraction list*/
			/*have all been traverse. So we shouldn't take their*/
			/*incident edges again.*/
			mark_tab[b_cour->id - 1] = NODE_LIST_TERMINATED;

			b_cour = b_cour->next;
		}

		break;

	case ALL_EDGES:
		b_cour = w;

		while (b_cour != NIL_BN)
		{
			n_u = b_cour->id;

			do
			{
				eptr = g->Nodes[n_u - 1].first_edge;

				while (eptr != NIL_E)
				{
					n_v = eptr->adjac->id;

					if ((g->Nodes[n_v - 1].n_sh->id != g->Nodes[n_u - 1].n_sh->id) /*Si u et v ne sont pas contractes ensembles*/
						&&
						(mark_tab[g->Nodes[n_v - 1].n_sh->id - 1] != NODE_NOT_MARKED) /*Si le sommet representant de v est marque*/
						&&
						(mark_tab[g->Nodes[n_v - 1].n_sh->id - 1] != NODE_LIST_TERMINATED)) /*Si le sommet v n'a pas deje ete parcouru*/
					{
						nod1 = g->Nodes[n_u - 1].n_sh->id;
						nod2 = g->Nodes[n_v - 1].n_sh->id;

						switch (cap_flag)
						{
						case CAP_VALUE_CAP:
							cap = eptr->cap;
							break;

						case CAP_VALUE_X:
							cap = eptr->X;
							break;

						case CAP_VALUE_1:
							cap = 1.0;
							break;

						case CAP_VALUE_0:
							cap = 0.0;
							break;

						case CAP_VALUE_0_FOR_FRAC_EDGE:
							if ((1.0 - eptr->X) > 0) /*Si arète fractionnaire*/
							{
								cap = 0.0;
							}
							else
							{
								cap = eptr->X;
							}
							break;

						case CAP_VALUE_M_FOR_FRAC_EDGE:
							if ((1.0 - eptr->X) > 0.0) /*Si arète fractionnaire*/
							{
								cap = (double)g->m_Edges;
							}
							else
							{
								cap = eptr->X;
							}
							break;

						default: /*Par defaut on prendra la valeur X*/
							puts("Erreur dans les FLAGS de GET_EDGE_LIST");
							exit(-1);
							cap = eptr->X;
							break;
						}

						edge_tab[k].node1 = nod1;
						edge_tab[k].node2 = nod2;
						edge_tab[k].cap = cap;

						k++;
					}

					eptr = eptr->next;
				}

				n_u = g->Nodes[n_u - 1].next_sh->id;

			} while (n_u != b_cour->id);

			/*When we finished the contraction list of b_cour*/
			/*we must set that nodes of the contraction list*/
			/*have all been traverse. So we shouldn't take their*/
			/*incident edges again.*/
			mark_tab[b_cour->id - 1] = NODE_LIST_TERMINATED;

			b_cour = b_cour->next;
		}

		break;
	}

	(*n_tab) = k;
}

/*Return a subset of nodes which constitutes a violated cut of the graph.*/
/*ghct represents the the Gomory-Hu cut tree given by the GHCutTree_Tree function.*/
/*b_cour represents an edge of the ghct that is violated*/
b_Node *get_cut_set(Tree *ghct, b_Edge *b_cour)
{
	Bool found = kECSP_False;
	b_Node *f_w;

	f_w = NIL_BN;

	create_cut_set(&(ghct->Tab[0]), &f_w, b_cour, &found);

	return f_w;
}

/*Return the set W_b of b_Node that should represent the complementary*/
/*of the set obtained with get_cut_set().*/
b_Node *get_complement_cut_set(Tree *ghct, b_Edge *b_cour)
{
	b_Node *f_w;
	long i;

	/*Search for node2 in ghct table*/
	i = 0;
	while ((i < ghct->n_nodes) && (ghct->Tab[i].id != b_cour->node2))
		i++;

	if (i >= ghct->n_nodes)
	{
		fprintf(stdout, "Erreur dans Get_Complement_Cut. b_cour == %ld %ld  %f\n", b_cour->node1, b_cour->node2, b_cour->cap);
		return NULL;
	}

	f_w = NIL_BN;

	create_complement_cut_set(&(ghct->Tab[i]), &f_w, b_cour);

	return f_w;
}

/*Checks if the cardinal of the set W is greater or equal to a.*/
/*You can use this function to check if |W| < b by checking if*/
/*the result of check_w_card(W,b) is false.*/
Bool check_w_card(b_Node *w, long a)
{
	b_Node *b_cour;
	long card;

	card = 0;
	b_cour = w;
	while ((b_cour != NIL_BN) && (card < a))
	{
		card++;
		b_cour = b_cour->next;
	}

	if (card >= a)
		return kECSP_True;
	else
		return kECSP_False;
}

/*Says if the set of edges contains a fractionnary edge*/
/*It returns the index of the first fractionnary edge.*/
/*It will return -1 if there isn't any frac edge in the table.*/
long frac_edge_index(simpleEdge *Tab, long m_tab)
{
	long i;
	Bool found;

	i = m_tab - 1;
	found = kECSP_False;
	while ((i >= 0) && (found == kECSP_False))
	{
		if ((Tab[i].cap == 0.0) || (Tab[i].cap == 1.0))
		{
			found = kECSP_False;

			i--;
		}
		else
		{
			found = kECSP_True;
		}
	}

	return i;
}

/*Makes the fusion of two cyclic linked lists*/
/*The lists must have at least 1 node.*/
Bool cyclic_set_merge(Node *f1, Node *f2)
{
	Node *aux;

	if ((f1 == NIL_N) || (f2 == NIL_N))
		return kECSP_False;

	aux = f2->next_sh;
	f2->next_sh = f1->next_sh;
	f1->next_sh = aux;

	return kECSP_True;
}

/*Gives the list of the edges of delta(v).*/
/*The method is similar to the method used in get_delta_w().*/
b_Edge *get_delta_v(Graph *g, long v, int edge_flag, int cap_flag)
{
	long k, nod1, nod2;
	double cap;
	Edge *eptr;
	long n_u, n_v, num;

	b_Edge *delta_v, *edge_cour;

	if (v > g->n_Nodes)
	{
		fprintf(stdout, "%s", "ERROR.Get_delta_v. Node v is out of bound.");
		return NULL;
	}

	/*Marking up all the nodes of the set W*/
	if (g->n_Nodes > N_TAB)
	{
		puts("N nodes superieur e N_TAB");

		exit(-1);
	}

	/*Construction of the list of edges*/
	if (g->m_Edges > M_TAB)
	{
		puts("M edges superieur e M_TAB");

		exit(-1);
	}

	/*puts("xxxx");
	for(i=0;i<2*g->m_Edges;i++)
	{
		printf("e == %x    %ld  %ld  %lx  %lx  %f\n",&(g->Edges[i]),g->Edges[i].adjac->id,g->Edges[i].back->adjac->id,g->Edges[i].next,g->Edges[i].back,g->Edges[i].cap);
		//printf("e == %ld %ld    %x\n",Gr.Edges[i].adjac->id,Gr.Edges[i].back->adjac->id,Gr.Edges[i].next);

	}*/

	k = 0;

	delta_v = NIL_BE;

	switch (edge_flag)
	{
	case NO_0_EDGES: /*Suppression des arètes nulles: opération theta_0*/

		n_u = v;
		do
		{
			eptr = g->Nodes[n_u - 1].first_edge;

			while (eptr != NIL_E)
			{
				n_v = eptr->adjac->id;

				if ((g->Nodes[n_v - 1].n_sh->id != g->Nodes[n_u - 1].n_sh->id) /*Si u et v ne sont pas contractes ensembles*/
					&&
					(eptr->X != 0.0)) /*Si l'arete est non nulle*/
				{
					nod1 = g->Nodes[n_u - 1].n_sh->id;
					nod2 = g->Nodes[n_v - 1].n_sh->id;
					num = eptr->num;

					switch (cap_flag)
					{
					case CAP_VALUE_CAP:
						cap = eptr->cap;
						break;

					case CAP_VALUE_X:
						cap = eptr->X;
						break;

					case CAP_VALUE_1:
						cap = 1L;
						break;

					case CAP_VALUE_0:
						cap = 0L;
						break;

					case CAP_VALUE_0_FOR_FRAC_EDGE:
						if ((1L - eptr->X) > 0) /*Si arète fractionnaire*/
						{
							cap = 0L;
						}
						else
						{
							cap = eptr->X;
						}
						break;

					case CAP_VALUE_M_FOR_FRAC_EDGE:
						if ((1L - eptr->X) > 0) /*Si arète fractionnaire*/
						{
							cap = (double)g->m_Edges;
						}
						else
						{
							cap = eptr->X;
						}
						break;

					default: /*Par defaut on prendra la valeur X*/
						puts("Erreur dans les FLAGS de GET_DELTA_V");
						exit(-1);
						cap = eptr->X;
						break;
					}

					edge_cour = (b_Edge *)malloc(sizeof(b_Edge));

					edge_cour->node1 = nod1;
					edge_cour->node2 = nod2;
					edge_cour->cap = cap;
					edge_cour->num = num;

					edge_cour->next = delta_v;
					delta_v = edge_cour;

					/*edge_tab[k].node1 = nod1;
					edge_tab[k].node2 = nod2;
					edge_tab[k].cap = cap;

					k++;*/
				}

				eptr = eptr->next;
			}

			n_u = g->Nodes[n_u - 1].next_sh->id;

		} while (n_u != v);

		break;

	case ALL_EDGES:
		n_u = v;
		do
		{
			eptr = g->Nodes[n_u - 1].first_edge;

			while (eptr != NIL_E)
			{
				n_v = eptr->adjac->id;

				if ((g->Nodes[n_v - 1].n_sh->id != g->Nodes[n_u - 1].n_sh->id)) /*Si u et v ne sont pas contractes ensembles*/
				{
					nod1 = g->Nodes[n_u - 1].n_sh->id;
					nod2 = g->Nodes[n_v - 1].n_sh->id;
					num = eptr->num;

					switch (cap_flag)
					{
					case CAP_VALUE_CAP:
						cap = eptr->cap;
						break;

					case CAP_VALUE_X:
						cap = eptr->X;
						break;

					case CAP_VALUE_1:
						cap = 1.0;
						break;

					case CAP_VALUE_0:
						cap = 0.0;
						break;

					case CAP_VALUE_0_FOR_FRAC_EDGE:
						if ((1.0 - eptr->X) > 0) /*Si arète fractionnaire*/
						{
							cap = 0.0;
						}
						else
						{
							cap = eptr->X;
						}
						break;

					case CAP_VALUE_M_FOR_FRAC_EDGE:
						if ((1.0 - eptr->X) > 0) /*Si arète fractionnaire*/
						{
							cap = (double)g->m_Edges;
						}
						else
						{
							cap = eptr->X;
						}
						break;

					default: /*Par defaut on prendra la valeur X*/
						puts("Erreur dans les FLAGS de GET_DELTA_V");
						exit(-1);
						cap = eptr->X;
						break;
					}

					edge_cour = (b_Edge *)malloc(sizeof(b_Edge));

					edge_cour->node1 = nod1;
					edge_cour->node2 = nod2;
					edge_cour->cap = cap;
					edge_cour->num = num;

					edge_cour->next = delta_v;
					delta_v = edge_cour;

					/*edge_tab[k].node1 = nod1;
					edge_tab[k].node2 = nod2;
					edge_tab[k].cap = cap;

					k++;*/
				}

				eptr = eptr->next;
			}

			n_u = g->Nodes[n_u - 1].next_sh->id;

		} while (n_u != v);

		break;
	}

	return delta_v;
}

/*Gives the list of the edges of delta(W).*/
/*The method is similar to the method used in get_edge_list().*/
b_Edge *get_delta_w(Graph *g, b_Node *w, int edge_flag, int cap_flag)
{
	b_Node *b_cour;
	long k, nod1, nod2;
	double cap;
	Edge *eptr;
	long n_u, n_v, i;
	long num;

	b_Edge *delta_w, *edge_cour;

	if (w == NIL_BN)
	{
		return NIL_BE;
	}

	/*Marking up all the nodes of the set W*/
	if (g->n_Nodes > N_TAB)
	{
		puts("N nodes superieur e N_TAB");

		exit(-1);
	}

	for (i = 0; i < g->n_Nodes; i++)
	{
		mark_tab[i] = NODE_NOT_MARKED;
	}

	b_cour = w;
	while (b_cour != NIL_BN)
	{
		i = b_cour->id - 1;

		mark_tab[i] = NODE_MARKED;

		b_cour = b_cour->next;
	}

	/*Construction of the list of edges*/
	if (g->m_Edges > M_TAB)
	{
		puts("M edges superieur e M_TAB");

		exit(-1);
	}

	/*We mark 0 on a node to say that this node is not*/
	/*in W. If we mark a non zero value on the node(1 or 2),*/
	/*it means that the node is in W. To build the list E(W),*/
	/*we follow the incidence list of each node of W.*/
	/*When we finish traversing the incidence list of*/
	/*a node u, we mark 2 in mark_tab to say that we*/
	/*have already traverse the incidence list of u.*/
	/*Then when we are traversing the incidence list*/
	/*of a node v(let vk be an edge of this list),we*/
	/*check if k is marked with 2 or not.If k is marked,*/
	/*it means that we have already traverse the incidence*/
	/*list of k, and we have already taken the edge kv.*/
	/*Then we don't take the edge vk.*/
	/*If the node k is not marked with 2, then we can*/
	/*take the edge vk.*/

	/*puts("xxxx");
	for(i=0;i<2*g->m_Edges;i++)
	{
		printf("e == %x    %ld  %ld  %lx  %lx  %f\n",&(g->Edges[i]),g->Edges[i].adjac->id,g->Edges[i].back->adjac->id,g->Edges[i].next,g->Edges[i].back,g->Edges[i].cap);
		//printf("e == %ld %ld    %x\n",Gr.Edges[i].adjac->id,Gr.Edges[i].back->adjac->id,Gr.Edges[i].next);

	}*/

	k = 0;

	delta_w = NIL_BE;

	switch (edge_flag)
	{
	case NO_0_EDGES: /*Suppression des arètes nulles: opération theta_0*/
		b_cour = w;

		while (b_cour != NIL_BN)
		{
			n_u = b_cour->id;
			do
			{
				eptr = g->Nodes[n_u - 1].first_edge;

				while (eptr != NIL_E)
				{
					n_v = eptr->adjac->id;

					if ((g->Nodes[n_v - 1].n_sh->id != g->Nodes[n_u - 1].n_sh->id) /*Si u et v ne sont pas contractes ensembles*/
						&&
						(mark_tab[g->Nodes[n_v - 1].n_sh->id - 1] == NODE_NOT_MARKED) /*Si le sommet representant de v n'est pas marque*/
						&&
						(eptr->X != 0.0)) /*Si l'arete est non nulle*/

					{
						nod1 = g->Nodes[n_u - 1].n_sh->id;
						nod2 = g->Nodes[n_v - 1].n_sh->id;
						num = eptr->num;

						switch (cap_flag)
						{
						case CAP_VALUE_CAP:
							cap = eptr->cap;
							break;

						case CAP_VALUE_X:
							cap = eptr->X;
							break;

						case CAP_VALUE_1:
							cap = 1.0;
							break;

						case CAP_VALUE_0:
							cap = 0.0;
							break;

						case CAP_VALUE_0_FOR_FRAC_EDGE:
							if ((1.0 - eptr->X) > 0) /*Si arète fractionnaire*/
							{
								cap = 0.0;
							}
							else
							{
								cap = eptr->X;
							}
							break;

						case CAP_VALUE_M_FOR_FRAC_EDGE:
							if ((1.0 - eptr->X) > 0.0) /*Si arète fractionnaire*/
							{
								cap = (double)g->m_Edges;
							}
							else
							{
								cap = eptr->X;
							}
							break;

						default: /*Par defaut on prendra la valeur X*/
							puts("Erreur dans les FLAGS de GET_DELTA_W");
							exit(-1);
							cap = eptr->X;
							break;
						}

						edge_cour = (b_Edge *)malloc(sizeof(b_Edge));

						edge_cour->node1 = nod1;
						edge_cour->node2 = nod2;
						edge_cour->cap = cap;
						edge_cour->num = num;

						edge_cour->next = delta_w;
						delta_w = edge_cour;

						/*edge_tab[k].node1 = nod1;
						edge_tab[k].node2 = nod2;
						edge_tab[k].cap = cap;

						k++;*/
					}

					eptr = eptr->next;
				}

				n_u = g->Nodes[n_u - 1].next_sh->id;

			} while (n_u != b_cour->id);

			/*When we finished the contraction list of b_cour*/
			/*we must set that nodes of the contraction list*/
			/*have all been traverse. So we shouldn't take their*/
			/*incident edges again.*/
			/*mark_tab[b_cour->id - 1] = NODE_LIST_TERMINATED;*/

			b_cour = b_cour->next;
		}

		break;

	case ALL_EDGES:
		b_cour = w;

		while (b_cour != NIL_BN)
		{
			n_u = b_cour->id;

			do
			{
				eptr = g->Nodes[n_u - 1].first_edge;

				while (eptr != NIL_E)
				{
					n_v = eptr->adjac->id;

					if ((g->Nodes[n_v - 1].n_sh->id != g->Nodes[n_u - 1].n_sh->id) /*Si u et v ne sont pas contractes ensembles*/
						&&
						(mark_tab[g->Nodes[n_v - 1].n_sh->id - 1] == NODE_NOT_MARKED)) /*Si le sommet representant de v n'est pas marque*/
					{
						nod1 = g->Nodes[n_u - 1].n_sh->id;
						nod2 = g->Nodes[n_v - 1].n_sh->id;
						num = eptr->num;

						switch (cap_flag)
						{
						case CAP_VALUE_CAP:
							cap = eptr->cap;
							break;

						case CAP_VALUE_X:
							cap = eptr->X;
							break;

						case CAP_VALUE_1:
							cap = 1.0;
							break;

						case CAP_VALUE_0:
							cap = 0.0;
							break;

						case CAP_VALUE_0_FOR_FRAC_EDGE:
							if ((1.0 - eptr->X) > 0.0) /*Si arète fractionnaire*/
							{
								cap = 0.0;
							}
							else
							{
								cap = eptr->X;
							}
							break;

						case CAP_VALUE_M_FOR_FRAC_EDGE:
							if ((1.0 - eptr->X) > 0.0) /*Si arète fractionnaire*/
							{
								cap = (double)g->m_Edges;
							}
							else
							{
								cap = eptr->X;
							}
							break;

						default: /*Par defaut on prendra la valeur X*/
							puts("Erreur dans les FLAGS de GET_DELTA_W");
							exit(-1);
							cap = eptr->X;
							break;
						}

						edge_cour = (b_Edge *)malloc(sizeof(b_Edge));

						edge_cour->node1 = nod1;
						edge_cour->node2 = nod2;
						edge_cour->cap = cap;
						edge_cour->num = num;

						edge_cour->next = delta_w;
						delta_w = edge_cour;

						/*edge_tab[k].node1 = nod1;
						edge_tab[k].node2 = nod2;
						edge_tab[k].cap = cap;

						k++;*/
					}

					eptr = eptr->next;
				}

				n_u = g->Nodes[n_u - 1].next_sh->id;

			} while (n_u != b_cour->id);

			/*When we finished the contraction list of b_cour*/
			/*we must set that nodes of the contraction list*/
			/*have all been traverse. So we shouldn't take their*/
			/*incident edges again.*/
			/*mark_tab[b_cour->id - 1] = NODE_LIST_TERMINATED;*/

			b_cour = b_cour->next;
		}

		break;
	}

	return delta_w;
}

/*This function returns a set of nodes which constitute a k-clique*/
void get_clique(Graph *gr, int k, int k_card_1, long k_card_2, b_Node *clique[], long *n_clique)
{
	long i, j, n, m;
	Bool nouveau_sommet, is_voisin, clique_trouver = kECSP_False;
	long *KFile, n_file, *deconnexion;
	long nb_clique;
	Node *nptr, *Ui, *new_node;
	b_Node *K, *n_cour;
	long val_max_deconnexion, max_deconnexion;
	b_Edge *delta_v, *be_cour, *delta_w;
	long deg;
	short *marquage_som;
	long n_p;

	n = gr->n_Nodes;
	m = gr->m_Edges;

	KFile = (long *)malloc(n * sizeof(long));

	deconnexion = (long *)malloc(n * sizeof(long));

	marquage_som = (short *)malloc(n * sizeof(short));

	nb_clique = 0;

	for (i = 0; i < n; i++)
	{
		marquage_som[i] = 0;
	}

	/*On compte le nombre de sommets du graphe reduit */
	n_p = 0;
	nptr = &(gr->Nodes[0]);
	do
	{
		n_p++;
		nptr = nptr->next_gr;

	} while (nptr != &(gr->Nodes[0]));

	nouveau_sommet = kECSP_True;
	nptr = &(gr->Nodes[0]);
	do
	{
		/*for(i=0;i<k+1;i++)*/
		for (i = 0; i < n; i++)
		{
			KFile[i] = 0;
		}

		/*On enfile NPTR et tous ses voisins qui*/
		/*sont relies e lui par une arete e 1.*/
		n_file = 1;
		KFile[n_file - 1] = nptr->id;

		delta_v = get_delta_v(gr, nptr->id, NO_0_EDGES, CAP_VALUE_X);
		be_cour = delta_v;
		while (be_cour != NIL_BE)
		{
			if (be_cour->cap == 1.0)
			{
				n_file++;
				KFile[n_file - 1] = be_cour->node2;
			}

			be_cour = be_cour->next;
		}

		/*Print_b_Edge_Set(delta_v,stdout);*/

		Delete_b_Edge_Set(&delta_v);

		/*e = nptr->first_edge;
		while(e != NIL_E)
		{
			if(e->X == 1.0)
			{
				n_file++;
				KFile[n_file-1] = e->adjac->id;
			}

			e = e->next;
		}*/

		do
		{
			for (i = 0; i < n; i++)
			{
				deconnexion[i] = 0;
			}

			/*On calcul le nombre de deconnexion pour chaque sommet*/
			for (i = 0; i < n_file - 1; i++)
			{
				Ui = &(gr->Nodes[KFile[i] - 1]);
				delta_v = get_delta_v(gr, Ui->id, NO_0_EDGES, CAP_VALUE_X);

				for (j = i + 1; j < n_file; j++)
				{
					/*Recheche si Uj est voisin de Ui*/
					is_voisin = kECSP_False;
					be_cour = delta_v;
					while (be_cour != NIL_BE)
					{
						if ((be_cour->node2 == KFile[j]) && (be_cour->cap == 1.0))
							is_voisin = kECSP_True;

						be_cour = be_cour->next;
					}

					/*e = Ui->first_edge;
					while((e != NIL_E) && (is_voisin == kECSP_False))
					{
						if( (e->adjac->id == KFile[j]) && (e->X == 1.0) )
							is_voisin = kECSP_True;

						e = e->next;
					}*/

					/*Si Ui et Uj ne sont pas voisin alors on*/
					/*incremente leur nombre de deconnexion.*/
					if (is_voisin == kECSP_False)
					{
						deconnexion[KFile[i] - 1]++;
						deconnexion[KFile[j] - 1]++;
					}
				}

				Delete_b_Edge_Set(&delta_v);
			}

			/*On enleve ensuite le sommet qui a la plus grande deconnexion*/
			max_deconnexion = 0;
			for (i = 0; i < n_file; i++)
			{
				if (deconnexion[KFile[max_deconnexion] - 1] < deconnexion[KFile[i] - 1])
				{
					max_deconnexion = i;
				}
			}

			val_max_deconnexion = deconnexion[KFile[max_deconnexion] - 1];

			/*Si tous les sommets sont connectes entre-eux alors on a une clique.*/
			/*On verifie qu'on ne va pas reduire le graphe e 2 sommets et*/
			/*on verifie donc la cardinalite de la clique.*/
			if (val_max_deconnexion == 0)
			{
				if ((n_file >= k) && (n_p - n_file >= 2))
				{
					/*nb_clique++;*/

					K = NIL_BN;

					/*On recopie la file dans le tableau des cliques*/
					/*et on marque chaque sommet de la clique et ceux*/
					/*qui sont dans leurs listes de contraction.*/
					for (i = 0; i < n_file; i++)
					{
						/*Il faut que le plus petit sommet soit au debut*/
						if (K == NIL_BN)
						{
							K = (b_Node *)malloc(sizeof(b_Node));
							K->id = KFile[i];
							K->next = NIL_BN;
						}
						else
						{
							n_cour = (b_Node *)malloc(sizeof(b_Node));
							n_cour->id = KFile[i];

							if (n_cour->id <= K->id)
							{
								n_cour->next = K;
								K = n_cour;
							}
							else
							{
								n_cour->next = K->next;
								K->next = n_cour;
							}
						}

						marquage_sommet(gr, KFile[i], mark_tab, NODE_MARKED);
						marquage_sommet(gr, KFile[i], marquage_som, 1);

						/*node_cour = KFile[i];
						do
						{
							mark_tab[node_cour-1] = NODE_MARKED;

							node_cour = gr->Nodes[node_cour-1].next_sh->id;
						}
						while(node_cour != KFile[i]);*/
					}

					/*puts("K == ");
					Print_b_Node_Set(K,stdout);*/

					/*clique[nb_clique-1] = K;*/

					/*On verifie la cardinalite de la coupe*/
					delta_w = get_delta_w(gr, K, NO_0_EDGES, CAP_VALUE_X);

					deg = 0;
					be_cour = delta_w;
					while (be_cour != NIL_BE)
					{
						deg++;

						be_cour = be_cour->next;
					}

					if ((k_card_1 <= deg) && (deg <= k_card_2))
					{
						nb_clique++;
						clique[nb_clique - 1] = K;

						clique_trouver = kECSP_True;
					}
					else
					{
						Delete_b_Node_Set(&K);

						clique_trouver = kECSP_False;
					}

					Delete_b_Edge_Set(&delta_w);

					/*clique_trouver = kECSP_True;*/
				}
				else
				{
					clique_trouver = kECSP_False;
				}
			}
			else /*Sinon on supprime de la file le sommet qui est le plus deconnecte*/
			{	 /*et on recompte les deconnexion*/
				for (j = max_deconnexion + 1; j < n_file; j++)
				{
					KFile[j - 1] = KFile[j];
				}

				n_file--;
			}
		} while (val_max_deconnexion != 0);

		if (clique_trouver == kECSP_False)
		{
			marquage_sommet(gr, nptr->id, marquage_som, 1);

			/*node_cour = nptr->id;
			do
			{
				mark_tab[node_cour-1] = NODE_MARKED;

				node_cour = gr->Nodes[node_cour-1].next_sh->id;
			}
			while(node_cour != nptr->id);*/

			/*Recherche d'un nouveau sommet e explorer*/
			new_node = &(gr->Nodes[0]);
			do
			{
				new_node = new_node->next_gr;
			} while ((new_node != &(gr->Nodes[0])) && (marquage_som[new_node->id - 1] != 0));

			if (marquage_som[new_node->id - 1] == 0)
			{
				nptr = new_node;

				nouveau_sommet = kECSP_True;
			}
			else
			{
				nouveau_sommet = kECSP_False;
			}
		}
		else
		{
			/*On s'arrete des qu'on trouve une clique*/
			nouveau_sommet = kECSP_False;
		}
	} while (nouveau_sommet == kECSP_True);

	(*n_clique) = nb_clique;

	free(KFile);
	free(deconnexion);
	free(marquage_som);
}

/*This function checks if all the edges given by Tab*/
/*statisfied X(e)=a.It will return -1 if all the edges*/
/*statisfied the condition X(e) = a and i if the edge ei*/
/*doesn't satisfied this constraint.*/
long find_edge(simpleEdge *Tab, long m_tab, float a)
{
	long i;

	i = m_tab - 1;
	while ((i >= 0) && (Tab[i].cap != a))
	{
		i--;
	}

	return i;
}

/*Create the complement of the set W that we have first built with create_cut_set()*/
/*If b_cour=uv is a bad edge of the ghct, with get_cut_set we create the set W above the*/
/*node above u.The complement of the set W will be the set under the node v.*/
void create_complement_cut_set(Tree_Node *t, b_Node **f_list, b_Edge *b_cour)
{
	b_Node *t_cour;
	Tree_Node *t_son; /*Son of the Tree_Node t*/

	if (t != NIL_TN)
	{
		/*Insertion of t as the following node of f_list in the set W.*/
		/*Insertion at the head of the list.*/
		/*We should be sure that the first node of the set is the node*/
		/*that has the lower Id (usefull for the set contraction).*/
		/*Note that this doesn't mean that we are sorting the set W*/

		t_cour = (b_Node *)malloc(sizeof(b_Node));
		t_cour->id = t->id;

		if (*f_list != NIL_BN)
		{
			if (t_cour->id < (*f_list)->id) /*We do the insertion at the head*/
			{
				t_cour->next = (*f_list);
				(*f_list) = t_cour;
			}
			else /*We do the insertion after the first element*/
			{
				t_cour->next = (*f_list)->next;
				(*f_list)->next = t_cour;
			}
		}
		else
		{
			t_cour->next = (*f_list);
			(*f_list) = t_cour;
		}

		/*Exploration of the sons of the node t in the tree*/
		/*We explore all the tree until we find the bad ghct edge b_cour.*/
		/*If we find the edge then we explore all the remaining tree without*/
		/*the part under the bad edge.*/

		t_son = t->lv; /*we start with the first son of t*/
		while (t_son != NIL_TN)
		{
			create_complement_cut_set(t_son, f_list, b_cour);

			t_son = t_son->lh;
		}
	}
}

/*Permet de marquer un sommet et tous les sommets*/
/*qui sont contractes avec lui. Le marquage se fait*/
/*dans le tableau mark_t*/
void marquage_sommet(Graph *gr, long v, Bool *mark_t, Bool mark_val)
{
	Node *nptr;

	nptr = &(gr->Nodes[v - 1]);
	do
	{
		mark_t[nptr->id - 1] = mark_val;

		nptr = nptr->next_sh;

	} while (nptr != &(gr->Nodes[v - 1]));
}

simpleEdge *GetReducedGraph_eList(Graph *g, long *m_edge)
{
	long first, node_cour, n_contr, v;
	Edge *e;
	long i, k;
	simpleEdge *eList, *e_aux;

	// eList = (simpleEdge *)malloc(g->m_Edges*sizeof(simpleEdge));
	eList = new simpleEdge[g->m_Edges];

	for (i = 0; i < g->n_Nodes; i++)
	{
		mark_tab[i] = kECSP_False;
	}

	k = 0;
	first = 1L;
	node_cour = 1L;
	do
	{
		/*Parcours de tous les sommets qui sont contractes sur node_cour*/
		n_contr = node_cour;

		do
		{
			e = g->Nodes[n_contr - 1].first_edge;
			while (e != NIL_E)
			{
				v = e->adjac->id;

				/*On ne prend pas les aretes e 0*/
				if ((mark_tab[v - 1] == kECSP_False) && (g->Nodes[v - 1].n_sh->id != g->Nodes[n_contr - 1].n_sh->id) && (e->X != 0.0))
				{

					eList[k].node1 = g->Nodes[n_contr - 1].n_sh->id;
					eList[k].node2 = g->Nodes[v - 1].n_sh->id;
					eList[k].cap = e->X;

					k++;
				}

				e = e->next;
			}

			mark_tab[n_contr - 1] = kECSP_True;

			n_contr = g->Nodes[n_contr - 1].next_sh->id;
		} while (n_contr != node_cour);

		node_cour = g->Nodes[node_cour - 1].next_gr->id;
	} while (node_cour != first);

	/*Optimisation de l'espace*/
	e_aux = eList;

	// eList = (simpleEdge *)malloc(k*sizeof(simpleEdge));
	eList = new simpleEdge[k];

	for (i = 0; i < k; i++)
	{
		eList[i].node1 = e_aux[i].node1;
		eList[i].node2 = e_aux[i].node2;
		eList[i].cap = e_aux[i].cap;
	}

	(*m_edge) = k;

	// free(e_aux);
	delete[] e_aux;

	return eList;
}

/*Serch a fractionnary cycle on the graph G*/
/*Option says if we are searching on the*/
/*shrunk graph or on the original graph.*/
/*Si Status == 1 alors on a trouve un cycle*/
/*Si status == 0 alors on n'eu qu'un chemin*/
/*Si status == -1 alors on a explore toutes les aretes fractionnaires*/
b_Node *find_frac_cycle(Graph *gr, long n, long m, long *p, Bool *marquage, Bool *cycle_type, short *status)
{
	b_Node *cycle, **b_cour;
	short *marquer;
	Edge *eptr;
	long i, j;
	long u, u_prec, cycle_head = 0;
	Bool cycle_trouver, is_not_a_cycle;
	Node *nptr;
	long nb_node_marquer;
	b_Node *b_last;
	b_Edge *e_cycle;
	b_Edge **be_cour, *bedge_cour, *dernier_entree;
	Edge *dernier_bon;
	Bool arete_trouvee;
 
	/*Sommets et aretes du cycle*/
	cycle = NIL_BN;
	e_cycle = NIL_BE; 
	/*On commence par rechercher une arete non marquee*/
	i = 0;
	while ((i < m) && ((marquage[gr->Edges[i].num - 1] == 2) || (marquage[gr->Edges[i].num - 1] == 3)))
		i++;

	/*Dans le cas oe on ne trouve pas d'arete*/
	/*alors on renvoie La liste des sommets fractionnaires*/
	if (i >= m)
	{
		(*status) = -1;

		(*cycle_type) = kECSP_False;

		nb_node_marquer = 0;
		for (i = 0; i < n; i++)
		{
			if (gr->Nodes[i].first_edge != NIL_E)
			{
				b_last = (b_Node *)malloc(sizeof(b_Node));
				b_last->id = gr->Nodes[i].id;
				b_last->next = cycle;
				cycle = b_last;

				nb_node_marquer++;
			}
		}

		(*p) = nb_node_marquer;

		return cycle;
	}

	/*Parcours et marquage des voisins de eptr*/
	marquer = (short *)malloc(n * sizeof(short));
	for (j = 0; j < n; j++)
		marquer[j] = 0; 

	u = gr->Edges[i].back->adjac->id;
	marquer[u - 1] = MARKED;

	/*Ajout d'un noeud dans le cycle*/
	b_cour = &cycle;
	(*b_cour) = (b_Node *)malloc(sizeof(b_Node));
	(*b_cour)->id = u;
	(*b_cour)->next = NIL_BN;
	b_cour = &((*b_cour)->next);

	u_prec = u;

	u = gr->Edges[i].adjac->id;
	marquer[u - 1] = MARKED; 

	/*ajout du deuxieme sommet dans le cycle*/
	(*b_cour) = (b_Node *)malloc(sizeof(b_Node));
	(*b_cour)->id = u;
	(*b_cour)->next = NIL_BN;
	b_last = (*b_cour);
	b_cour = &((*b_cour)->next);

	/*ajout de l'arete dans le cycle*/
	be_cour = &e_cycle;

	(*be_cour) = (b_Edge *)malloc(sizeof(b_Edge));
	(*be_cour)->node1 = u_prec;
	(*be_cour)->node2 = u;
	(*be_cour)->cap = gr->Edges[i].X;
	(*be_cour)->num = gr->Edges[i].num;
	(*be_cour)->next = NIL_BE;

	dernier_entree = (*be_cour);

	be_cour = &((*be_cour)->next);

	marquage[gr->Edges[i].num - 1] = 1;

	nptr = &(gr->Nodes[u - 1]);
	nb_node_marquer = 2; /*Correspond au nombre de sommets marques*/
	cycle_trouver = kECSP_False;
	is_not_a_cycle = kECSP_False;
	while ((nb_node_marquer < n) && ((!cycle_trouver) && (!is_not_a_cycle)))
	{
		arete_trouvee = kECSP_False;
		dernier_bon = NIL_E;

		eptr = nptr->first_edge;
		while ((eptr != NIL_E) && (arete_trouvee == kECSP_False))
		{
			/*On prend une arete qui n'a pas encore ete*/
			/*prise dans un cycle ou qui ne mene nul part*/
			if ((eptr->adjac->id != u_prec) && (marquage[eptr->num - 1] != 2) && (marquage[eptr->num - 1] != 3))
			{
				/*On cherche d'abord l'arete qui complete la precedente e 1*/
				if ((eptr->X + dernier_entree->cap) == 1.0)
				{
					arete_trouvee = kECSP_True;
				}
				else
				{
					/*On retient la derniere bonne arete*/
					dernier_bon = eptr;
					eptr = eptr->next;
				}
			}
			else
			{
				eptr = eptr->next;
			}
		} 

		if (arete_trouvee == kECSP_False)
		{
			eptr = dernier_bon;
		}

		if (eptr != NIL_E)
		{
			is_not_a_cycle = kECSP_False;

			u_prec = u;
			u = eptr->adjac->id;

			if (marquer[u - 1] != MARKED)
			{
				marquer[u - 1] = MARKED;

				nptr = &(gr->Nodes[u - 1]);

				nb_node_marquer++;

				/*Ajout d'un noeud dans le cycle*/
				(*b_cour) = (b_Node *)malloc(sizeof(b_Node));
				(*b_cour)->id = u;
				(*b_cour)->next = NIL_BN;
				b_last = (*b_cour); /*On retient le dernier entre dans le cycle*/
				b_cour = &((*b_cour)->next);

				(*be_cour) = (b_Edge *)malloc(sizeof(b_Edge));
				(*be_cour)->node1 = eptr->back->adjac->id;
				(*be_cour)->node2 = eptr->adjac->id;
				(*be_cour)->cap = eptr->X;
				(*be_cour)->num = eptr->num;
				(*be_cour)->next = NIL_BE;

				dernier_entree = (*be_cour);

				be_cour = &((*be_cour)->next);
			}
			else
			{
				cycle_trouver = kECSP_True;
				cycle_head = u; /*On retient la tete du cycle:ce n'est pas*/
								/*forcement le premier de la liste.*/

				(*status) = 1;
			}

			/*On marque l'arete pour etre ser de ne pas la reprendre*/
			marquage[eptr->num - 1] = 1;
		}
		else
		{
			is_not_a_cycle = kECSP_True;

			/*On marque la derniere arete comme sterile*/
			marquage[dernier_entree->num - 1] = 3;

			(*status) = 0;
		}
	}

	if (cycle_trouver == kECSP_True)
	{
		/*Il faut verifier que le cycle*/
		/*ne contient pas de branche isolee.*/
		/*On verifie pour cela si le dernier sommet*/
		/*du cycle est adjacent au premier sommet du cycle.*/
		/*S'ils sont voisins alors le cycle ne contient pas*/
		/*de branche isolee. Sinon on supprime de la liste*/
		/*le dernier sommet du cycle.*/
		/*On supprime aussi l'arete de la liste des aretes*/
		u = cycle->id;
		while (u != cycle_head)
		{
			b_last = cycle;
			cycle = cycle->next;
			free(b_last);

			bedge_cour = e_cycle;
			e_cycle = e_cycle->next;
			free(bedge_cour);

			nb_node_marquer--;

			u = cycle->id;
		}

		/*On marquage toutes les aretes du cycle par 2*/
		bedge_cour = e_cycle;
		while (bedge_cour != NIL_BE)
		{
			marquage[bedge_cour->num - 1] = 2;

			bedge_cour = bedge_cour->next;
		}

		Delete_b_Edge_Set(&e_cycle);
	}
	else /*suppression de la liste chainee*/
	{ 
		Delete_b_Node_Set(&cycle);
		Delete_b_Edge_Set(&e_cycle); 
	}

	/*DeleteGraph(&gr);*/
	free(marquer);

	if (cycle_trouver == kECSP_True)
		(*cycle_type) = kECSP_True;
	else
		(*cycle_type) = kECSP_False;

	(*p) = nb_node_marquer;
	return cycle;
}

/*Recherche des motifs:les motifs sont des triangles*/
void get_sp_partition_chaine(Graph *gr, b_Node **chaine, long *chaine_sz, long *nb_chaine)
{
	b_Node *Triangle[2];
	long i, nb_triangle;
	long u_cour = 0, u = 0;
	b_Node *b_cour;
	b_Edge *delta_v, *be_cour;
	Node *nptr, *new_node;
	long n, m, nb_voisin, prec;
	Bool nouveau_sommet = kECSP_True;
	short *marquage;
	b_Node *dernier;

	n = gr->n_Nodes;
	m = gr->m_Edges;

	marquage = (short *)malloc(n * sizeof(short));

	for (i = 0; i < 2; i++)
	{
		Triangle[i] = NIL_BN;
	}

	for (i = 0; i < n; i++)
	{
		mark_tab[i] = 0;
	}

	/*On commence par contracter tous les triangles*/
	do
	{
		nb_triangle = 0;
		get_clique(gr, 3, 1, 100, Triangle, &nb_triangle);
		contract_set_w(gr, Triangle[0]); 
		Delete_b_Node_Set(&Triangle[0]); 
	} while (nb_triangle >= 1);

	/*Recherche des motifs*/
	/*On met le resultat dans chaine[i]*/
	(*nb_chaine) = 0;
	nptr = &(gr->Nodes[0]);
	do
	{
		for (i = 0; i < n; i++)
		{
			marquage[i] = 0;
		}

		delta_v = get_delta_v(gr, nptr->id, NO_0_EDGES, CAP_VALUE_X);
		be_cour = delta_v;
		while (be_cour != NIL_BE)
		{
			marquage[be_cour->node2 - 1]++;

			be_cour = be_cour->next;
		}

		Delete_b_Edge_Set(&delta_v);

		/*Comptage du nombre de voisin par arete double*/
		nb_voisin = 0;
		for (i = 0; i < n; i++)
		{
			if ((marquage[i] >= 2) && (mark_tab[i] == 0))
			{
				nb_voisin++;
				u_cour = i + 1;
			}
		}

		if (nb_voisin >= 1)
		{ 
			/*U contient le sommet en question*/
			chaine[*nb_chaine] = (b_Node *)malloc(sizeof(b_Node));
			chaine[*nb_chaine]->id = nptr->id;

			chaine[*nb_chaine]->next = (b_Node *)malloc(sizeof(b_Node));
			chaine[*nb_chaine]->next->id = u_cour;

			chaine[*nb_chaine]->next->next = NIL_BN;

			dernier = chaine[*nb_chaine]->next;

			mark_tab[nptr->id - 1] = 1;
			mark_tab[u_cour - 1] = 1;

			chaine_sz[*nb_chaine] = 2;

			prec = nptr->id;
			do
			{
				delta_v = get_delta_v(gr, u_cour, NO_0_EDGES, CAP_VALUE_X);

				for (i = 0; i < n; i++)
				{
					marquage[i] = 0;
				}

				be_cour = delta_v;
				while (be_cour != NIL_BE)
				{
					marquage[be_cour->node2 - 1]++;

					be_cour = be_cour->next;
				}

				Delete_b_Edge_Set(&delta_v);

				/*Comptage du nombre de voisins par arète double*/
				nb_voisin = 0;
				for (i = 0; i < n; i++)
				{
					if ((marquage[i] >= 2) && (mark_tab[i] == 0))
					{
						nb_voisin++;
						u = i + 1; 
					}
				}

				if (nb_voisin >= 1)
				{ 
					b_cour = (b_Node *)malloc(sizeof(b_Node));
					b_cour->id = u;

					/*Il faut inserer en fin de liste*/
					dernier->next = b_cour;
					b_cour->next = NIL_BN;

					dernier = b_cour; 
					mark_tab[u - 1] = 1;

					/*Longueur de la chaine: pour calculer les RHS*/
					chaine_sz[*nb_chaine]++;
				}

				prec = u_cour;
				u_cour = u;
			} while (nb_voisin != 0); 
			(*nb_chaine)++;
		}

		/*Recherche d'un nouveau sommet*/
		new_node = nptr;
		do
		{
			new_node = new_node->next_gr;
		} while ((new_node != &(gr->Nodes[0])) && (mark_tab[new_node->id - 1] != 0));

		/*nptr = nptr->next_gr;*/
		if (new_node != &(gr->Nodes[0]))
		{
			nouveau_sommet = kECSP_True;

			nptr = new_node;
		}
		else
		{
			nouveau_sommet = kECSP_False;
		}

	} while (nouveau_sommet == kECSP_True);

	free(marquage);

/*	std::cout <<"CH == ";
	for(i=0;i<(*nb_chaine);i++)
	{
		std::cout << " ( ";
		Print_b_Node_Set(chaine[i],stdout);
		std::cout << " ) ";
	}*/
}

/*Calculate the connected components of a graph using X or cap as capcity*/
/*This function is not used in the reduction operation.It's used*/
/*in the Abacus program while checking the connexity of a solution X*/
b_Node **graph_connected_components(Graph *gr, long *n_comp, int edge_flag)
{
	long non_explorer, i, u;
	Fifo *fifo;
	Bool *marquer, *explorer;
	long node_cour;
	Edge *e;
	long n, m;
	b_Node *new_node;
	b_Node **comp_tab;

	n = gr->n_Nodes;
	m = gr->m_Edges;

	comp_tab = new b_Node *[(n / 2)];

	/*marquer = (Bool *)malloc(n*sizeof(Bool));*/
	marquer = new Bool[n];

	if (marquer == NULL)
	{
		fprintf(stdout, "%s", "ERROR:Connected components.Unable to allocate table MARQUE.");
		return kECSP_False;
	}

	for (i = 0; i < n; i++)
	{
		marquer[i] = kECSP_False;
	}

	/*explorer = (Bool *)malloc(n*sizeof(Bool));*/
	explorer = new Bool[n];
	if (explorer == NULL)
	{
		fprintf(stdout, "%s", "ERROR:Connected components.Unable to allocate table EXPLORER.");

		/*free(marquer);*/
		delete[] marquer;
		return kECSP_False;
	}

	for (i = 0; i < n; i++)
	{
		explorer[i] = kECSP_False;
	}

	(*n_comp) = 0;

	fifo_init(&fifo, n);

	non_explorer = 1L;

	fifo_insert(fifo, 1L);

	while (non_explorer <= n)
	{
		node_cour = fifo_get_first(fifo);

		comp_tab[*n_comp] = NIL_BN;

		while (node_cour != -1L)
		{
			i = node_cour - 1;
			explorer[i] = kECSP_True;

			if (marquer[i] == kECSP_False)
			{
				marquer[i] = kECSP_True;

				/*Ajout d'une composante connexe*/
				// printf("U == %d\n",node_cour);

				/*new_node = (b_Node *)malloc(sizeof(b_Node));

				new_node->id = node_cour;

				new_node->next = comp_tab[*n_comp];
				comp_tab[*n_comp] = new_node;*/

				add_to_connected_comp(comp_tab, n_comp, node_cour);
			}

			/*Parcours de la liste d'incidence*/
			e = gr->Nodes[i].first_edge;

			while (e != NIL_E)
			{
				switch (edge_flag)
				{
				case ALL_EDGES:
					u = e->adjac->id;

					if (marquer[u - 1] == kECSP_False)
					{
						marquer[u - 1] = kECSP_True;
						fifo_insert(fifo, u);
						add_to_connected_comp(comp_tab, n_comp, u);
					}

					break;

				case NO_0_EDGES: /*Sans les aretes de capacite cap_value*/
					if (e->X != 0.0)
					{
						u = e->adjac->id;

						if (marquer[u - 1] == kECSP_False)
						{
							// printf("U == %d\n",u);

							marquer[u - 1] = kECSP_True;
							fifo_insert(fifo, u);

							/*new_node = (b_Node *)malloc(sizeof(b_Node));

							new_node->id = u;

							new_node->next = comp_tab[*n_comp];
							comp_tab[*n_comp] = new_node;*/

							add_to_connected_comp(comp_tab, n_comp, u);
						}
					}

					break;

				default: /*We ommit the edges with X(e)=0 in the default case*/
					if (e->X != 0.0)
					{
						u = e->adjac->id;

						if (marquer[u - 1] == kECSP_False)
						{
							marquer[u - 1] = kECSP_True;
							fifo_insert(fifo, u);
							add_to_connected_comp(comp_tab, n_comp, u);
						}
					}

					break;
				}

				e = e->next;
			}

			node_cour = fifo_get_first(fifo);
		}

		/*When the queue is empty then we have build a connected component*/
		/*We search for the last node that hasn't been explored yet.*/
		while ((non_explorer <= n) && (explorer[non_explorer - 1] == kECSP_True))
		{
			non_explorer++;
		}

		fifo_insert(fifo, non_explorer);
		(*n_comp)++;
	}

	/*free(marquer);
	free(explorer);*/
	delete[] marquer;
	delete[] explorer;
	fifo_delete_fifo(&fifo);

	/*puts("ccccc");
	printf("NCOMP == %d\n",*n_comp);
	Print_b_Node_Set(comp_tab[0],stdout);
	puts("dddd");*/

	return comp_tab;
}

Bool separation_sp_partition_2(Graph *gr, b_Node *chaine[], long *chaine_sz, long nb_chaine, long ***sp_part_list, 
long *sp_p, long *nb_sp_part, long *rhs)
{
	b_Edge *sp_part;
	short *mark_edge;
	long i, j, k, n, m;
	b_Edge **delta_Vi;
	b_Edge *be_cour;

	b_Node *b_cour, *b_cour2;
	b_Edge *nouveau;
	Node *nptr;
	long *renum1, *renum2;
	Bool trouver;
	long **pcc;
	double X_D;
	long rhs_sp;
	long n1, n2;
	Bool resultat, is_serie_paralelle;
	b_Edge *sp_mop;
	long **partition;

	n = gr->n_Nodes;
	m = gr->m_Edges;

	(*nb_sp_part) = 0;

	mark_edge = (short *)malloc(m * sizeof(short));

	for (j = 0; j < nb_chaine; j++)
	{
		delta_Vi = (b_Edge **)malloc(chaine_sz[j] * sizeof(b_Edge *));

		for (i = 0; i < n; i++)
		{
			mark_tab[i] = 0;
		}

		/*On commence par marquer les sommets de la chaine*/
		b_cour = chaine[j];
		while (b_cour != NIL_BN)
		{
			mark_tab[b_cour->id - 1] = 1;

			b_cour = b_cour->next;
		}

		/*Calcul des delta_Vi*/
		k = 0;
		b_cour = chaine[j];
		while (b_cour != NIL_BN)
		{
			delta_Vi[k] = get_delta_v(gr, b_cour->id, ALL_EDGES, CAP_VALUE_X);

			k++;
			b_cour = b_cour->next;
		}

		/*Il faut verifier que le graphe support est serie-paralelle*/
		if (chaine_sz[j] >= 3)
		{

			is_serie_paralelle = kECSP_True;
			b_cour = chaine[j];
			k = 0;
			while ((k < chaine_sz[j] - 1) && (is_serie_paralelle))
			{
				b_cour2 = b_cour->next->next;
				while ((b_cour2 != NIL_BN) && (is_serie_paralelle))
				{
					be_cour = delta_Vi[k];
					while ((be_cour != NIL_BE) && ((be_cour->cap == 0) || ((be_cour->cap > 0) && (be_cour->node2 != b_cour2->id))))
					{
						be_cour = be_cour->next;
					}

					if (be_cour != NIL_BE)
					{
						is_serie_paralelle = kECSP_False;
					}

					b_cour2 = b_cour2->next;
				}

				b_cour = b_cour->next;
				k++;
			}
		}
		else
		{
			is_serie_paralelle = kECSP_True;
		}

		if (is_serie_paralelle == kECSP_True)
		{
			/*Calcul de la liste des aretes*/
			for (k = 0; k < m; k++)
			{
				mark_edge[k] = 0;
			}

			sp_part = NIL_BE;
			for (k = 0; k < chaine_sz[j]; k++)
			{
				be_cour = delta_Vi[k];
				while (be_cour != NIL_BE)
				{
					if (mark_edge[be_cour->num - 1] == 0)
					{
						nouveau = (b_Edge *)malloc(sizeof(b_Edge));
						nouveau->node1 = be_cour->node1;
						nouveau->node2 = be_cour->node2;
						nouveau->cap = be_cour->cap;
						nouveau->num = be_cour->num;
						nouveau->next = sp_part;

						sp_part = nouveau;

						mark_edge[be_cour->num - 1] = 1;
					}

					be_cour = be_cour->next;
				}
			}

			/*suppression de la liste des aretes*/
			/*for(i=0;i<chaine_sz[j];i++)
			{
				Delete_b_Edge_Set(&delta_Vi[i]);
			}
			free(delta_Vi);*/

			/*Calcul des coef*/
			/*p >= 4. Donc il y aura des coef 2*/
			/*p == 3. Tous les coef sont e 1*/
			if (chaine_sz[j] >= 3)
			{
				/*Renumerotation de la liste d'arete*/
				renum1 = (long *)malloc(n * sizeof(long));
				renum2 = (long *)malloc((chaine_sz[j] + 1) * sizeof(long));

				b_cour = chaine[j];
				k = 0;
				while (b_cour != NIL_BN)
				{
					k++;
					renum1[b_cour->id - 1] = k;
					renum2[k - 1] = b_cour->id;

					b_cour = b_cour->next;
				}

				/*MAJ de la numerotation de tous les sommets de W*/
				k++;

				nptr = &(gr->Nodes[0]);
				trouver = kECSP_False;
				do
				{
					if (mark_tab[nptr->id - 1] == 0)
					{
						renum1[nptr->id - 1] = k;
						renum2[k - 1] = nptr->id;
					}

					nptr = nptr->next_gr;

				} while (nptr != &(gr->Nodes[0]));

				/*Calcul des aretes faisant partie du MOP: y compris les aretes e 0*/

				sp_mop = NIL_BE;
				be_cour = sp_part;
				while (be_cour != NIL_BE)
				{
					/*On ne prend pas une arete entre les sommets i et i+2 au moins*/
					if (((renum1[be_cour->node1 - 1] != k) && (renum1[be_cour->node2 - 1] != k) &&
						 (abs(renum1[be_cour->node2 - 1] - renum1[be_cour->node1 - 1]) <= 1)) ||
						(renum1[be_cour->node1 - 1] == k) || (renum1[be_cour->node2 - 1] == k))
					{
						nouveau = (b_Edge *)malloc(sizeof(b_Edge));

						nouveau->node1 = renum1[be_cour->node1 - 1];
						nouveau->node2 = renum1[be_cour->node2 - 1];
						nouveau->cap = be_cour->cap;
						nouveau->num = be_cour->num;

						nouveau->next = sp_mop;
						sp_mop = nouveau;
					}

					be_cour->node1 = renum1[be_cour->node1 - 1];
					be_cour->node2 = renum1[be_cour->node2 - 1];

					be_cour = be_cour->next;
				}

				/*puts("MOP == ");
				Print_b_Edge_Set(sp_mop,stdout);*/

				pcc = Floyd(sp_mop, k);

				/*Calcul de X_D*/
				X_D = 0;
				be_cour = sp_part;
				while (be_cour != NIL_BE)
				{
					n1 = be_cour->node1;
					n2 = be_cour->node2;

					/*n1 = renum2[be_cour->node1-1];
					n2 = renum2[be_cour->node2-1];*/

					X_D = X_D + pcc[n1 - 1][n2 - 1] * be_cour->cap;

					/*On stocke le coef de l'arete dans son champ cap*/
					be_cour->cap = pcc[n1 - 1][n2 - 1];

					be_cour = be_cour->next;
				}

				free(renum1);
				free(renum2);

				Delete_b_Edge_Set(&sp_mop);

				for (i = 0; i < k; i++)
				{
					free(pcc[i]);
				}
				free(pcc);
			}
			else
			{
				/*Calcul de X_D*/
				X_D = 0;
				be_cour = sp_part;
				while (be_cour != NIL_BE)
				{
					X_D = X_D + be_cour->cap;

					be_cour = be_cour->next;
				}
			}

			/*Calcul du RHS*/
			rhs_sp = 0;//sup_part(k_ordre, 2) * (chaine_sz[j] + 1) - 1;

			// printf("RHS == %li  X_D %f\n",rhs_sp,X_D);

			// Print_b_Edge_Set(sp_part,stdout);

			//if (X_D < (rhs_sp - EPSILON))
			{
				(*sp_part_list)[*nb_sp_part] = new long[n];

				for (i = 0; i < n; i++)
				{
					(*sp_part_list)[*nb_sp_part][i] = chaine_sz[j] + 1;
				}

				i = 0;
				b_cour = chaine[j];
				while (b_cour != NIL_BN)
				{
					nptr = &(gr->Nodes[b_cour->id - 1]);
					do
					{
						(*sp_part_list)[*nb_sp_part][nptr->id - 1] = i + 1;

						nptr = nptr->next_sh;

					} while (nptr != &(gr->Nodes[b_cour->id - 1]));

					i++;
					b_cour = b_cour->next;
				}

				Delete_b_Edge_Set(&sp_part);

				// puts("SP-CHAINE violee");
				rhs[*nb_sp_part] = rhs_sp;
				sp_p[*nb_sp_part] = chaine_sz[j] + 1;

				(*nb_sp_part)++;
			}
			/*else
			{
				Delete_b_Edge_Set(&sp_part);

				// puts("SP-CHAINE non violee");
			}*/
		}

		/*suppression de la liste des aretes*/
		for (i = 0; i < chaine_sz[j]; i++)
		{
			Delete_b_Edge_Set(&delta_Vi[i]);
		}
		free(delta_Vi);
	}

	if ((*nb_sp_part) >= 1)
	{
		resultat = kECSP_True;
	}
	else
	{
		resultat = kECSP_False;
	}

	free(mark_edge);

	return resultat;
}

/*Calcule le plus court chemin entre*/
/*tous les couples de sommets dans le graphe G*/
/*On utilise l'algorithme de Floyd*/
long **Floyd(b_Edge *sh_eList, long n)
{
	long **pcc;
	long i, j, k;
	long n1, n2;
	b_Edge *be_cour;

	// printf("N == %d\n",n);

	/*Initialisation du tableau des pcc*/
	pcc = (long **)malloc(n * sizeof(long *));

	for (i = 0; i < n; i++)
	{
		pcc[i] = (long *)malloc(n * sizeof(long));
	}

	/*Initialisation du tableau des pcc*/
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			pcc[i][j] = INFINI;
		}
	}

	/*for(i=0;i<m;i++)
	{
		n1 = sh_eList[i].node1;
		n2 = sh_eList[i].node2;

		pcc[n1-1][n2-1] = 1;
		pcc[n2-1][n1-1] = 1;
	}*/

	be_cour = sh_eList;
	while (be_cour != NIL_BE)
	{
		n1 = be_cour->node1;
		n2 = be_cour->node2;

		pcc[n1 - 1][n2 - 1] = 1;
		pcc[n2 - 1][n1 - 1] = 1;

		be_cour = be_cour->next;
	}

	for (i = 0; i < n; i++)
	{
		pcc[i][i] = 0;
	}

	/*nptr1 = &(gr->Nodes[0]);
	do
	{
		k = nptr1->id-1;

		nptr2 = &(gr->Nodes[0]);
		do
		{
			i = nptr2->id-1;

			nptr3 = &(gr->Nodes[0]);
			do
			{
				j = nptr3->id-1;

				pcc[i][j] = min(pcc[i][j],pcc[i][k] + pcc[k][j]);

				nptr3 = nptr3->next_gr;

			}while(nptr3 != &(gr->Nodes[0]));

			nptr2 = nptr2->next_gr;

		}while(nptr2 != &(gr->Nodes[0]));

		nptr1 = nptr1->next_gr;

	}while(nptr1 != &(gr->Nodes[0]));*/

	for (k = 0; k < n; k++)
	{
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				pcc[i][j] = (long)min(pcc[i][j], pcc[i][k] + pcc[k][j]);
			}
		}
	}

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			// printf("pcc[%d][%d] == %d\n",i,j,pcc[i][j]);
		}
	}

	/*nptr1 = &(gr->Nodes[0]);
	do
	{
		i = nptr1->id-1;

		nptr2 = &(gr->Nodes[0]);
		do
		{
			j = nptr2->id-1;

			printf("i = %d  j = %d  %d\n",i+1,j+1,pcc[i][j]);

			nptr2 = nptr2->next_gr;

		}while(nptr2 != &(gr->Nodes[0]));

		nptr1 = nptr1->next_gr;

	}while(nptr1 != &(gr->Nodes[0]));*/

	return pcc;
}

/*Permet de reinitialiser la table des contractions*/
/*en rappelant la derniere contraction sauvegardee*/
void recall_contraction_info(Graph *gr, int niveau)
{
	long n, i;

	n = gr->n_Nodes;

	if (niveau == 1)
	{
		for (i = 0; i < n; i++)
		{
			gr->Nodes[i].next_gr = gr->Nodes[i].a_next_gr;
			gr->Nodes[i].n_sh = gr->Nodes[i].a_n_sh;
			gr->Nodes[i].next_sh = gr->Nodes[i].a_next_sh;
		}
	}
	else if (niveau == 2)
	{
		for (i = 0; i < n; i++)
		{
			gr->Nodes[i].next_gr = gr->Nodes[i].b_next_gr;
			gr->Nodes[i].n_sh = gr->Nodes[i].b_n_sh;
			gr->Nodes[i].next_sh = gr->Nodes[i].b_next_sh;
		}
	}
}

/*Separation des F-partition*/
/*Retourne Vrai si F-Partition trouvee et Faux sinon*/
Bool separation_f_partition(Graph *gr, int k, b_Node *frac_cycle, long cycle_sz, long **f_partition, long **V0_F, long *V0_F_sz, long *rhs, FILE *sortie_frac_gk)
{
	long nb_fp;
	b_Node *V0;
	Bool *mark_tab2;
	long rhs_f, f_card, delta_V0_card;
	b_Edge **delta_Vi;
	b_Edge *F_part, *be_cour, *fe_cour;
	float X_F, X_D, X_D0;
	b_Node *b_cour;
	long i, s, j;
	b_Node *b_vi;
	b_Edge *f_part_list = NIL_BE, *e_aux, *D0_frac_edge;
	long n, m;
	Bool resultat = kECSP_False;
	Node *nptr;
	b_Edge **be_ptr;
	double minpoids;

	// puts("F-part");

	n = gr->n_Nodes;
	m = gr->m_Edges;

	V0 = NIL_BN;
	nb_fp = 0;

	/*Calacul de V0:toujours dans le graphe reduit*/
	// mark_tab2 = (Bool *)malloc(n*sizeof(Bool));
	mark_tab2 = new Bool[n];

	for (i = 0; i < n; i++)
	{
		mark_tab2[i] = kECSP_False;
	}

	/*Marquage des sommets du cycle*/
	b_cour = frac_cycle;
	double sum_Con = 0.0;
	while (b_cour != NIL_BN)
	{
		sum_Con += (gr)->con[b_cour->id - 1];
		mark_tab2[b_cour->id - 1] = kECSP_True;

		b_cour = b_cour->next;
	}
 

	/*Parcours du graphe reduit*/
	V0 = NIL_BN;
	nptr = &(gr->Nodes[0]);
	do
	{
		if (mark_tab2[nptr->id - 1] == kECSP_False)
		{ 
			b_cour = (b_Node *)malloc(sizeof(b_Node));
			b_cour->id = nptr->id;
			b_cour->next = V0;
			V0 = b_cour;
		}

		nptr = nptr->next_gr;

	} while (nptr != &(gr->Nodes[0]));

	// free(mark_tab2);
	delete[] mark_tab2;

	/*Si V0 est vide alors ce n'est pas la peine de continuer*/
	if (V0 == NIL_BN)
	{
		return kECSP_False;
	}

	/*********************/
	/*printf("%s","V0 == ");
	Print_b_Node_Set(V0,stdout);
	puts("");
	Print_b_Node_Set(frac_cycle,stdout);
	puts("\niiiiiiiiiiiii");*/
	/********************/

	/*Calcul des delta(Vi):on prend en compte les aretes nulles*/
	// delta_Vi = (b_Edge **)malloc((cycle_sz+1)*sizeof(b_Edge *));
	delta_Vi = new b_Edge *[cycle_sz + 1];
	delta_Vi[0] = get_delta_w(gr, V0, ALL_EDGES, CAP_VALUE_X);

	i = 1;
	b_cour = frac_cycle;
	while (b_cour != NIL_BN)
	{
		/*b_vi = (b_Node *)malloc(sizeof(b_Node));
		b_vi->id = b_cour->id;
		b_vi->next = NIL_BN;

		delta_Vi[i] = get_delta_w(gr,b_vi,ALL_EDGES,CAP_VALUE_X);

		free(b_vi);*/

		delta_Vi[i] = get_delta_v(gr, b_cour->id, ALL_EDGES, CAP_VALUE_X);

		i++;
		b_cour = b_cour->next;
	}

	/*Calcul de X_D*/
	X_D = 0;
	for (s = 1; s <= cycle_sz; s++)
	{
		be_cour = delta_Vi[s];
		while (be_cour != NIL_BE)
		{
			X_D = X_D + be_cour->cap;

			be_cour = be_cour->next;
		}
	}

	/*Recherche du F*/
	F_part = NIL_BE;
	f_card = 0;
	delta_V0_card = 0;
	be_cour = delta_Vi[0];
	X_D0 = 0;
	X_F = 0;
	D0_frac_edge = NIL_BE;
	minpoids = 100;
	while (be_cour != NIL_BE)
	{
		delta_V0_card++;
		X_D0 = X_D0 + be_cour->cap;

		/*On commence par prendre toutes les aretes e 1*/
		if (be_cour->cap > 0)
		{
			fe_cour = (b_Edge *)malloc(sizeof(b_Edge));
			fe_cour->node1 = be_cour->node1;
			fe_cour->node2 = be_cour->node2;
			fe_cour->cap = be_cour->cap;
			fe_cour->num = be_cour->num;

			fe_cour->next = F_part;
			F_part = fe_cour;

			f_card++;
			X_F = X_F + fe_cour->cap;

			if (fe_cour->cap < minpoids)
			{
				D0_frac_edge = fe_cour;
				minpoids = fe_cour->cap;
			}
		}

		be_cour = be_cour->next;
	}

	/*Verification de la parite de F et p:*/
	/*S'ils sont de meme parite alors on augmente ou diminue F de 1*/
	if ((((is_pair(k * cycle_sz) == kECSP_True) && (is_pair(f_card) == kECSP_True))) ||
		(((is_pair(k * cycle_sz) == kECSP_False) && (is_pair(f_card) == kECSP_False))))
	{
		if ((D0_frac_edge != NIL_BE) && (f_card >= 1))
		{
			be_ptr = &F_part;
			while ((*be_ptr) != D0_frac_edge)
			{
				be_ptr = &((*be_ptr)->next);
			}

			(*be_ptr) = (*be_ptr)->next;

			f_card--;
			X_F = X_F - D0_frac_edge->cap;

			free(D0_frac_edge);
		}
	}

	/*puts("CYCLE == ");
	Print_b_Node_Set(frac_cycle,stdout);

	puts("F == ");
	Print_b_Edge_Set(F_part,stdout);

	puts("");*/

	/*Calcul du RHS*/
	rhs_f = sum_Con - f_card;
	rhs_f = sup_part(rhs_f, 2L);

	/*puts("V0 == ");
	Print_b_Node_Set(V0,stdout);
	fprintf(stdout,"P == %ld |F| == %ld  |V0| == %ld RHS ==	%d\n",cycle_sz,f_card,delta_V0_card,rhs_f);
	fprintf(stdout,"X_D == %f  X_D0 == %f   X_F == %f\n",X_D,X_D0,X_F);*/

	resultat = kECSP_False;

	/*On verifie si la F-partition est violee*/
	if ((((is_pair(k * cycle_sz) == kECSP_True) && (is_pair(f_card) == kECSP_False))) ||
		(((is_pair(k * cycle_sz) == kECSP_False) && (is_pair(f_card) == kECSP_True))))
	{
		/*Verification de l'inegalite*/
		if (((X_D + X_D0) / 2.0 - X_F) < rhs_f - EPSILON)
		{

			resultat = kECSP_True;
			// fputs("F-Partition violee\n",sortie_frac_gk);

			(*f_partition) = new long[n];
			for (i = 0; i < n; i++)
			{
				(*f_partition)[i] = 0;
			}

			i = 0;
			b_cour = V0;
			while (b_cour != NIL_BN)
			{
				nptr = &(gr->Nodes[b_cour->id - 1]);
				do
				{
					(*f_partition)[nptr->id - 1] = i;

					nptr = nptr->next_sh;

				} while (nptr != &(gr->Nodes[b_cour->id - 1]));

				b_cour = b_cour->next;
			}

			i++;

			b_cour = frac_cycle;
			while (b_cour != NIL_BN)
			{
				nptr = &(gr->Nodes[b_cour->id - 1]);
				do
				{
					(*f_partition)[nptr->id - 1] = i;

					nptr = nptr->next_sh;

				} while (nptr != &(gr->Nodes[b_cour->id - 1]));

				b_cour = b_cour->next;
				i++;
			}

			/*Copy des aretes de F*/

			(*V0_F_sz) = 0;
			be_cour = F_part;
			while (be_cour != NIL_BE)
			{
				(*V0_F_sz)++;

				be_cour = be_cour->next;
			}

			(*V0_F) = new long[(*V0_F_sz)];

			j = 0;
			be_cour = F_part;
			while (be_cour != NIL_BE)
			{
				(*V0_F)[j] = be_cour->num;
				j++;

				be_cour = be_cour->next;
			}

			/*for(i=0;i<n;i++)
			printf("P[%d] == %d\n",i,(*f_partition)[i]);

			for(i=0;i<(*V0_F_sz);i++)
			printf("F[%d] == %d\n",i,(*V0_F)[i]);*/
		}
		else
		{
			/*Sinon on essaie de rajouter e F deux aretes*/

			// fputs("F-Partition non-violee\n",sortie_frac_gk);
			resultat = kECSP_False;
		}
	}

	/*Desallocation de la memoire*/
	Delete_b_Node_Set(&V0);
	Delete_b_Edge_Set(&F_part);

	for (i = 0; i < cycle_sz + 1; i++)
	{
		Delete_b_Edge_Set(&delta_Vi[i]);
	}

	// free(delta_Vi);
	delete[] delta_Vi;

	if (resultat == kECSP_True)
	{
		/*puts("TRUE");
		Print_b_Edge_Set(f_part_list,stdout);*/

		(*rhs) = rhs_f;
	}

	return resultat;
}

/*Separation des Partition*/
/*Retourne Vrai si Partition trouvee et Faux sinon*/
Bool separation_partition(Graph *gr, int k, b_Node *frac_cycle, long cycle_sz, long **partition, FILE *sortie_frac_gk)
{
	b_Node *V0;
	Bool *mark_tab2;
	long rhs_p, p_card;
	b_Edge **delta_Vi;
	b_Edge *be_cour;
	float X_D;
	b_Node *b_cour;
	long i, s;
	b_Node *b_vi;
	long n, m;
	Bool resultat = kECSP_False;
	Node *nptr;

	/*puts("Part");*/

	n = gr->n_Nodes;
	m = gr->m_Edges;

	V0 = NIL_BN;

	/*Calacul de V0:toujours dans le graphe reduit*/
	// mark_tab2 = (Bool *)calloc(n,sizeof(Bool));
	mark_tab2 = new Bool[n];

	for (i = 0; i < n; i++)
	{
		mark_tab2[i] = kECSP_False;
	}

	/*Marquage des sommets du cycle*/
	b_cour = frac_cycle;
	while (b_cour != NIL_BN)
	{
		mark_tab2[b_cour->id - 1] = kECSP_True;

		b_cour = b_cour->next;
	}

	/*Parcours du graphe reduit*/
	V0 = NIL_BN;
	nptr = &(gr->Nodes[0]);
	do
	{
		if (mark_tab2[nptr->id - 1] == kECSP_False)
		{
			b_cour = (b_Node *)malloc(sizeof(b_Node));
			b_cour->id = nptr->id;
			b_cour->next = V0;
			V0 = b_cour;
		}

		nptr = nptr->next_gr;

	} while (nptr != &(gr->Nodes[0]));

	// free(mark_tab2);
	delete[] mark_tab2;

	/*Calcul des delta(Vi):on prend en compte les aretes nulles*/
	if (V0 != NIL_BN)
	{
		p_card = cycle_sz + 1;

		// delta_Vi = (b_Edge **)malloc((p_card)*sizeof(b_Edge *));
		delta_Vi = new b_Edge *[p_card];
		delta_Vi[0] = get_delta_w(gr, V0, ALL_EDGES, CAP_VALUE_X);
		i = 1;
	}
	else
	{
		p_card = cycle_sz;

		// delta_Vi = (b_Edge **)malloc((p_card)*sizeof(b_Edge *));
		delta_Vi = new b_Edge *[p_card];
		i = 0;
	}

	/*Calcul des delta(Vi):on prend en compte les aretes nulles*/
	b_cour = frac_cycle;
	while (b_cour != NIL_BN)
	{
		/*b_vi = (b_Node *)malloc(sizeof(b_Node));
		b_vi->id = b_cour->id;
		b_vi->next = NIL_BN;

		delta_Vi[i] = get_delta_w(gr,b_vi,ALL_EDGES,CAP_VALUE_X);

		free(b_vi);*/

		delta_Vi[i] = get_delta_v(gr, b_cour->id, ALL_EDGES, CAP_VALUE_X);

		i++;
		b_cour = b_cour->next;
	}

	/*Calcul de X_D*/
	X_D = 0;
	for (s = 0; s < p_card; s++)
	{
		be_cour = delta_Vi[s];
		while (be_cour != NIL_BE)
		{
			X_D = X_D + be_cour->cap;

			be_cour = be_cour->next;
		}
	}

	/*Calcul du RHS*/
	/*rhs_p = k * p_card;
	rhs_p = sup_part(rhs_p, 2L);*/

	/*Verification de la parite de p:il faut que p soit impair*/
	/*puts("V0 == ");
	Print_b_Node_Set(V0,stdout);
	fprintf(sortie_frac_gk,"P == %ld RHS == %ld  X_D == %f\n",p_card,rhs_p,X_D);*/

	resultat = kECSP_False;

	if (is_pair(k * p_card) == kECSP_False)
	{
		/*Verification de l'inegalite*/
		if ((X_D / 2.0) < rhs_p - EPSILON)
		{
			resultat = kECSP_True;
			// fputs("Partition violee\n",sortie_frac_gk);

			(*partition) = new long[n];
			for (i = 0; i < n; i++)
			{
				(*partition)[i] = 0;
			}

			i = 0;
			b_cour = frac_cycle;
			while (b_cour != NIL_BN)
			{
				nptr = &(gr->Nodes[b_cour->id - 1]);
				do
				{
					(*partition)[nptr->id - 1] = i + 1;

					nptr = nptr->next_sh;

				} while (nptr != &(gr->Nodes[b_cour->id - 1]));

				b_cour = b_cour->next;
				i++;
			}

			b_cour = V0;
			while (b_cour != NIL_BN)
			{
				nptr = &(gr->Nodes[b_cour->id - 1]);
				do
				{
					(*partition)[nptr->id - 1] = i + 1;

					nptr = nptr->next_sh;

				} while (nptr != &(gr->Nodes[b_cour->id - 1]));

				b_cour = b_cour->next;
			}
		}
		else
		{
			// fputs("Partition non-violee\n",sortie_frac_gk);
			resultat = kECSP_False;
		}
	}

	/*Desallocation de la memoire*/
	Delete_b_Node_Set(&V0);

	for (i = 0; i < p_card; i++)
	{
		Delete_b_Edge_Set(&delta_Vi[i]);
	}

	// free(delta_Vi);
	delete[] delta_Vi;

	if (resultat == kECSP_True)
	{
		//(*rhs) = rhs_p;
	}

	return resultat;
}

/*Retourne Vrai si l'entier z est pair et Faux sinon*/
Bool is_pair(long z)
{
	ldiv_t parite;

	parite = ldiv(z, 2);

	if (parite.rem == 1)
		return kECSP_False;
	else
		return kECSP_True;
}

/*Recursive function for searching a violated cut in the graph.*/
/*It needs first to calculate the ghct associated to the graph.*/
/*b_cour is the ghct edge that doesn't satisfied the criterion*/
/*we have fixed and that we shouldn't explore.*/
/*f_list is first b_Node of the final list that will represent*/
/*the set W.*/
void create_cut_set(Tree_Node *t, b_Node **f_list, b_Edge *b_cour, Bool *found)
{
	b_Node *t_cour;
	Tree_Node *t_son; /*Son of the Tree_Node t*/

	if (t != NIL_TN)
	{
		/*Insertion of t as the following node of f_list in the set W.*/
		/*Insertion at the head of the list.*/
		/*We should be sure that the first node of the set is the node*/
		/*that has the lower Id (usefull for the set contraction).*/
		/*Note that it's not necessary that the list is sorted.*/
		t_cour = (b_Node *)malloc(sizeof(b_Node));
		t_cour->id = t->id;

		if (*f_list != NIL_BN)
		{
			if (t_cour->id < (*f_list)->id) /*We do the insertion at the head*/
			{
				t_cour->next = (*f_list);
				(*f_list) = t_cour;
			}
			else /*We do the insertion after the first element*/
			{
				t_cour->next = (*f_list)->next;
				(*f_list)->next = t_cour;
			}
		}
		else
		{
			t_cour->next = (*f_list);
			(*f_list) = t_cour;
		}

		/*We explore all the tree until we find the bad ghct edge b_cour.*/
		/*If we find the edge then we explore all the remaining tree without*/
		/*the part under the bad edge.*/

		t_son = t->lv; /*we start with the first son of t*/
		while (t_son != NIL_TN)
		{
			/*if we find the bad edge we don't explore the following part of the tree*/
			/*else we explore the remaining part of the tree*/
			if (t_son->id != b_cour->node2)
			{
				create_cut_set(t_son, f_list, b_cour, found);
			}
			else
			{
				if (!(*found))
				{
					(*found) = kECSP_True;
				}
				else
				{
					create_cut_set(t_son, f_list, b_cour, found);
				}
			}

			t_son = t_son->lh;
		}
	}
}

void add_to_connected_comp(b_Node **tab_comp, long *n_comp, long u)
{
	b_Node *new_node;

	new_node = (b_Node *)malloc(sizeof(b_Node));

	new_node->id = u;

	/*On insere l'element u de sorte que le premier sommet*/
	/*de la liste soit celui avec le plus petit incide.*/
	if (tab_comp[*n_comp] != NIL_BN)
	{
		if (tab_comp[*n_comp]->id > u)
		{
			/*Dans ce cas on fait une insertion en tete*/
			new_node->next = tab_comp[*n_comp];
			tab_comp[*n_comp] = new_node;
		}
		else
		{
			/*Dans ce cas on insere le nouveau noeud apres le premier noeud*/
			new_node->next = tab_comp[*n_comp]->next;
			tab_comp[*n_comp]->next = new_node;
		}
	}
	else
	{
		/*Dans ce cas on fait une insertion en tete*/
		new_node->next = tab_comp[*n_comp];
		tab_comp[*n_comp] = new_node;
	}
}

long sup_part(long num, long denom)
{
	ldiv_t parite;

	parite = ldiv(num, denom);

	return parite.quot + parite.rem;
}

/*Calculate the connecteed components of the graph given*/
/*by the edge list Te.*/
b_Node **connected_components(simpleEdge *Te, long n_Te, long m_Te, long *n_comp, int edge_flag, double cap_value)
{
	long non_explorer, i, u;
	Fifo *fifo;
	Bool *marquer, *explorer;
	long node_cour;
	Edge *e;
	Graph gr;
	simpleEdge *T_p; /*Table of edges with new labels*/
	simpleNode *T1, *T2;
	long n1, n2, n, m;
	long v_u, v_node_cour;
	b_Node **comp_tab;

	n1 = n_Te;
	n2 = 0;

	/*We work first on the numbered graph*/
	T_p = renumber_edge_list(Te, m_Te, &T1, n1, &T2, &n2);

	if (!InitGraph_from_list(T_p, n2, m_Te, &gr))
	{
		free(T_p);
		free(T1);
		free(T2);

		return kECSP_False;
	}

	n = n2;
	m = m_Te;

	comp_tab = new b_Node *[(n / 2)];

	marquer = (Bool *)malloc(n * sizeof(Bool));
	if (marquer == NULL)
	{
		fprintf(stdout, "%s", "ERROR:Connected components.Unable to allocate table MARQUE.");
		return kECSP_False;
	}

	for (i = 0; i < n; i++)
	{
		marquer[i] = kECSP_False;
	}

	explorer = (Bool *)malloc(n * sizeof(Bool));
	if (explorer == NULL)
	{
		fprintf(stdout, "%s", "ERROR:Connected components.Unable to allocate table EXPLORER.");

		free(marquer);
		return kECSP_False;
	}

	for (i = 0; i < n; i++)
	{
		explorer[i] = kECSP_False;
	}

	*n_comp = 0;

	fifo_init(&fifo, n);

	non_explorer = 1L;

	fifo_insert(fifo, 1L);

	while (non_explorer <= n)
	{
		node_cour = fifo_get_first(fifo);
		comp_tab[*n_comp] = NIL_BN;

		while (node_cour != -1L)
		{
			i = node_cour - 1;
			explorer[i] = kECSP_True;

			if (marquer[i] == kECSP_False)
			{
				marquer[i] = kECSP_True;

				/*Ajout d'une composante connexe*/
				v_node_cour = T2[node_cour - 1];
				add_to_connected_comp(comp_tab, n_comp, v_node_cour);
			}

			/*Parcours de la liste d'incidence*/
			e = gr.Nodes[i].first_edge;

			while (e != NIL_E)
			{
				switch (edge_flag)
				{
				case ALL_EDGES:
					u = e->adjac->id;

					if (marquer[u - 1] == kECSP_False)
					{
						marquer[u - 1] = kECSP_True;
						fifo_insert(fifo, u);

						v_u = T2[u - 1];

						add_to_connected_comp(comp_tab, n_comp, v_u);
					}

					break;

				case NO_CAP_GREATER_THAN: /*Sans les aretes de capacite >= cap_value*/
					if (e->cap < cap_value)
					{
						u = e->adjac->id;

						if (marquer[u - 1] == kECSP_False)
						{
							marquer[u - 1] = kECSP_True;
							fifo_insert(fifo, u);

							v_u = T2[u - 1];

							add_to_connected_comp(comp_tab, n_comp, v_u);
						}
					}

					break;

				case NO_0_EDGES: /*Sans les aretes de capacite nulle*/
					if (e->cap != 0.0)
					{
						u = e->adjac->id;

						if (marquer[u - 1] == kECSP_False)
						{
							marquer[u - 1] = kECSP_True;
							fifo_insert(fifo, u);

							v_u = T2[u - 1];

							add_to_connected_comp(comp_tab, n_comp, v_u);
						}
					}

					break;

				default: /*Par defaut on ne prend pas les aretes e 0*/
					if (e->cap != 0.0)
					{
						u = e->adjac->id;

						if (marquer[u - 1] == kECSP_False)
						{
							marquer[u - 1] = kECSP_True;
							fifo_insert(fifo, u);

							v_u = T2[u - 1];

							add_to_connected_comp(comp_tab, n_comp, v_u);
						}
					}

					break;
				}

				e = e->next;
			}

			node_cour = fifo_get_first(fifo);
		}

		/*When the queue is empty then we have build a connected component*/
		/*We search for the last node that hasn't been explored yet.*/
		while ((non_explorer <= n) && (explorer[non_explorer - 1] == kECSP_True))
		{
			non_explorer++;
		}

		fifo_insert(fifo, non_explorer);
		(*n_comp)++;
	}

	free(marquer);
	free(explorer);
	fifo_delete_fifo(&fifo);

	/*Rerenumerotation des sommets*/
	/*for(i=0;i<*n_comp;i++)
	{
		printf("I == %ld\n",i);

		b_cour = comp_tab[i];
		while(b_cour != NIL_BN)
		{
			printf("%ld ",b_cour->id);

			n_id = T2[b_cour->id-1];
			b_cour->id = n_id;

			b_cour = b_cour->next;
		}
		puts("");
	}*/

	free(T1);
	free(T2);
	free(T_p);
	DeleteGraph(&gr);

	// return kECSP_True;
	return comp_tab;
}

/*Initialize the graph from the edges table given by TSrc*/
BOOL InitHO(simpleEdge *Te, HO_graph **GDest, long n_Nodes, long m_Edges)
{
	long n, m, m0, i, j, nod1, nod2;
	double cap;
	HO_node *nptr, *nptr1, *nptr2;
	HO_edge *eptr1, *eptr2;
	simpleEdge *TSrc;

	n = n_Nodes;
	m = m_Edges;
	m0 = 0;

	TSrc = (simpleEdge *)malloc(m * sizeof(simpleEdge));
	m0 = 0;
	for (i = 0; i < m; i++)
	{
		if (Te[i].cap > 0)
		{
			TSrc[m0].node1 = Te[i].node1;
			TSrc[m0].node2 = Te[i].node2;
			TSrc[m0].cap = Te[i].cap;

			m0++;
		}
	}

	m = m0;

	if (!alloc_graph_ho(n, m, GDest))
	{
		NO_MEM;
		return FALSE;
	}

	(*GDest)->n_nodes = n;
	(*GDest)->n_edges = m;

	for (i = n, nptr = &((*GDest)->nodes[n - 1L]); i > 0L; --nptr, --i)
	{
		nptr->id = i;
		nptr->first_edge = HO_NILE;
	}

	eptr1 = &((*GDest)->edges[0L]);
	eptr2 = &((*GDest)->edges[m]);

	for (j = 0L; j < m; j++)
	{
		nod1 = TSrc[j].node1 - 1;
		nod2 = TSrc[j].node2 - 1;
		cap = TSrc[j].cap;

		if (cap > 0) /*L'algorithme ne supporte que*/
		{			 /*le graphe contiennent des arètes nulles */

			/* put edge into internal graph representation */
			nptr1 = &((*GDest)->nodes[nod1]);
			nptr2 = &((*GDest)->nodes[nod2]);
			eptr1->adjac = nptr2;
			eptr2->adjac = nptr1;
			eptr1->cap = cap;
			eptr2->cap = cap;
			eptr1->rcap = cap;
			eptr2->rcap = cap;
			eptr1->back = eptr2;
			eptr2->back = eptr1;

			if (nptr1->first_edge == HO_NILE)
			{
				nptr1->first_edge = eptr1;
				eptr1->next = eptr1;
			}
			else
			{
				eptr1->next = nptr1->first_edge->next;
				nptr1->first_edge->next = eptr1;
			}

			if (nptr2->first_edge == HO_NILE)
			{
				nptr2->first_edge = eptr2;
				eptr2->next = eptr2;
			}
			else
			{
				eptr2->next = nptr2->first_edge->next;
				nptr2->first_edge->next = eptr2;
			}
			++eptr1;
			++eptr2;
		}
		else
			m0++;
	}

	/*(*GDest)->n_nodes = n;
	(*GDest)->n_edges = m - m0;*/

	free(TSrc);

	return TRUE;
}

/*Calculate the global Mincut of the graph gr*/
/*and returns pointer to a structure which*/
/*represents the set W and its complementary W_b.*/
/*Its also contains the value of the global mincut.*/
Mincut *global_mincut_ho(simpleEdge *TSrc, long n, long m)
{
	b_Node *n_cour;
	long n_shore1;		 /*shore1_nodes gives the number nodes which constitute a global mincut*/
	Mincut *gmcu = NULL; /*Mincut data structure for the result*/
	double mincap;

	HO_graph *gr;
	long i;

	simpleEdge *T_p; /*Table of edges with new labels*/
	simpleNode *T1, *T2;
	long n1, n2, n_id;

	n1 = n;
	n2 = 0;

	/*We work first on the numbered graph*/
	T_p = renumber_edge_list(TSrc, m, &T1, n1, &T2, &n2);

	/*Initialisation of the graph under the NOI format*/
	if (!InitHO(T_p, &gr, n2, m))
	{
		free(T_p);
		free(T1);
		free(T2);

		return gmcu;
	}

	if (!gmincut_ho(gr, &mincap, &n_shore1))
	{
		/*free(gr->nodes);
		free(gr->edges);*/
		dealloc_graph_ho(gr);

		free(T_p);
		free(T1);
		free(T2);

		return gmcu;
	}

	/*After calculating the gmincut with the NOI function we*/
	/*build the linked list of b_Node representing the set W and*/
	/*the set W_b, according to the given option.*/

	gmcu = (Mincut *)malloc(sizeof(Mincut));
	gmcu->f_w = NIL_BN;
	gmcu->f_w_b = NIL_BN;
	gmcu->mincap = 0.0;

	/*Getting W and W_b*/
	for (i = n2 - 1; i >= 0; --i)
	{
		n_cour = (b_Node *)malloc(sizeof(b_Node));

		n_id = gr->nodes[i].id - 1;
		n_cour->id = T2[n_id];

		/*Insertion of n_cour in the list gmcu->f_w or*/
		/*in the list gmcu->f_w_b.*/
		/*We should be sure that the first node of the*/
		/*set is the node that has the lower Id(usefull*/
		/*for node contraction).Note that the list W and*/
		/*W_b will be sorted at the end of the procedure.*/
		if (gr->nodes[i].shore1 == TRUE)
		{
			/*n_cour->next = gmcu->f_w;
			gmcu->f_w = n_cour;*/
			insert_b_Node(&gmcu->f_w, &n_cour);
		}
		else
		{
			/*n_cour->next = gmcu->f_w_b;
			gmcu->f_w_b = n_cour;*/

			insert_b_Node(&gmcu->f_w_b, &n_cour);
		}
	}

	gmcu->mincap = mincap;

	/*Deleting the graph*/
	dealloc_graph_ho(gr);

	free(T_p);
	free(T1);
	free(T2);

	/*free(gr->nodes);
	free(gr->edges);
	free(gr);*/

	return gmcu;
}

static long *number;

BOOL initialize(HO_graph *gr)
{
	HO_node *nptr;
	long n, id;

	/* initial working set W contains nodes 2 to n,
	   W represented as doubly linked list           */

	n = gr->n_nodes;
	gr->nodes[n - 1L].right_link = &(gr->nodes[1L]);
	gr->nodes[1L].left_link = &(gr->nodes[n - 1L]);
	for (nptr = &(gr->nodes[n - 1L]), id = n; id >= 1L; --nptr, --id)
	{
		if (nptr->first_edge == HO_NILE)
		{ // fprintf (stderr,"Error in input graph - no incident edges for node %d\n",nptr->id);
			return (FALSE);
		}
		else
			nptr->scan_ptr = nptr->first_edge;
		nptr->excess = 0.0L;
		nptr->unmarked = TRUE;
		nptr->in_S = FALSE;
		nptr->in_W = TRUE;
		nptr->stack_link = HO_NILN;
		if (id > 2L)
			nptr->left_link = nptr - 1L;
		if (id < n)
			nptr->right_link = nptr + 1L;
	}
	return (TRUE);
}

long bfs0(HO_node *t)
{
	long level, count;
	HO_node *q_rear, *q_front, *nptr;
	HO_edge *eptr;

	t->dist = 0L;
	count = 1L;
	t->unmarked = FALSE;
	q_front = t;
	q_rear = q_front;

bfs_next:
	level = q_rear->dist + 1L;
	eptr = q_rear->first_edge;
	do
	{
		if (eptr->adjac->unmarked && eptr->back->rcap > EPS)
		{
			nptr = eptr->adjac;
			nptr->unmarked = FALSE;
			nptr->dist = level;
			++number[level];
			++count;
			q_front->bfs_link = nptr;
			q_front = nptr;
		}
		eptr = eptr->next;
	} while (eptr != q_rear->first_edge);
	if (q_rear == q_front)
		goto bfs_ready;

	q_rear = q_rear->bfs_link;
	goto bfs_next;

bfs_ready:;
	return (count);
}

void bfs1(HO_node *t)
{
	HO_node *q_front, *q_rear, *nptr;
	HO_edge *eptr;
	long level;

	t->unmarked = FALSE;
	t->dist = 0L;
	q_front = t;
	q_rear = q_front;

bfs_next:
	level = q_rear->dist + 1L;
	eptr = q_rear->first_edge;
	do
	{
		if (eptr->adjac->in_W && eptr->adjac->unmarked &&
			eptr->back->rcap > EPS)
		{
			nptr = eptr->adjac;
			nptr->unmarked = FALSE;
			nptr->dist = level;
			q_front->bfs_link = nptr;
			q_front = nptr;
		}
		eptr = eptr->next;
	} while (eptr != q_rear->first_edge);
	if (q_rear == q_front)
		goto bfs_ready;

	q_rear = q_rear->bfs_link;
	goto bfs_next;

bfs_ready:
	return;
}

BOOL gmincut_ho(HO_graph *gr, double *mincap, long *n_shore)
{
	/* Determines global minimum cut in an undirected graph,
	   i.e. a cut of minimum capacity with respect to cuts
	   between all pairs of nodes.

	   References:
	   ----------
	   J. Hao/ J.B. Orlin: "A Faster Algorithm for Finding
	   the Minimum Cut in a Graph", Proc. of the 3rd Annual
	   ACM-SIAM Symposium on Discrete Algorithms, Orlando,
	   Florida, 1992

	*/
	HO_node *s_ptr, *t_ptr, *W_ptr, *w_end_ptr,
		*aptr, *nptr, *nnptr;
	HO_node **dormant, **dor_ptr, **active, **ptr;
	HO_edge *eptr;
	long n, k, card_S, max_dor, max_dist;
	BOOL found;
	double incre, cap;
	long adist, dmin = 0, i;

	if (!initialize(gr))
		return (FALSE);

	*mincap = DBL_MAX;
	n = gr->n_nodes;
	dormant = (HO_node **)malloc(n * sizeof(HO_node *));
	if (dormant == (HO_node **)0)
	{
		NO_MEM;
		return (FALSE);
	}

	active = (HO_node **)calloc(n + 1L, sizeof(HO_node *));
	/* holds stacks of active nodes arranged by distances */
	if (active == (HO_node **)0)
	{
		NO_MEM;
		return (FALSE);
	}

	number = (long *)calloc(n, sizeof(long));
	/* counts ocurrences of distances for nodes in W */
	if (number == (long *)0)
	{
		NO_MEM;
		return (FALSE);
	}

	s_ptr = gr->nodes;
	s_ptr->in_S = TRUE;
	s_ptr->in_W = FALSE;
	card_S = 1;
	t_ptr = &(gr->nodes[n - 1L]);

	/* breadth first search to get exact distances from first
	   sink, exact distances used in test of graph connectivity */

	k = bfs0(t_ptr);
	if (k < n)
	{ /* input graph not connected */
		for (nptr = &(gr->nodes[n - 1]); nptr >= gr->nodes; nptr--)
		{
			// nptr->shore1 = ! nptr->unmarked;
			if (nptr->unmarked == TRUE)
				nptr->shore1 = FALSE;
			else
				nptr->shore1 = TRUE;
		}

		*n_shore = k;
		*mincap = 0;
		return (TRUE);
	}

	number[0L] = 1L;
	W_ptr = &(gr->nodes[2L]);

	/* initialize set of dormant nodes */
	dormant[0L] = s_ptr;
	max_dor = 0L;
	s_ptr->left_link = s_ptr;
	s_ptr->right_link = s_ptr;

	/* initial preflow push from node s = gr->nodes[0] */

	max_dist = 0L; /* = max_dist of active nodes */
	eptr = s_ptr->first_edge;
	do
	{
		nptr = eptr->adjac;
		cap = eptr->rcap;
		nptr->excess += cap;
		s_ptr->excess -= cap;
		eptr->back->rcap += cap;
		eptr->rcap = 0.0L;

		if (nptr != t_ptr && nptr->excess <= cap + EPS)
		{ /* push node nptr onto stack for nptr->dist,
		 but only once in case of double dges      */
			nptr->stack_link = active[nptr->dist];
			active[nptr->dist] = nptr;
			if (nptr->dist > max_dist)
				max_dist = nptr->dist;
		}
		eptr = eptr->next;
	} while (eptr != s_ptr->first_edge);

	/* main loop */

next_cut:

	do
	{ /* get maximum distance active node */
		aptr = active[max_dist];
		while (aptr != HO_NILN)
		{ /* remove node *aptr from stack */
			active[max_dist] = aptr->stack_link;
			eptr = aptr->scan_ptr;

			/* node *aptr will not be put back onto stack again
		   in current mincut computation, either it is
		   processed until its excess becomes zero or
		   else it will go into the set of dormant nodes */

		edge_scan: /* for current active node */
			nptr = eptr->adjac;
			if (nptr->in_W && nptr->dist == aptr->dist - 1L &&
				eptr->rcap > EPS)
			{
				incre = aptr->excess;
				if (incre <= eptr->rcap)
				{ /* perform a non saturating push */
					eptr->rcap -= incre;
					eptr->back->rcap += incre;
					aptr->excess = 0.0L;
					nptr->excess += incre;
					if (nptr != t_ptr && nptr->excess <= incre + EPS)
					{
						nptr->stack_link = active[nptr->dist];
						active[nptr->dist] = nptr;
					}
					aptr->scan_ptr = eptr;
					goto node_ready;
				}
				else
				{ /* perform a saturating push */
					incre = eptr->rcap;
					eptr->back->rcap += incre;
					aptr->excess -= incre;
					nptr->excess += incre;
					eptr->rcap = 0.0L;
					if (nptr != t_ptr && nptr->excess <= incre + EPS)
					{
						nptr->stack_link = active[nptr->dist];
						active[nptr->dist] = nptr;
					}
					if (aptr->excess <= EPS)
					{
						aptr->scan_ptr = eptr;
						goto node_ready;
					}
				}
			}
			if (eptr->next == aptr->first_edge)
			{ /* all admissable edges of current active node
			 scanned, relabel or update set of dormant
				 nodes now */
				adist = aptr->dist;
				if (number[adist] == 1L)
				{
					/* dist[j] != dist[i] for all j in W-{i},
						   extend dormant set by another layer */
					dor_ptr = &(dormant[++max_dor]);
					*dor_ptr = HO_NILN;
					nptr = W_ptr;
					w_end_ptr = W_ptr->left_link;

				transfer:
					nnptr = nptr->right_link;
					if (nptr->dist >= adist)
					{
						/* remove node nptr from set W */
						nptr->left_link->right_link =
							nptr->right_link;
						nnptr->left_link = nptr->left_link;
						if (W_ptr == nptr)
							W_ptr = nnptr;
						/* W_ptr != NILN since t_ptr
				   is  contained in W       */
						nptr->in_W = FALSE;
						--number[nptr->dist];

						/* clear stack for nptr->dist */
						active[nptr->dist] = HO_NILN;
						nptr->scan_ptr = nptr->first_edge;

						/* put node nptr into linked list
				   dormant[max_dor]              */
						if (*dor_ptr == HO_NILN)
						{
							*dor_ptr = nptr;
							nptr->right_link = nptr;
							nptr->left_link = nptr;
						}
						else
						{
							nptr->right_link = *dor_ptr;
							nptr->left_link = (*dor_ptr)->left_link;
							(*dor_ptr)->left_link = nptr;
							nptr->left_link->right_link = nptr;
						}
					}
					if (nptr == w_end_ptr)
						goto node_ready;

					nptr = nnptr;
					goto transfer;

				} /* dist[j] != dist[i] for all j in W-{i} */
				else
				{
					/* check if there is an edge (u, v), u=*aptr,
						  such that v in W and rcap(u,v) > 0    */
					eptr = aptr->first_edge;
					found = FALSE;
					do
					{
						if (eptr->adjac->in_W && eptr->rcap > EPS)
						{
							found = TRUE;
							dmin = eptr->adjac->dist;
							break;
						}
						else
							eptr = eptr->next;
					} while (eptr != aptr->first_edge);
					if (found)
					{
						aptr->scan_ptr = eptr;

						/* get new distance label for *aptr */
						while (eptr->next != aptr->first_edge)
						{
							eptr = eptr->next;
							if (eptr->adjac->in_W &&
								eptr->adjac->dist < dmin &&
								eptr->rcap > EPS)
								dmin = eptr->adjac->dist;
						}
						--number[adist];
						aptr->dist = dmin + 1L;
						++number[dmin + 1];
						max_dist = dmin + 1L;
						eptr = aptr->scan_ptr;
						goto edge_scan;
					}
					else
					{ /* extend dormant set by another
						 layer containing node *aptr only,
						 remove aptr from W nodes first    */

						aptr->in_W = FALSE;
						--number[adist];
						aptr->scan_ptr = aptr->first_edge;
						aptr->left_link->right_link =
							aptr->right_link;
						aptr->right_link->left_link =
							aptr->left_link;
						if (W_ptr == aptr)
							W_ptr = aptr->right_link;
						dormant[++max_dor] = aptr;
						aptr->right_link = aptr;
						aptr->left_link = aptr;
						goto node_ready;
					}
				}
			}
			else
			{
				eptr = eptr->next;
				goto edge_scan;
			}

		node_ready:
			aptr = active[max_dist];
		} /* aptr != NILN */
		--max_dist;
	} while (max_dist > 0L);

check_min:

	if (*mincap > t_ptr->excess)
	{
		*mincap = t_ptr->excess;
		*n_shore = 0L;
		for (nptr = &(gr->nodes[n - 1L]); nptr >= gr->nodes; --nptr)
			if (nptr->in_W)
				nptr->shore1 = FALSE;
			else
			{
				nptr->shore1 = TRUE;
				++(*n_shore);
			}
	}

	/* preparations for next cut step */

	/* delete t_ptr from W */
	t_ptr->in_W = FALSE;
	--number[t_ptr->dist];
	if (t_ptr->right_link == t_ptr)
		W_ptr = HO_NILN;
	else
	{
		t_ptr->right_link->left_link = t_ptr->left_link;
		t_ptr->left_link->right_link = t_ptr->right_link;
		if (W_ptr == t_ptr)
			W_ptr = t_ptr->right_link;
	}

	/* put t_ptr into source set S and dormant[0] set */
	t_ptr->in_S = TRUE;
	++card_S;
	t_ptr->right_link = dormant[0L]->right_link;
	t_ptr->left_link = dormant[0L];
	dormant[0L]->right_link->left_link = t_ptr;
	dormant[0L]->right_link = t_ptr;

	if (card_S == n)
		goto mincut_ready;

	/* saturate all arcs from *t_ptr to nodes not in S */
	eptr = t_ptr->first_edge;
	do
	{
		nptr = eptr->adjac;
		if (!nptr->in_S && eptr->rcap > EPS)
		{
			t_ptr->excess -= eptr->rcap;
			nptr->excess += eptr->rcap;
			eptr->back->rcap += eptr->rcap;
			eptr->rcap = 0.0L;
			nptr->scan_ptr = nptr->first_edge;
		}
		eptr = eptr->next;
	} while (eptr != t_ptr->first_edge);

	if (W_ptr == HO_NILN)
	{ /* set of W nodes empty, dormant[max_dor]
		 taken as next set of W nodes
	   */
		W_ptr = dormant[max_dor--];
		nptr = W_ptr;
		do
		{
			nptr->in_W = TRUE;
			nptr = nptr->right_link;
		} while (nptr != W_ptr);
	}

	/* get node from W with minimum distance as new sink */
	dmin = W_ptr->dist;
	t_ptr = W_ptr;
	nptr = t_ptr;
	while (nptr->right_link != W_ptr)
	{
		nptr = nptr->right_link;
		if (nptr->dist < dmin)
		{
			dmin = nptr->dist;
			t_ptr = nptr;
		}
	}

	/* breadth first search to get exact distances
	   for nodes of W with respect to new sink,
	   nodes of W with positive excess will be pushed
	   onto stack of active nodes, not all nodes of W
	   are reachable in the residual graph by breadth
	   first search, however, all such nodes are put
	   into another dormant set                       */

	if (W_ptr->right_link == W_ptr)
		/* only one node left in W = new sink */
		goto check_min;

	for (i = n - 1L, ptr = &(active[n - 1L]); i >= 0L; i--, ptr--)
	{
		number[i] = 0L;
		*ptr = HO_NILN;
	}
	nptr = t_ptr;
	while (nptr->right_link != t_ptr)
	{
		nptr = nptr->right_link;
		nptr->unmarked = TRUE;
		nptr->scan_ptr = nptr->first_edge;
	}

	max_dist = 0L;
	number[0L] = 1L;

	bfs1(t_ptr);

	/*   check next set W for nodes to be transferred
	 to another dormant set and for active nodes
	 to be pushed onto stack                      */

	dor_ptr = &(dormant[max_dor + 1L]);
	*dor_ptr = HO_NILN;
	nptr = W_ptr;
	w_end_ptr = W_ptr->left_link;

check_W:
	nnptr = nptr->right_link;
	if (nptr->unmarked)
	{
		/* remove node *nptr from set W */
		nptr->in_W = FALSE;
		nptr->right_link->left_link = nptr->left_link;
		nptr->left_link->right_link = nptr->right_link;
		if (W_ptr == nptr)
			W_ptr = nnptr;

		/* put node *nptr into new set dormant[dor_max+1] */
		if (*dor_ptr == HO_NILN)
		{
			*dor_ptr = nptr;
			nptr->right_link = nptr;
			nptr->left_link = nptr;
		}
		else
		{
			nptr->right_link = (*dor_ptr)->right_link;
			nptr->left_link = (*dor_ptr);
			nptr->right_link->left_link = nptr;
			(*dor_ptr)->right_link = nptr;
		}
	}
	else if (nptr != t_ptr)
	{
		++number[nptr->dist];
		if (nptr->excess > EPS)
		{
			nptr->stack_link = active[nptr->dist];
			active[nptr->dist] = nptr;
			if (nptr->dist > max_dist)
				max_dist = nptr->dist;
		}
	}
	if (nptr == w_end_ptr)
		goto end_check;

	nptr = nnptr;
	goto check_W;

end_check:
	if (*dor_ptr != HO_NILN)
		++max_dor;

	goto next_cut;

mincut_ready:

	/***********/
	free(number);
	free(active);
	free(dormant);
	/***********/

	return (TRUE);
}

BOOL alloc_graph_ho(long n, long m, HO_graph **gr)
{
	if ((*gr = (HO_graph *)malloc(sizeof(HO_graph))) == (HO_graph *)0)
		return (FALSE);
	if (((*gr)->nodes = (HO_node *)malloc(n * sizeof(HO_node))) == (HO_node *)0)
	{
		free(*gr);
		return (FALSE);
	}
	if (((*gr)->edges = (HO_edge *)malloc(2L * m * sizeof(HO_edge))) == (HO_edge *)0)
	{
		free(*gr);
		free((*gr)->nodes);
		return (FALSE);
	}
	return (TRUE);
}

void dealloc_graph_ho(HO_graph *gr)
{
	free(gr->nodes);
	free(gr->edges);
	free(gr);
}

/*Initialize the graph from the edges table given by TSrc*/
BOOL InitGusfield(simpleEdge *TSrc, GUS_graph **GDest, long n_Nodes, long m_Edges)
{
	long i;
	long n, m, m0, nod1, nod2;
	double cap;
	GUS_node *nptr, *nptr1, *nptr2;
	GUS_edge *eptr1, *eptr2;

	n = n_Nodes;
	m = m_Edges;

	if (!alloc_graph(n, m, GDest))
	{
		return FALSE;
	}

	(*GDest)->n_nodes = n;
	(*GDest)->n_edges0 = m;
	m0 = 0;

	for (i = n, nptr = &((*GDest)->nodes[n - 1]); i > 0L; --nptr, --i)
	{
		nptr->id = i;
		nptr->first_edge = GUS_NILE;
	}

	eptr1 = &((*GDest)->edges[0L]);
	eptr2 = &((*GDest)->edges[m]);

	for (i = 0; i < m; i++)
	{
		nod1 = TSrc[i].node1 - 1;
		nod2 = TSrc[i].node2 - 1;
		cap = TSrc[i].cap;

		if (cap > EPS)
		{ /* put edge into internal graph representation */

			nptr1 = &((*GDest)->nodes[nod1]);
			nptr2 = &((*GDest)->nodes[nod2]);
			eptr1->adjac = nptr2;
			eptr2->adjac = nptr1;
			eptr1->cap = cap;
			eptr2->cap = cap;
			eptr1->back = eptr2;
			eptr2->back = eptr1;

			if (nptr1->first_edge == GUS_NILE)
			{
				nptr1->first_edge = eptr1;
				eptr1->next = GUS_NILE;
			}
			else
			{
				eptr1->next = nptr1->first_edge;
				nptr1->first_edge = eptr1;
			}

			if (nptr2->first_edge == GUS_NILE)
			{
				nptr2->first_edge = eptr2;
				eptr2->next = GUS_NILE;
			}
			else
			{
				eptr2->next = nptr2->first_edge;
				nptr2->first_edge = eptr2;
			}

			++eptr1;
			++eptr2;
		}
		else
		{
			/* zero capacity edge not put into edge lists
					of its incident nodes, just counted        */
			m0++;
		}
	}

	(*GDest)->n_edges = m - m0;

	return TRUE;
}

/*Calculate the Gomory-Hu cut tree for the graph g and */
/*return the liste of the edges of the tree*/
simpleEdge *G(simpleEdge *TSrc, long n, long m, long *n_ghct)
{
	simpleEdge *Tab;
	long i;
	GUS_graph *g;
	simpleEdge *T_p; /*Table of edges with new labels*/
	simpleNode *T1, *T2;
	long n1, n2, n_id;

	*n_ghct = 0;

	n1 = n;
	n2 = 0;

	/*We work first on the numbered graph*/
	T_p = renumber_edge_list(TSrc, m, &T1, n1, &T2, &n2);

	/*Initialisation of Gusfield Graph Data Structure*/
	if (!InitGusfield(T_p, &g, n2, m))
	{
		dealloc_graph(g);
		return NULL;
	}

	if (ghc_tree(g))
	{
		Tab = (simpleEdge *)malloc((g->n_nodes - 1) * sizeof(simpleEdge));

		for (i = g->n_nodes; i >= 2; i--)
		{
			Tab[i - 2].node1 = g->nodes[i - 1].id;
			Tab[i - 2].node2 = g->nodes[i - 1].parent->id;
			Tab[i - 2].cap = g->nodes[i - 1].mincap;
		}
	}
	else
	{
		free(T_p);
		free(T1);
		free(T2);

		/*dealloc_graph(g);
		free(g->nodes);
		free(g->edges);*/

		return NULL;
	}

	/*Return to the old numerotation*/
	for (i = 0; i < n2 - 1; i++)
	{
		n_id = Tab[i].node1 - 1;
		Tab[i].node1 = T2[n_id];

		n_id = Tab[i].node2 - 1;
		Tab[i].node2 = T2[n_id];
	}

	*n_ghct = n2;

	free(T_p);
	free(T1);
	free(T2);

	/*dealloc_graph(g);
	free(g->nodes);
	free(g->edges);*/

	return Tab;
}

/*Calculate the GH cut tree and return the tree under*/
/*the Tree data structure representation(cf simpleEdge.h).*/
/*min and max are the lower and upper bounds in which*/
/*the capacity of a GHCT edge must be in to be a good edge.*/
Tree *GHCutTree(simpleEdge *TSrc, long n, long m, double min, double max)
{
	Tree *ghct = NULL;
	long i, node1, node2;
	double cap;
	GUS_graph *g;
	b_Edge *b_cour;	 /*Pointer to the current b_edge founded in the GHCT*/
	simpleEdge *T_p; /*Table of edges with new labels*/
	simpleNode *T1, *T2;
	long n1, n2, n_id;

	n1 = n;
	n2 = 0;

	/*We work first on the numbered graph*/
	T_p = renumber_edge_list(TSrc, m, &T1, n1, &T2, &n2);

	/*Initialisation of Gusfield Graph Data Structure*/
	if (!InitGusfield(T_p, &g, n2, m))
	{
		dealloc_graph(g);

		free(T_p);
		free(T1);
		free(T2);

		return ghct;
	}

	/*Calculation of the GH cut tree and creation the tree data structure*/
	if (ghc_tree(g))
	{
		ghct = (Tree *)malloc(sizeof(Tree));

		ghct->n_nodes = n2;

		ghct->Tab = (Tree_Node *)malloc(n2 * sizeof(Tree_Node));

		ghct->b_List = NIL_BE;

		for (i = 0; i < n2; i++)
		{
			ghct->Tab[i].id = g->nodes[i].id;
			/*ghct->Tab[i].id = i+1;*/
			ghct->Tab[i].lv = NULL;
			ghct->Tab[i].lh = NULL;
			ghct->Tab[i].cap = 0;
		}

		for (i = 1; i < n2; i++)
		{
			node1 = g->nodes[i].parent->id;
			node2 = g->nodes[i].id;
			cap = g->nodes[i].mincap;

			--node1;
			--node2;

			/*Insertion in the tree*/
			if (ghct->Tab[node1].lv == NULL)
			{
				ghct->Tab[node1].lv = &(ghct->Tab[node2]);
				ghct->Tab[node2].cap = cap;
			}
			else
			{
				ghct->Tab[node2].lh = ghct->Tab[node1].lv;
				ghct->Tab[node1].lv = &(ghct->Tab[node2]);
				ghct->Tab[node2].cap = cap;
			}

			/*Checking if the edge satisfies the criterion or not*/
			/*We keep an edge if its capacity is in the bounds.*/

			if ((min <= cap) && (cap <= max))
			{
				/*insertion in the head of the list*/
				b_cour = (b_Edge *)malloc(sizeof(b_Edge));

				b_cour->node1 = node1 + 1;
				b_cour->node2 = node2 + 1;
				b_cour->cap = cap;

				b_cour->next = ghct->b_List;
				ghct->b_List = b_cour;
			}
		} 
		/*Return to the old numerotation*/
		for (i = 0; i < n2; i++)
		{
			n_id = ghct->Tab[i].id - 1;
			ghct->Tab[i].id = T2[n_id];
		}

		b_cour = ghct->b_List;
		while (b_cour != NIL_BE)
		{
			n_id = b_cour->node1 - 1;
			b_cour->node1 = T2[n_id];
			n_id = b_cour->node2 - 1;
			b_cour->node2 = T2[n_id];

			b_cour = b_cour->next;
		} 
	}

	free(T_p);
	free(T1);
	free(T2);

	dealloc_graph(g); 

	return ghct;
}

/*Compute the st-maxflow with the Goldberg-Tarjan*/
/*algorithm implemented by Gusfield.The graph is*/
/*the simpleEdge list Te. m is the length of Te.*/
Mincut *st_maxflow(simpleEdge *Te, long n, long m, long s, long t)
{
	GUS_graph *gr;
	Mincut *gmcu;
	simpleEdge *Te_p;
	simpleNode *T1, *T2;
	long i, n1, n2, new_s, new_t;
	b_Node **b_cour_w, **b_cour_w_b;
	double maxfl;

	/*Relabelling of the edge list*/
	n1 = n;
	n2 = 0;
	Te_p = renumber_edge_list(Te, m, &T1, n1, &T2, &n2);

	/*Initializing the Gusfield Graph Data Structure*/
	InitGusfield(Te_p, &gr, n2, m);

	/*Computation of the st maxflow with the Gusfield implementation*/
	new_s = T1[s - 1];
	new_t = T1[t - 1];

	init_maxflow(n2);
	maxfl = maxflow(gr, &(gr->nodes[new_s - 1]), &(gr->nodes[new_t - 1]));
	fini_maxflow();

	maxfl = (&(gr->nodes[new_t - 1]))->excess - 1.0L;

	/*printf("Mincap == %f\n",maxfl);*/

	/*After calculation of the mincut we build the two set w and w_b*/
	if (maxfl >= (double)0)
	{
		gmcu = (Mincut *)malloc(sizeof(Mincut));

		gmcu->mincap = maxfl;

		b_cour_w = &gmcu->f_w;
		b_cour_w_b = &gmcu->f_w_b;

		for (i = 0; i < n2; i++)
		{
			if (gr->nodes[i].alive == FALSE) /*Not alive means that the node is in the s side*/
			{
				(*b_cour_w) = (b_Node *)malloc(sizeof(b_Node));

				(*b_cour_w)->id = T2[i]; /*We go back to the old numerotation*/
				(*b_cour_w)->next = NIL_BN;

				b_cour_w = &((*b_cour_w)->next);
			}
			else
			{
				(*b_cour_w_b) = (b_Node *)malloc(sizeof(b_Node));

				(*b_cour_w_b)->id = T2[i]; /*We go back to the old numerotation*/
				(*b_cour_w_b)->next = NIL_BN;

				b_cour_w_b = &((*b_cour_w_b)->next);
			}
		}

		free(T1);
		free(T2);
		free(Te_p);
		dealloc_graph(gr);

		return gmcu;
	}
	else
	{
		free(T1);
		free(T2);
		free(Te_p);
		dealloc_graph(gr);

		return NULL;
	}
}

#define NO_MEM fprintf(stderr, "Unable to allocate memory\n")

/* protoype definitions */

void print_time(float t)
{
	/*long t_sec, t_hun;
	float t1;

	if (t < 1.0)
	  printf ("0'%d [s'10ms]\n", (long) (t * 100.0));
	else
	 { t1 = t * 100.0;
	   t_sec = (long) t;
	   t_hun = (long) t1 - (long) 100 * t_sec;
	   printf ("%d'%d [s'10ms]\n", t_sec, t_hun);
	 }*/
}

BOOL alloc_graph(long n, long m, GUS_graph **gr)
{
	if ((*gr = (GUS_graph *)malloc(sizeof(GUS_graph))) == (GUS_graph *)0)
		return (FALSE);
	if (((*gr)->nodes = (GUS_node *)malloc(n * sizeof(GUS_node))) == (GUS_node *)0)
	{
		free(*gr);
		return (FALSE);
	}
	if (((*gr)->edges = (GUS_edge *)malloc(2L * m * sizeof(GUS_edge))) == (GUS_edge *)0)
	{
		free(*gr);
		free((*gr)->nodes);
		return (FALSE);
	}
	return (TRUE);
}

void dealloc_graph(GUS_graph *gr)
{
	free(gr->nodes);
	free(gr->edges);
	/*gr = NULL;*/
	// free (gr);
}

char skip_comment(FILE *fd)
{
	char any;
	char *c;
	static char buf[256];

	do
	{
		c = fgets(buf, 256, fd);
		fscanf(fd, "%c", &any);
	} while (any == 'c');
	return (any);
}

GUS_graph *get_graph(FILE *fd)
{
	long n, m, m0, i, j, nod1, nod2;
	double cap;
	GUS_node *nptr, *nptr1, *nptr2;
	GUS_edge *eptr1, *eptr2;
	GUS_graph *gr;
	char inp[4], any;

	/* Scans input file and constructs internal representation
	   of graph.

	   Undirected edges are represented by pairs of edge desrip-
	   tors associated to incident nodes, the incidence list of
	   a node is represented as a one way linked list of edge de-
	   scriptors each of which contains a pointer to an adjacent
	   node and a "back" pointer to the descriptor of the other
	   edge in the pair contained in the incidence list of the
	   adjacent node.
	*/

	fscanf(fd, "%c", &any);
	if (any == 'c')
		any = skip_comment(fd);
	if (any != 'p')
	{
		fprintf(stderr, "Problem line missing in input file\n");
		exit(1);
	}
	fscanf(fd, "%s", inp);
	if (inp[0] != 'g' || inp[1] != 'h' ||
		inp[2] != 'c' || inp[3] != 't')
	{
		fprintf(stderr,
				"problem in input file is not Gomory/Gu cut tree\n");
		return ((GUS_graph *)0);
	}

	fscanf(fd, "%ld", &n);
	fscanf(fd, "%ld", &m);

	if (!alloc_graph(n, m, &gr))
	{
		NO_MEM;
		return ((GUS_graph *)0);
	}

	gr->n_nodes = n;
	gr->n_edges0 = m;
	m0 = 0;

	// printf ("Number of nodes and edges: %d, %d\n", n, m);

	for (i = n, nptr = &(gr->nodes[n - 1]); i > 0L; --nptr, --i)
	{
		nptr->id = i;
		nptr->first_edge = GUS_NILE;
	}

	eptr1 = &(gr->edges[0L]);
	eptr2 = &(gr->edges[m]);
	for (j = 0L; j < m; j++)
	{
		if (fscanf(fd, "%s %ld %ld %lg", inp, &nod1, &nod2, &cap) == EOF)
		{ // fprintf (stderr,"EOF reached in input file, %d edges read\n", j);
			dealloc_graph(gr);
			return ((GUS_graph *)0);
		}
		if (inp[0] != 'e')
		{
			printf("syntax error - edge descriptor line expected\n");
			dealloc_graph(gr);
			return ((GUS_graph *)0);
		}
		if (nod1 < 1 || nod1 > n)
		{ // fprintf (stderr,"invalid node id %d in input graph\n", nod1);
			dealloc_graph(gr);
			return ((GUS_graph *)0);
		}
		if (nod2 < 1 || nod2 > n)
		{ // fprintf (stderr,"invalid node id %d in input graph\n", nod2);
			dealloc_graph(gr);
			return ((GUS_graph *)0);
		}
		if (nod1 == nod2)
		{ // fprintf (stderr,"loop on node %d in input graph\n", nod1);
			dealloc_graph(gr);
			return ((GUS_graph *)0);
		}
		if (cap > EPS)
		{ /* put edge into internal graph representation */
			--nod1;
			--nod2;
			nptr1 = &(gr->nodes[nod1]);
			nptr2 = &(gr->nodes[nod2]);
			eptr1->adjac = nptr2;
			eptr2->adjac = nptr1;
			eptr1->cap = cap;
			eptr2->cap = cap;
			eptr1->back = eptr2;
			eptr2->back = eptr1;
			if (nptr1->first_edge == GUS_NILE)
			{
				nptr1->first_edge = eptr1;
				eptr1->next = GUS_NILE;
			}
			else
			{
				eptr1->next = nptr1->first_edge;
				nptr1->first_edge = eptr1;
			}
			if (nptr2->first_edge == GUS_NILE)
			{
				nptr2->first_edge = eptr2;
				eptr2->next = GUS_NILE;
			}
			else
			{
				eptr2->next = nptr2->first_edge;
				nptr2->first_edge = eptr2;
			}
			++eptr1;
			++eptr2;
		}
		else
		{ /* zero capacity edge not put into edge lists
			 of its incident nodes, just counted        */
			m0++;
		}
		gr->n_edges = m - m0;
	}
	return (gr);
}

static GUS_node **active;
// static long *number;
static long max_dist, bound;
static BOOL co_check;

/* This maxflow version is for undirected graphs and
   computes maximum flow only, i.e. the resulting flow
   is not computed for all edges, the graph structure
   is expected to be initialized already, the function
   "init_maxflow" must be called first, then "maxflow"
   may be called any number of times, the function
   "fini_maxflow" should be called after the final
   maxflow call.                                      */

BOOL init_maxflow(long n)
{
	active = (GUS_node **)malloc((n + 1L) * sizeof(GUS_node *));
	/* holds stacks of active nodes arranged by distances */
	if (active == (GUS_node **)0)
	{
		printf("Unable to allocate memory\n");
		return (FALSE);
	}
	number = (long *)malloc((n + 1L) * sizeof(long));
	/* counts occurences of node distances in set
	   of alive nodes, i.e. nodes not contained in
	   the set of nodes disconnected from the sink */
	if (number == (long *)0)
	{
		printf("Unable to allocate memory\n");
		return (FALSE);
	}
	co_check = TRUE;
	return (TRUE);

} /* init_maxflow */

void fini_maxflow()
{
	free(active);
	free(number);

} /* fini_maxflow */

void global_relabel(GUS_graph *gr, GUS_node *tptr)
{
	/* breadth first search to get exact distance labels
	   from sink with reordering of stack of active nodes */

	GUS_node *front, *rear, *nptr, **ptr;
	GUS_edge *eptr;
	long n, level, count, i;

	n = gr->n_nodes;
	for (nptr = &(gr->nodes[n - 1L]); nptr >= gr->nodes; nptr--)
	{
		nptr->unmarked = TRUE;
		nptr->stack_link = GUS_NILN;
		nptr->scan_ptr = nptr->first_edge;
	}
	tptr->unmarked = FALSE;
	/* initialize stack of active nodes */
	for (ptr = &(active[n]); ptr >= active; ptr--)
		*ptr = GUS_NILN;
	for (i = 0L; i <= n; i++)
		number[i] = 0L;
	max_dist = 0L;
	count = 1L; /* number of alive nodes */
	front = tptr;
	rear = front;

bfs_next:
	level = rear->dist + 1L;
	eptr = rear->first_edge;
	while (eptr != GUS_NILE)
	{
		nptr = eptr->adjac;
		if (nptr->alive && nptr->unmarked && eptr->back->rcap > EPS)
		{
			nptr->unmarked = FALSE;
			nptr->dist = level;
			++count;
			++number[level];
			if (nptr->excess > EPS)
			{
				nptr->stack_link = active[level];
				active[level] = nptr;
				max_dist = level;
			}
			front->bfs_link = nptr;
			front = nptr;
		}
		eptr = eptr->next;
	}
	if (front == rear)
		goto bfs_ready;

	rear = rear->bfs_link;
	goto bfs_next;

bfs_ready:

	if (count < bound)
	{ /* identify nodes that are marked alive but have
	 not been reached by BFS and mark them as dead  */
		for (nptr = &(gr->nodes[n - 1L]); nptr >= gr->nodes; nptr--)
			if (nptr->unmarked && nptr->alive)
			{
				nptr->dist = n;
				nptr->alive = FALSE;
			}
		bound = count;
	}

} /* global_relabel */

double maxflow(GUS_graph *gr, GUS_node *s_ptr, GUS_node *t_ptr)
{
	/* Determines maximum flow and minimum cut between nodes
	   s (= *s_ptr) and t (= *t_ptr) in an undirected graph

	   References:
	   ----------
	   A. Goldberg/ E. Tarjan: "A New Approach to the
	   Maximum Flow Problem", in Proc. 18th ACM Symp.
	   on Theory of Computing, 1986.
	*/
	GUS_node *aptr, *nptr, *q_front, *q_rear;
	GUS_edge *eptr;
	long n, m, m0, level, i, n_discharge;
	double incre;
	long dmin;
	double cap;

	/* node ids range from 1 to n, node array indices
	   range from 0 to n-1                             */

	n = gr->n_nodes;
	for (nptr = &(gr->nodes[n - 1L]); nptr >= gr->nodes; nptr--)
	{
		nptr->scan_ptr = nptr->first_edge;
		if (nptr->scan_ptr == GUS_NILE)
		{
			fprintf(stderr, "isolated node in input graph\n");
			return (FALSE);
		}
		nptr->excess = 0.0L;
		nptr->stack_link = GUS_NILN;
		nptr->alive = TRUE;
		nptr->unmarked = TRUE;
	}
	m = gr->n_edges;
	m0 = gr->n_edges0;
	for (eptr = &(gr->edges[m - 1L]); eptr >= gr->edges; eptr--)
		eptr->rcap = eptr->cap;
	for (eptr = &(gr->edges[m0 + m - 1L]); eptr >= &(gr->edges[m0]); eptr--)
		eptr->rcap = eptr->cap;

	for (i = n; i >= 0L; i--)
	{
		number[i] = 0L;
		active[i] = GUS_NILN;
	}
	t_ptr->dist = 0L;

	/* breadth first search to get exact distances
	   from sink and for test of graph connectivity */

	t_ptr->unmarked = FALSE;
	q_front = t_ptr;
	q_rear = q_front;
bfs_next:
	level = q_rear->dist + 1L;
	eptr = q_rear->first_edge;
	while (eptr != GUS_NILE)
	{
		if (eptr->adjac->unmarked && eptr->back->rcap > EPS)
		{
			nptr = eptr->adjac;
			nptr->unmarked = FALSE;
			nptr->dist = level;
			++number[level];
			q_front->bfs_link = nptr;
			q_front = nptr;
		}
		eptr = eptr->next;
	}
	if (q_rear == q_front)
		goto bfs_ready;

	q_rear = q_rear->bfs_link;
	goto bfs_next;

bfs_ready:
	if (co_check)
	{
		co_check = FALSE;
		for (nptr = &(gr->nodes[n - 1]); nptr >= gr->nodes; --nptr)
			if (nptr->unmarked)
			{
				//fprintf(stderr, "Input graph not connected\n");
				return (-1.0L);
			}
	}

	s_ptr->dist = n; /* number[0] and number[n] not required */
	t_ptr->dist = 0L;
	t_ptr->excess = 1.0L; /* to be subtracted again */

	/* initial preflow push from source node */

	max_dist = 0L; /* = max_dist of active nodes */
	eptr = s_ptr->first_edge;
	while (eptr != GUS_NILE)
	{
		nptr = eptr->adjac;
		cap = eptr->rcap;
		nptr->excess += cap;
		s_ptr->excess -= cap;
		eptr->back->rcap += cap;
		eptr->rcap = 0.0L;

		if (nptr != t_ptr && nptr->excess <= cap + EPS)
		{ /* push node nptr onto stack for nptr->dist,
		 but only once in case of double edges     */
			nptr->stack_link = active[nptr->dist];
			active[nptr->dist] = nptr;
			if (nptr->dist > max_dist)
				max_dist = nptr->dist;
		}
		eptr = eptr->next;
	}

	s_ptr->alive = FALSE;
	bound = n;
	n_discharge = 0L;

	/* main loop */

	do
	{ /* get maximum distance active node */
		aptr = active[max_dist];
		while (aptr != GUS_NILN)
		{
			active[max_dist] = aptr->stack_link;
			eptr = aptr->scan_ptr;

		edge_scan: /* for current active node  */
			nptr = eptr->adjac;
			if (nptr->dist == aptr->dist - 1L &&
				eptr->rcap > EPS)
			{
				incre = aptr->excess;
				if (incre <= eptr->rcap)
				{ /* perform a non saturating push */
					eptr->rcap -= incre;
					eptr->back->rcap += incre;
					aptr->excess = 0.0L;
					nptr->excess += incre;
					if (nptr->excess <= incre + EPS)
					{ /* push nptr onto active stack */
						nptr->stack_link = active[nptr->dist];
						active[nptr->dist] = nptr;
					}
					aptr->scan_ptr = eptr;
					goto node_ready;
				}
				else
				{ /* perform a saturating push */
					incre = eptr->rcap;
					eptr->back->rcap += incre;
					aptr->excess -= incre;
					nptr->excess += incre;
					eptr->rcap = 0.0L;
					if (nptr->excess <= incre + EPS)
					{ /* push nptr onto active stack */
						nptr->stack_link = active[nptr->dist];
						active[nptr->dist] = nptr;
					}
					if (aptr->excess <= EPS)
					{
						aptr->scan_ptr = eptr;
						goto node_ready;
					}
				}
			}
			if (eptr->next == GUS_NILE)
			{ /* all incident arcs scanned, but node still
		 has positive excess, check if for all nptr
		 nptr->dist != aptr->dist                  */

				if (number[aptr->dist] == 1L)
				{ /* put all nodes v with dist[v] >= dist[a]
					 into the set of "dead" nodes since they
					 are disconnected from the sink          */

					for (nptr = &(gr->nodes[n - 1L]);
						 nptr >= gr->nodes; nptr--)
						if (nptr->alive &&
							nptr->dist > aptr->dist)
						{
							--number[nptr->dist];
							active[nptr->dist] = GUS_NILN;
							nptr->alive = FALSE;
							nptr->dist = n;
							--bound;
						}
					--number[aptr->dist];
					active[aptr->dist] = GUS_NILN;
					aptr->alive = FALSE;
					aptr->dist = n;
					--bound;
					goto node_ready;
				}
				else
				{ /* determine new label value */
					dmin = n;
					aptr->scan_ptr = GUS_NILE;
					eptr = aptr->first_edge;
					while (eptr != GUS_NILE)
					{
						if (eptr->adjac->dist < dmin && eptr->rcap > EPS)
						{
							dmin = eptr->adjac->dist;
							if (aptr->scan_ptr == GUS_NILE)
								aptr->scan_ptr = eptr;
						}
						eptr = eptr->next;
					}
					if (++dmin < bound)
					{ /* ordinary relabel operation */
						--number[aptr->dist];
						aptr->dist = dmin;
						++number[dmin];
						max_dist = dmin;
						eptr = aptr->scan_ptr;
						goto edge_scan;
					}
					else
					{
						aptr->alive = FALSE;
						--number[aptr->dist];
						aptr->dist = n;
						--bound;
						goto node_ready;
					}
				}
			}
			else
			{
				eptr = eptr->next;
				goto edge_scan;
			}

		node_ready:
			++n_discharge;
			if (n_discharge == n)
			{
				n_discharge = 0L;
				global_relabel(gr, t_ptr);
			}
			aptr = active[max_dist];
		} /* aptr != NILN */
		--max_dist;
	} while (max_dist > 0L);

	return (t_ptr->excess - 1.0L);
}

BOOL ghc_tree(GUS_graph *gr)
{
	/* Determines Gomory/Hu cut tree for input graph with
	   capacitated edges, the tree structures is represented
	   by parent pointers which are part of the node structure,
	   the capacity of a tree edge is stored at the child node,
	   the root of the cut tree is the first node in the list
	   of graph nodes (&gr->nodes[0]). The implementation is
	   described in [1].

	   References:
	   ----------
	   1) D. Gusfield: "Very Simple Algorithms and Programs for
		  All Pairs Network Flow Analysis", Computer Science Di-
		  vision, University of California, Davis, 1987.

	   2) R.E. Gomory and T.C. Hu: "Multi-Terminal Network Flows",
	  SIAM J. Applied Math. 9 (1961), 551-570.

	 */

	GUS_node *nptr, *nptr1, *nptrn, *sptr, *tptr;
	long n, m;
	double maxfl;

	n = gr->n_nodes;
	m = gr->n_edges;

	if (!init_maxflow(n))
		return (FALSE);

	nptr1 = gr->nodes;
	nptrn = &(gr->nodes[n - 1L]);
	for (nptr = nptrn; nptr >= nptr1; nptr--)
		nptr->parent = nptr1;

	for (sptr = &(gr->nodes[1L]); sptr <= nptrn; sptr++)
	{
		tptr = sptr->parent;
		maxfl = maxflow(gr, sptr, tptr);

		if (maxfl < 0L)
			return (FALSE);

		sptr->mincap = maxfl;
		for (nptr = &(gr->nodes[1L]); nptr <= nptrn; nptr++)
			if (nptr != sptr &&
				!nptr->alive && nptr->parent == tptr)
				nptr->parent = sptr;
		if (!tptr->parent->alive)
		{
			sptr->parent = tptr->parent;
			tptr->parent = sptr;
			sptr->mincap = tptr->mincap;
			tptr->mincap = maxfl;
		}
	}

	fini_maxflow();

	return TRUE;
}
