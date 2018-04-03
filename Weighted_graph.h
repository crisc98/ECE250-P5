/*****************************************
 * Instructions
 *  - Replace 'uwuserid' with your uWaterloo User ID
 *  - Select the current calendar term and enter the year
 *  - List students with whom you had discussions and who helped you
 *
 * uWaterloo User ID:  uwuserid @uwaterloo.ca
 * Submitted for ECE 250
 * Department of Electrical and Computer Engineering
 * University of Waterloo
 * Calender Term of Submission:  (Winter|Spring|Fall) 201N
 *
 * By submitting this file, I affirm that
 * I am the author of all modifications to
 * the provided code.
 *
 * The following is a list of uWaterloo User IDs of those students
 * I had discussions with in preparing this project:
 *    -
 *
 * The following is a list of uWaterloo User IDs of those students
 * who helped me with this project (describe their help; e.g., debugging):
 *    -
 *****************************************/

#ifndef WEIGHTED_GRAPH_H
#define WEIGHTED_GRAPH_H

#ifndef nullptr
#define nullptr 0
#endif

#include <iostream>
#include <limits>
#include "Exception.h"

// include whatever classes you want

class Weighted_graph {
	private:
		int numVerts;
		int numEdges;

		double *graph_mat;
		int *degreeArr;

		static const double INF;

	public:
		Weighted_graph( int = 50 );
		~Weighted_graph();

		int degree( int ) const;
		int edge_count() const;
		double adjacent( int, int ) const;
		double distance( int, int );
		int getIndex( int, int ) const;
		int minDistance(double*, bool*) const;

		void insert( int, int, double );

	// Friends

	friend std::ostream &operator<<( std::ostream &, Weighted_graph const & );
};

const double Weighted_graph::INF = std::numeric_limits<double>::infinity();

// Your implementation here
Weighted_graph::Weighted_graph(int n){
	if(n<=0){
		n=1;
	}
	graph_mat = new double[n*n];
	degreeArr = new int[n];
	numVerts = n;
	numEdges = 0;
	for(int i=0;i<n;i++){
		degreeArr[i] = 0;
	}
	for(int i=0;i<n*n;i++){
		graph_mat[i] = INF;
	}

}

Weighted_graph::~Weighted_graph(){
	delete[] graph_mat;
	delete[] degreeArr;
}

int Weighted_graph::degree(int n) const {
	if (n>= numVerts||n<0){
		throw illegal_argument();
	}
	return degreeArr[n];
}

int Weighted_graph::edge_count() const {
	return numEdges;
}

int Weighted_graph::getIndex( int n, int m) const{
	return n*numVerts+m;
}

double Weighted_graph::adjacent(int m, int n) const {
	if(m >= numVerts || n>= numVerts|| n<0 ||m <0){
		throw illegal_argument();
	}
	if(graph_mat[getIndex(m,n)] == INF && graph_mat[getIndex(n,m)] == INF){
		return 0.0;
	}else {
		return graph_mat[getIndex(m, n)];
	}

}
int Weighted_graph::minDistance(double dist[], bool visited[]) const
{
   //Initialize min value and min index
   double min = INF;
	 int min_pos;

	 //Loop through all vertices to find the one with the min distance
	 //Update the min if a smaller vert distance is found
   for (int v = 0; v < numVerts; v++)
     if (visited[v] == 0 && dist[v] <= min)
         min = dist[v], min_pos = v;

   return min_pos;
}

double Weighted_graph::distance(int m, int n) {
	//If the distance requested is between a node and itself, return 0
	if(m==n){
		return 0.0;
	}

	//If the nodes being acted on do not correspond to realistic nodes on the graph, throw an illegal_argument Exception
	if(m >= numVerts || n >= numVerts || n < 0 || m < 0){
		throw illegal_argument();
	}

	//Initialize the distance and visited arrays.
	//Dist stores the distances from the source node to the node index of the arrays
	//visited stores whether or not the node at the index has been visited in the algo yet
	double dist[numVerts];
	bool visited[numVerts];

	//Initialize all indexes of dist and visited to infinity and false respectively
	for (int i = 0; i < numVerts; i++)
        dist[i] = INF, visited[i] = false;

	//Define the distance from to current node to itself as 0
	dist[m] = 0.0;

	for (int i = 0; i < numVerts; i++)
     {
       //Select the minimum vertex from the vertices that have not yet been visited
       int u = minDistance(dist, visited);

			 //Indicate that the vertex u has been visited
       visited[u] = true;

			 //If the selected vertex is out target vertex, n, return its distance and the algo is finished
			 if(u == n)
			 	return dist[u];

       //Update distance values for the verticies adjacent to u
       for (int v = 0; v < numVerts; v++){

				 //Update vertex distance if the vertex is connected to u, has not yet been visited,
				 // and has a path weight smaller than the current distance
				 //TLDR; if there is a better distance, use it
         if (!visited[v] && graph_mat[getIndex(u,v)] && dist[u] != INF && dist[u]+graph_mat[getIndex(u,v)] < dist[v])
            dist[v] = dist[u] + graph_mat[getIndex(u,v)];
				}
     }
		 //Return the distance from m to n
	 return dist[n];
}

void Weighted_graph::insert(int m, int n, double w){

	//Throw an illegal_argument Exception if erroneous input is given to the function
	if(w <0 ||n == m || m>=numVerts || n >= numVerts ||n < 0 || m <0){
		throw illegal_argument();
	}

	//If the the distance between the vertices is infinity,
	//increment the number of edges and the degrees of each vertex
	if(graph_mat[getIndex(m,n)] == INF){
		numEdges++;
		degreeArr[m] += 1;
		degreeArr[n] += 1;
	}

	//Update the graph matrix the argument weight between n and m and vice-versa
	graph_mat[getIndex(m,n)] = w;
	graph_mat[getIndex(n,m)] = w;
}
// You can modify this function however you want:  it will not be tested

std::ostream &operator<<( std::ostream &out, Weighted_graph const &graph ) {
	return out;
}

// Is an error showing up in ece250.h or elsewhere?
// Did you forget a closing '}' ?

#endif
