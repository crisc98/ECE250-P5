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
		// your implementation here
		//  you can add both private member variables and private member functions
		int numVerts;
		int numEdges;

		double **graph_mat;
		double *dist;
		int *visited;
		int *degreeArr;

		static const double INF;

	public:
		Weighted_graph( int = 50 );
		~Weighted_graph();

		int degree( int ) const;
		int edge_count() const;
		double adjacent( int, int ) const;
		double distance( int, int );

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
	graph_mat = new double * [n];
	dist = new double[n];
	visited = new int[n];
	degreeArr = new int[n];
	numVerts = n;
	numEdges = 0;
	for(int i=0;i<n;i++){
		graph_mat[i]= new double[n];
		dist[i] = INF;
		visited[i] = 0;
		degreeArr[i] = 0;
	}
	for(int i=0;i<n;i++){
		for(int j=0; j<n;j++){
			graph_mat[i][j] = INF;
		}
	}

}

Weighted_graph::~Weighted_graph(){
	for(int i = 0; i < numVerts; i++){
		delete[] graph_mat[i];
	}
	delete[] graph_mat;
	delete[] dist;
	delete[] visited;
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

double Weighted_graph::adjacent(int m, int n) const {
	if(m >= numVerts || n>= numVerts|| n<0 ||m <0){
		throw illegal_argument();
	}
	if(graph_mat[m][n] == INF && graph_mat[n][m] == INF){
		return 0.0;
	}else {
		return graph_mat[m][n];
	}

}

double Weighted_graph::distance(int m, int n) {
	if(m==n){
		return 0.0;
	}
	if(m >= numVerts || n >= numVerts || n < 0 || m < 0 || degree(m) == 0 || degree(n) == 0){
		throw illegal_argument();
	}


	int done = 0;
	for(int i = 0; i <numVerts; i++){
		dist[i] = INF;
		visited[i]=0;
	}

	visited[m] = 1;
	dist[m] = 0.0;
	 int minDistPos = m;

	 while(done == 0){
		 for(int i = 0; i < numVerts; i++){
			 if(graph_mat[minDistPos][i] != 0 && graph_mat[minDistPos][i] != INF && graph_mat[minDistPos][i] + dist[minDistPos]< dist[i]){
				 dist[i] = graph_mat[i][minDistPos] + dist[minDistPos];
			 }
			 visited[minDistPos] = 1;
		 }
		 double minDist = INF;
		 for(int i = 0; i< numVerts; i++){
			 if(dist[i] <= minDist && visited[i] == 0){
				 minDist = dist[i];
				 minDistPos = i;
			 }
		 }
		 done = 1;
		 for(int i = 0; i< numVerts; i++){
			 done = done*visited[i];
		 }
	 }
	 return dist[n];
}

void Weighted_graph::insert(int m, int n, double w){
	if(w <0 ||n == m || m>=numVerts || n >= numVerts){
		throw illegal_argument();
	}
	if(graph_mat[m][n] == INF){
		numEdges++;
		degreeArr[m] += 1;
		degreeArr[n] += 1;
	}
	graph_mat[m][n] = w;
	graph_mat[n][m] = w;
}
// You can modify this function however you want:  it will not be tested

std::ostream &operator<<( std::ostream &out, Weighted_graph const &graph ) {
	return out;
}

// Is an error showing up in ece250.h or elsewhere?
// Did you forget a closing '}' ?

#endif
