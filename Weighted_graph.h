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
#include <bits/stdc++.h>
#include "Exception.h"

using namespace std;
typedef pair<double, double> iPair;

// include whatever classes you want

class Weighted_graph {
	private:
		// your implementation here
		//  you can add both private member variables and private member functions
		int numVerts;
		int numEdges;

		list< pair <double, double> > *adj;

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
	numVerts = n;
	numEdges = 0;
	adj = new list<iPair> [n];
}

Weighted_graph::~Weighted_graph(){
	for(int i = 0; i < numVerts; i++){
		 adj[i].clear();
	}
	delete[] adj;
}

int Weighted_graph::degree(int n) const {
	if (n>= numVerts||n<0){
		throw illegal_argument();
	}
	return adj[n].size();
}

int Weighted_graph::edge_count() const {
	return numEdges;
}

double Weighted_graph::adjacent(int m, int n) const {
	if(m >= numVerts || n>= numVerts|| n<0 ||m <0){
		throw illegal_argument();
	}
	for (pair<double, double> p : adj[m]){
		if(p.first == n){
			return p.second;
		}
	}
	// if(graph_mat[m][n] == INF && graph_mat[n][m] == INF){
	// 	return 0.0;
	// }else {
	// 	return graph_mat[m][n];
	// } //TODO fix this for pq implementation
return 0.0;
}

double Weighted_graph::distance(int m, int n) {
	if(m==n){
		return 0.0;
	}
	if(m >= numVerts || n >= numVerts || n < 0 || m < 0 || degree(m) == 0 || degree(n) == 0){
		throw illegal_argument();
	}

	vector<double> dist(numVerts, INF);
	priority_queue< iPair, vector<iPair>, greater<iPair> > pq;
	// for(int i = 0; i <numVerts; i++){
	// 	dist[i] = INF;
	// }
	pq.push(make_pair(0,m));
	dist[m] = 0;

	while(!pq.empty()){
		int u = pq.top().second;
		pq.pop();

		list<pair<double, double> >::iterator i;
		for(i = adj[u].begin(); i!= adj[u].end(); i++){
			double v = (*i).first;
			double weight = (*i).second;

			if(dist[v] > dist[u] + weight){
				dist[v] = dist[u]+weight;
				pq.push(make_pair(dist[v],v));

			}
		}
	}
	//
	// int done = 0;
	// for(int i = 0; i <numVerts; i++){
	// 	dist[i] = INF;
	// 	visited[i]=0;
	// }
	//
	// visited[m] = 1;
	// dist[m] = 0.0;
	//  int minDistPos = m;
	//
	//  while(done == 0){
	// 	 for(int i = 0; i < numVerts; i++){
	// 		 if(graph_mat[minDistPos][i] != 0 && graph_mat[minDistPos][i] != INF && graph_mat[minDistPos][i] + dist[minDistPos]< dist[i]){
	// 			 dist[i] = graph_mat[i][minDistPos] + dist[minDistPos];
	// 		 }
	// 		 visited[minDistPos] = 1;
	// 	 }
	// 	 double minDist = INF;
	// 	 for(int i = 0; i< numVerts; i++){
	// 		 if(dist[i] <= minDist && visited[i] == 0){
	// 			 minDist = dist[i];
	// 			 minDistPos = i;
	// 		 }
	// 	 }
	// 	 done = 1;
	// 	 for(int i = 0; i< numVerts; i++){
	// 		 done = done*visited[i];
	// 	 }
	//  }
	 return dist[n];
}

void Weighted_graph::insert(int m, int n, double w){
	if(w <0 ||n == m || m>=numVerts || n >= numVerts){
		throw illegal_argument();
	}
	//if(graph_mat[m][n] == INF){
		bool found = false;
		for(pair <double, double> p : adj[m]){
			if(p.first == n){
				found =true;
				p.second = w;
			}
		}
		if (!found){
			adj[m].push_back(make_pair(n, w));
			numEdges++;
		}
		found = false;
		for(pair <double, double> p : adj[n]){
			if(p.first == m){
				found = true;
				p.second = w;
			}
		}
		if(!found){
			adj[n].push_back(make_pair(m, w));
			numEdges++;
		}
	//}
// 	graph_mat[m][n] = w;
// 	graph_mat[n][m] = w;
 }


// You can modify this function however you want:  it will not be tested

std::ostream &operator<<( std::ostream &out, Weighted_graph const &graph ) {
	return out;
}

// Is an error showing up in ece250.h or elsewhere?
// Did you forget a closing '}' ?

#endif
