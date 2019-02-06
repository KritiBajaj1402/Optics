#include <bits/stdc++.h>
#include "nanoflann.hpp"
#include "KDTreeVectorOfVectorsAdaptor.h"

using namespace std;

#define endl "\n"

// const int eps_for_precision = 1e-7;	// for comparing with floats/double
const double eps_for_precision = 1e-10;	// for comparing with floats/double
const double undefined = 1e9;

vector< vector<double> > points;
vector<double> reachabilityDistance;
vector<bool> isProcessed;
vector<int> orderedList;

double epsilon;
int minPts;
int numPts;
int dims;

vector<string> split(string const &input) { 
    istringstream buffer(input);
    vector<string> ret((istream_iterator<string>(buffer)), istream_iterator<string>());
    return ret;
}

void readFile(string filename = "data.tsv"){
	// probably a faster way of reading from a file
/*
	std::vector<char> v;
	if (FILE *fp = fopen("filename", "r"))
	{
		char buf[1024];
		while (size_t len = fread(buf, 1, sizeof(buf), fp))
			v.insert(v.end(), buf, buf + len);
		fclose(fp);
	}
*/

	ifstream input;
	input.open(filename);
	string pointLine;
	vector<string> splittedLine;
	// vector<int> point;
	numPts = 0;
	// dims;
	while(getline(input, pointLine)) {
		// point.clear();
		vector<double> point;
		splittedLine = split(pointLine);
		for (string coord : splittedLine){
			double point_coord = stod(coord);
			point.push_back(point_coord);
		}
		if (numPts == 0) dims = point.size();
		points.push_back(point);
		reachabilityDistance.push_back(undefined);
		isProcessed.push_back(0);
		numPts++;
	}
	cout << "read " << numPts << " points, each of dimension " << dims << endl;
}

double getSquaredDistance(int pointId1, int pointId2){
	double squaredDistance = 0;
	for (int i = 0; i < dims; ++i){
		double diff = (points[pointId1][i] - points[pointId2][i]);
		squaredDistance += diff*diff;
	}
	return squaredDistance;
}

// returns the ids of the points which are in neighborhood of the point whose id is pointId
vector< pair<double, int> > getNeighbors(int pointId){
	double squaredEpsilon = epsilon*epsilon;
	vector< pair<double, int> > neighbors;
	for (int i = 0; i < numPts; ++i){
		if (i != pointId){
			double squaredDistance = getSquaredDistance(pointId, i);
			// if (squaredDistance <= epsilon){
			if (squaredDistance <= squaredEpsilon){
				neighbors.push_back({sqrt(squaredDistance), i});
			}
		}
	}
	sort(neighbors.begin(), neighbors.end());
	return neighbors;
}

bool equality_for_floats(double d1, double d2){
	return ((d1 > d2 - eps_for_precision) && (d1 < d2 + eps_for_precision));
}

void updateSeeds(vector< pair<double, int> >& neighbors, int p, set< pair<double, int> >& seeds){
	double coreDistance = neighbors[minPts - 2].first;
	for (pair<double, int> p1 : neighbors){
		double distance = p1.first;
		int o = p1.second;
		if (!isProcessed[o]){
			double newReachDist = max(coreDistance, distance);
			// if (equality_for_floats(reachabilityDistance[o], undefined)){	// o is not in seeds; reachabilityDistance[o] is undefined
			if (reachabilityDistance[o] > undefined - 1){
				reachabilityDistance[o] = newReachDist;
				seeds.insert({newReachDist, o});
			}
			else{	// o in seeds, check for improvement
				if (newReachDist < reachabilityDistance[o] - eps_for_precision){
					// seeds.erase({reachabilityDistance[o], o});
					// set< pair<double, int> >::iterator it1 = seeds.lower_bound({(reachabilityDistance[0] - eps_for_precision), o});
					// set< pair<double, int> >::iterator it1 = seeds.lower_bound({(reachabilityDistance[o] - eps_for_precision), o});
					// set< pair<double, int> >::iterator it1 = seeds.lower_bound({(reachabilityDistance[o]), o});
					set< pair<double, int> >::iterator it1 = seeds.find({(reachabilityDistance[o]), o});
					if (it1 != seeds.end() ) seeds.erase(it1);
					else {
						cout << "not found----" << endl;
						for (auto p1 : seeds) cout << p1.first << ", " << p1.second << ";" << endl;
						cout << "reachabilityDistance[o] : " << reachabilityDistance[o] << " ,o: " << o << ", " << (reachabilityDistance[o] - eps_for_precision)
						<< ", newReachDist: " << newReachDist << endl;
					}
					// reachabilityDistance[0] = newReachDist;
					reachabilityDistance[o] = newReachDist;
					seeds.insert({newReachDist, o});
				}
			}
		}
	}
}

void optics(){
	// int p = 0;
	// for each unprocessed point p:
	for (int p = 0; p < numPts; ++p){
		if (!isProcessed[p]){
			vector< pair<double, int> > neighbors = getNeighbors(p);
			isProcessed[p] = 1;
			orderedList.push_back(p);
			
			if (neighbors.size() + 1 >= minPts){	// if p is a core point
				// double coreDistance = neighbors[minPts - 2];
				set< pair<double, int> > seeds;
				updateSeeds(neighbors, p, seeds);

				while(!seeds.empty()){
					int q = (seeds.begin())->second;
					seeds.erase(seeds.begin());
					vector< pair<double, int> > neighbors_q = getNeighbors(q);
					if (isProcessed[q]) cout << "gadbad : q : " << q << endl;
					isProcessed[q] = 1;
					orderedList.push_back(q);
					if (neighbors_q.size() + 1 >= minPts){	// if q is a core point
						updateSeeds(neighbors_q, q, seeds);
					}
				}
			}
		cout << "Ordered list size : " << orderedList.size() << endl;
		}
	}
}

void updateSeeds2(vector< pair<size_t, double> >& neighbors, int p, set< pair<double, int> >& seeds, double coreDistance){
	// double coreDistance = neighbors[minPts - 2].first;
	for (pair<int, double> p1 : neighbors){
		int o = p1.first;
		double distance = sqrt(p1.second);
		if (!isProcessed[o]){
			double newReachDist = max(coreDistance, distance);
			if (reachabilityDistance[o] > undefined - 1){
				reachabilityDistance[o] = newReachDist;
				seeds.insert({newReachDist, o});
			}
			else{	// o in seeds, check for improvement
				if (newReachDist < reachabilityDistance[o] - eps_for_precision){
					set< pair<double, int> >::iterator it1 = seeds.find({(reachabilityDistance[o]), o});
					if (it1 != seeds.end() ) seeds.erase(it1);
					else cout << "not found----" << endl;
					reachabilityDistance[o] = newReachDist;
					seeds.insert({newReachDist, o});
				}
			}
		}
	}
}

void optics2(){
	KDTreeVectorOfVectorsAdaptor< vector< vector<double> >, double > mat_index(dims, points, 10);
	mat_index.index->buildIndex();

	for (int p = 0; p < numPts; ++p){
		if (!isProcessed[p]){
			vector< pair<double, int> > neighbors = getNeighbors(p);
			isProcessed[p] = 1;
			orderedList.push_back(p);

//
			vector<double> pointQueryP(dims);
			for (size_t i = 0; i < dims; ++i)
				pointQueryP[i] = points[p][i];
			const size_t num_results = minPts;
			vector<size_t> ret_indexes(num_results);
			vector<double> squaredDistFromP(num_results);
			nanoflann::KNNResultSet<double> resultSet(num_results);
			resultSet.init(&ret_indexes[0], &squaredDistFromP[0]);
			mat_index.index->findNeighbors(resultSet, &pointQueryP[0], nanoflann::SearchParams(10));
//

			double coredist = sqrt(squaredDistFromP[num_results-1]);
			double squaredEpsilon = epsilon*epsilon;
			if (coredist < epsilon) {
			// if (neighbors.size() + 1 >= minPts){	// if p is a core point
				// double coreDistance = neighbors[minPts - 2];
				set< pair<double, int> > seeds;

				// radius search
				const double search_radius = squaredEpsilon;
				vector<pair<size_t,double> > ret_matches;
				nanoflann::SearchParams params;
				const size_t nMatches = mat_index.index->radiusSearch(&pointQueryP[0], search_radius, ret_matches, params);

				updateSeeds2(ret_matches, p, seeds, coredist);

				while(!seeds.empty()){
					int q = (seeds.begin())->second;
					seeds.erase(seeds.begin());
					
					if (isProcessed[q]) cout << "gadbad : q : " << q << endl;
					isProcessed[q] = 1;
					orderedList.push_back(q);
//
					vector<double> pointQueryQ(dims);
					for (size_t i = 0; i < dims; ++i)
						pointQueryQ[i] = points[q][i];

					// do a knn search
					vector<size_t> ret_indexes2(num_results);
					vector<double> squaredDistFromQ(num_results);
					nanoflann::KNNResultSet<double> resultSet2(num_results);
					resultSet2.init(&ret_indexes2[0], &squaredDistFromQ[0]);
					mat_index.index->findNeighbors(resultSet2, &pointQueryQ[0], nanoflann::SearchParams(10));
					//

					double coredistq = sqrt(squaredDistFromQ[num_results-1]);
					if (coredistq < epsilon) {
						// radius search
						vector<pair<size_t,double> > ret_matches2;
						nanoflann::SearchParams params2;
						const size_t nMatches2 = mat_index.index->radiusSearch(&pointQueryQ[0], search_radius, ret_matches2, params2);

						// vector< pair<double, int> > neighbors_q = getNeighbors(q);
						// if (neighbors_q.size() + 1 >= minPts){	// if q is a core point
						updateSeeds2(ret_matches2, q, seeds, coredistq);
						// }
					}
				}
			}
		cout << "Ordered list size : " << orderedList.size() << endl;
		}
	}
}

void printReachabilityDist(){
	ofstream outf;
	outf.open("2014CS10220_reachabilityDist.txt");
	for (int pointId : orderedList){
		outf << reachabilityDistance[pointId] << "\n";
	}
	outf.close();

}

int main(int argc, char* argv[]){
	cout.precision(11);
	if(argc != 3){
		cerr << "Enter 2 inputs: <minPts> <epsilon>\n";
		return 0;
	}
	minPts = stoi(argv[1]);
	epsilon = stod(argv[2]);
	cout << "minPts: " << minPts << ", epsilon: " << epsilon << endl;
	readFile();

	optics2();

	printReachabilityDist();
	
	// cout << "duplicates: "<< endl;
	// sort(orderedList.begin(), orderedList.end());
	// int len = orderedList.size();
	// for(int i = 1; i < orderedList.size(); ++i){
	// 	if (orderedList[i-1] == orderedList[i]) cout << orderedList[i] << ", ";
	// }
	// cout << endl;

	return 0;
}