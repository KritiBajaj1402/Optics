This is an implementation of the OPTICS (Ordering points to identify the clustering structure) algorithm which was done as part of an assignment in the Data Mining course at IIT Delhi (COL 870) 

Usage:
1. Compile using : sh compile.sh
2. Run file using : sh optics.sh <minPts> <epsilon>
(reads input from file "data.tsv" and generates a file "reachabilityDist.txt" which contains the reachability distances of the points in the order in which they were processed)
3. Plot the bar graph(reachability distance) using : sh plot.sh
(uses the file "reachabilityDist.txt" generated in previous step)

KD-Tree has been used for finding the nearest neighbours and the core-distances of the points.
KD-Tree library source : https://github.com/jlblancoc/nanoflann
Functions used(from the KD-Tree library): findNeighbors() and radiusSearch()

