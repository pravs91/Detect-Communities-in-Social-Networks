# Detect-Communities-in-Social-Networks
Detected communities in a social network graph based on Girvan-Newmann's betweenness centrality algorithm. A community is defined as 
a sub-graph in which the nodes have more edges within the sub-graph than edges linking to the rest of the graph.

# Data set
* Facebook data from the [SNAP project](https://snap.stanford.edu/data/egonets-Facebook.html)
* The input was fed from a text file. The file consists of lines with 2 integers each, corresponding 
to an edge from "from" vertex to a "to" vertex.
* The graph has 4039 nodes and 88234 edges

# Betweenness
The betweenness of an edge is a measure of the total flow it carries in the graph. If a large number of shortest paths between
any two vertices ```u``` and ```v``` go through an edge ```E```, then edge ```E``` will have a high betweenness value.

# Algorithm to find N communities
```
while numCommunitiesDetected < N: 
  do BFS for all vertices in the graph to get number of shortest paths between all vertex pairs
  Distribute flow and compute betweenness of all edges
  Sort the edges by betweenness
  Remove edge with highest betweenness
  set numCommunitiesDetected as number of Strongly-Connected-Components in graph
  reset the betweenness of all edges to zero 
  ```
# Description of classes
* ```SocialGraph``` - Interface that has the methods exposed to the user
* ```GraphNode``` - represents a vertex in the graph
* ```GraphEdge``` - represents an edge in the graph
* ```GraphLoader``` - used to load a graph from a text file
* ```CommunityDetection``` - implements the ```SocialGraph``` interface and the community detection algorithm. 
