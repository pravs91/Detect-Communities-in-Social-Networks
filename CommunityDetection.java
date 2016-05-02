package capstone;
/**
 * 
 * @author praveens
 * Class that represents a node in the graph
 */
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Stack;

public class CommunityDetection implements SocialGraph {
	private int numVertices;
	private int numEdges;
	
	// map between each vertex and its GraphNode
	private HashMap<Integer, GraphNode> nodeMap;
	
	// list of edges in the graph (i.e. GraphEdge objects)
	private List<GraphEdge> edges;
	
	public CommunityDetection(){
		this.numVertices = 0;
		this.numEdges = 0;
		this.nodeMap = new HashMap<Integer, GraphNode>();
		this.edges = new LinkedList<GraphEdge>();
	}
	
	/**
	 * Add vertex with id 'num'
	 * @param num
	 */
	public void addVertex(int num){
		if(!nodeMap.containsKey(num)){
			GraphNode node = new GraphNode(num);
			nodeMap.put(num, node);
			numVertices++;
		}
	}
	
	/**
	 * Add edge from 'src' to 'dest' vertex
	 * @param src
	 * @param dest
	 */
	public void addEdge(int src, int dest){
		if(!nodeMap.containsKey(src) || !nodeMap.containsKey(dest)){
			throw new IllegalArgumentException("Either 'from' or 'to' not present in the graph");
		}		
		GraphEdge edge = new GraphEdge(src,dest);
		if(!edges.contains(edge)){
			edges.add(edge);
			GraphNode srcNode = nodeMap.get(src);
			srcNode.addDirectedEdge(dest);				
			srcNode.addEdgeObject(dest, edge); // map between dest and corresponding GraphEdge object
			numEdges++;			
		}
	}
	
	/**
	 * remove the edge given by GraphEdge based on src and dest
	 * @param edge
	 */
	private void removeEdge(GraphEdge edge){
		GraphNode src = nodeMap.get(edge.getSrc());
		src.removeEdgeTo(edge.getDest());
		numEdges--;
	}
	
	public int getNumVertices(){
		return numVertices;
	}
	
	public int getNumEdges(){
		return numEdges;
	}
	
	/** 
	 * @return a map between each vertex and its edges
	 */
	public HashMap<Integer, HashSet<Integer>> exportGraph() {	
		HashMap<Integer, HashSet<Integer>> graphMap = new HashMap<Integer, HashSet<Integer>>();
		for(int num: nodeMap.keySet()){
			GraphNode node = nodeMap.get(num);
			graphMap.put(num, node.getEdges());
		}
		return graphMap;
	}	

	/**
	 * Helper method to print a graph for debugging
	 */
	private void printGraph(){
		HashMap<Integer, HashSet<Integer>> graphMap = exportGraph();
		for(Integer i: graphMap.keySet()){
			System.out.println(i + "-> " + graphMap.get(i));
		}
	}
	
	/**
	 * Helper method to print betweenness of all edges for debugging
	 */
	private void printBetweenness(){
		for(GraphEdge edge: edges){
			if(edge.getBetweenness()>0)
				System.out.println(edge.getSrc() + "->" + edge.getDest() + " betweenness: " + edge.getBetweenness());
		}

	}
	
	/**
	 * Find 'n' or more communities in the graph. A community is a strongly-connected component.
	 * Algorithm:
	 * 1) Find the BFS for all vertices in the graph to get betweenness of all edges
	 * 2) Remove the edge(s) with highest betweenness
	 * 3) Repeat 1 and 2 until the graph splits into desired number of communities (SCCs)
	 * @param n
	 */
	public List<SocialGraph> findCommunities(int n){
		int numCommunities = getSCCs().size();
		List<SocialGraph> communities = new ArrayList<SocialGraph>();
		System.out.println("Initial #SCCs: " + numCommunities);
		if(numCommunities>=n){
			System.out.println("Graph already has " + n + " connected components (communities)");
			return new ArrayList<SocialGraph>();
		}
		int edgesRemoved=0;
		// Loop until desired number of communities are found
		while(numCommunities<n && !edges.isEmpty()){
			// BFS for all verticies
			doBfsAllVertices();
			// sort the edges by betweenness
			Collections.sort(edges);
//			printBetweenness();
			edgesRemoved += removeEdgesWithHighestBetweennness();
			communities = getSCCs();
			numCommunities = communities.size();
			// reset the betweenness of all edges to 0 for next iteration
			resetBetweenness();
		}
		System.out.println("Removed " + edgesRemoved + " edges to get " + n + " communities");
//		printBetweenness();
		return communities;
	}
	
	/**
	 * BFS for all vertices in the graph
	 */
	private void doBfsAllVertices(){
		for(int v: nodeMap.keySet()){
			BFS(v);
		}		
	}
	
	/**
	 * reset the betweenness of all edges to 0
	 */
	private void resetBetweenness(){
		for(GraphEdge e: edges){
			e.resetBetweenness();
		}
	}
	
	/**
	 * Remove edge(s) with highest betweenness in edges list
	 * Use Math.round() during comparison to avoid errors due to double precision accuracy
	 * @return
	 */
	private int removeEdgesWithHighestBetweennness(){
		long maxBetweenness = Math.round(edges.get(edges.size()-1).getBetweenness()); // highest betweenness is last element
//		System.out.println(maxBetweenness);
		int index = edges.size()-1;
		int edgesRemoved = 0;
		while(index>=0 && Math.round(edges.get(index).getBetweenness())==maxBetweenness){
//			System.out.println("removing edge " + edges.get(index).getSrc() + "->" + edges.get(index).getDest()
//					+ " with betweenness " + edges.get(index).getBetweenness());
			removeEdge(edges.get(index));
			edges.remove(index);
			index--;
			edgesRemoved++;
		}		
		return edgesRemoved;
	}
	
	/**
	 * BFS from node 'start' to calculate shortest paths to all other nodes
	 * @param start
	 * 
	 */
	private void BFS(int start){
		if(!nodeMap.containsKey(start)){
			throw new IllegalArgumentException("node not present in graph");
		}
		Queue<Integer> queue = new LinkedList<Integer>();
		HashSet<Integer> visited = new HashSet<Integer>(); // visited nodes set
		// map between each node and the number of shortests paths to it from 'start'
		HashMap<Integer, Integer> numShortest = new HashMap<Integer, Integer>();
		
		// map between each node and the shortest distance to it from 'start'
		HashMap<Integer, Integer> shortestDist = new HashMap<Integer, Integer>();
		
		// map between each node and its list of parents in its shortest paths
		HashMap<Integer, List<Integer>> parentMap = new HashMap<Integer, List<Integer>>();
		
		// stack to store the order in which nodes are visited
		Stack<Integer> visitOrder = new Stack<Integer>();
		
		queue.add(start);
		visited.add(start);
		numShortest.put(start, 1);
		shortestDist.put(start, 0);
		while(!queue.isEmpty()){
			GraphNode curr = nodeMap.get(queue.remove());
			int currID = curr.getID();
			visitOrder.add(currID);
			for(int neighbor: curr.getEdges()){
				// visiting neighbor the first time
				if(!visited.contains(neighbor)){ 
					// currently, the #paths to neighbor is same as curr (first time visit)
					numShortest.put(neighbor, numShortest.get(currID));
					// currently, the shortest dist to neighbor is curr+1 (first time visit)
					shortestDist.put(neighbor, shortestDist.get(currID)+1);
					queue.add(neighbor);
					visited.add(neighbor);
					List<Integer> parents = new ArrayList<Integer>();
					parents.add(currID);// set curr as neighbor's parent
					parentMap.put(neighbor, parents);

				} else {
					// distance to neighbor through curr is same as shortest distance
					if(shortestDist.get(neighbor) == shortestDist.get(currID)+1){
						// add # shortest paths through curr to neighbor
						numShortest.put(neighbor, numShortest.get(neighbor) + numShortest.get(currID));
						parentMap.get(neighbor).add(currID); // additional parent to neighbor
					} 
					// distance through curr is shorter
					else if(shortestDist.get(neighbor) > shortestDist.get(currID) + 1){
						// update #paths and distance based on currID
						shortestDist.put(neighbor,shortestDist.get(currID)+1);
						numShortest.put(neighbor, numShortest.get(currID));
						parentMap.get(neighbor).clear();
						parentMap.get(neighbor).add(currID);
					}
					// no need to do anything if distance through curr is longer
				}
				
			}
		}
		// debug messages
//		for(Integer n: numShortest.keySet()){
//			System.out.println(n + " #paths: " + numShortest.get(n) + ";dist: " + shortestDist.get(n) 
//								+ ";has children: " + hasChildren.get(n));
//			System.out.println(n + "'s parents: " + parentMap.get(n));
//		}
//		System.out.println("visit order: " + visitOrder);
		distributeFlow(parentMap,numShortest,visitOrder);
	}
	
	/**
	 * Distribute the flow to all edges based on number of shortest paths found in BFS.
	 * The flow is distributed in bottom-up manner i.e. reverse order of visit from BFS
	 * @param parentMap - map of each node to its shortest path parents
	 * @param numShortest - number of shortest paths to each node from 'start'
	 * @param visitOrder - stack that holds visit order of BFS
	 */
	private void distributeFlow(HashMap<Integer,List<Integer>> parentMap,HashMap<Integer,Integer> numShortest, 
								Stack<Integer> visitOrder){
		HashMap<Integer, Double> outFlowMap = new HashMap<Integer, Double>();
		while(!visitOrder.empty()){
			int currID = visitOrder.pop();
//			System.out.println("visiting " + currID);
			GraphNode curr = nodeMap.get(currID);
			if(!outFlowMap.containsKey(currID)){// leaf node, outward flow is 0
				outFlowMap.put(currID, 0.0);
			}
			double inFlow = 1.0 + outFlowMap.get(currID); // inflow is 1 + outflow
			// total # shortest paths to curr node
			double numPathsToCurr = (double) numShortest.get(currID);
			if(parentMap.containsKey(currID)){
				for(int parent: parentMap.get(currID)){
					// distribute flow based on number of paths from this parent
					double betweenness = ((double) numShortest.get(parent)/numPathsToCurr) * inFlow; 
//					System.out.println("flow for " + currID + " from " + parent + ": " 
//							+ numShortest.get(parent) + "/" + numPathsToCurr + " * " + inFlow + "= " +betweenness);
					if(outFlowMap.containsKey(parent)){
						outFlowMap.put(parent,betweenness + outFlowMap.get(parent));
					} else {
						outFlowMap.put(parent,betweenness);
					}
//					System.out.println("curr: " + currID + "parent: " + parent);
					GraphEdge edgeToParent = curr.getEdgeObjects().get(parent);
					edgeToParent.addBetweenness(betweenness);
				}				
			}
		}
//		for(int n: outFlowMap.keySet()){
//			System.out.println(n + "'s out flow: " + outFlowMap.get(n));
//		}
	}
	
	/*
	 * get the strongly connected components of a directed graph
	 * Use algorithm from Coursera lecture:
	 * Step 1: Do DFS(G) keeping track of the order in which the nodes finish
	 * Step 2: Transpose G (reverse all edges)
	 * Step 3: Do DFS(G-Transpose) exploring in the reverse order of finish from step 1
	 */
	
	private List<SocialGraph> getSCCs() {
		Stack<Integer> vertices = new Stack<Integer>();
		Stack<Integer> finished = new Stack<Integer>(); // track the nodes that have finished visiting all neighbors
		HashSet<Integer> visited = new HashSet<Integer>();
		
		// create a stack of all the vertices in the graph
		for(Integer v: nodeMap.keySet()){
			vertices.push(v);
		}
		// Step 1
		DFS(vertices,finished,visited);
		
		// Step 2
		CommunityDetection transpose = getTranspose();
		
		// Step 3
		vertices = finished;
		finished = new Stack<Integer>();
		visited = new HashSet<Integer>();
		List<List<Integer>> sccList = transpose.DFS(vertices, finished, visited);
//		System.out.println(sccList);
		List<SocialGraph> SCCGraphList = new ArrayList<SocialGraph>();
		for(List<Integer> scc: sccList){
			SocialGraph sccGraph = buildSCCFromList(scc);
			SCCGraphList.add(sccGraph);
		}
		return SCCGraphList;
	}
	
	/**
	 * Do DFS on all vertices in the graph, return a list of SCCs
	 * @param vertices
	 * @param finished
	 * @param visited
	 * @return A list of SCCs in the graph
	 * Note: The list of SCCs is valid only in the second call to this function in step 3 of the algorithm
	 */
	private List<List<Integer>> DFS(Stack<Integer> vertices,Stack<Integer> finished,HashSet<Integer> visited){
		List<List<Integer>> sccList = new ArrayList<List<Integer>>();
		while(!vertices.empty()){
			int v = vertices.pop();
			if(!visited.contains(v)){
				List<Integer> scc = new ArrayList<Integer>();// contains the SCC rooted at this vertex
				DFSVisit(v,finished,visited,scc);
				sccList.add(scc);
			}
		}
		return sccList;
	}

	/**
	 * Do DFS starting at the given vertex v
	 * @param v
	 * @param finished
	 * @param visited
	 */
	private void DFSVisit(int v, Stack<Integer> finished,HashSet<Integer> visited,List<Integer> scc){
		visited.add(v);
		scc.add(v);
		GraphNode node = nodeMap.get(v);
		for(int neighbor: node.getEdges()){
			if(!visited.contains(neighbor)){
				DFSVisit(neighbor, finished, visited, scc);
			}
		}
		finished.push(v);
	}
	
	/**
	 * get the transpose of a graph i.e. all edges reversed
	 * @return transpose of graph
	 */
	private CommunityDetection getTranspose(){
		CommunityDetection transpose = new CommunityDetection();
		for(int v: nodeMap.keySet()){
			transpose.addVertex(v);
			GraphNode curr = nodeMap.get(v);
			for(int neighbor: curr.getEdges()){
				transpose.addVertex(neighbor);
				transpose.addEdge(neighbor, v);
			}
		}
		return transpose;
	}
	
	/**
	 * Takes a list of nodes in scc, and builds a sub-graph
	 * @param scc - the list containing the SCC
	 * @return a sub-graph built from scc
	 */
	private SocialGraph buildSCCFromList(List<Integer> scc){
		SocialGraph sccGraph = new CommunityDetection();
		for(int v: scc){
			sccGraph.addVertex(v);
			GraphNode curr = nodeMap.get(v);
			for(int neighbor: curr.getEdges()){
				if(scc.contains(neighbor)){
					sccGraph.addVertex(neighbor);
					sccGraph.addEdge(v, neighbor);
				}
			}
		}
		return sccGraph;
	}
	public static void main(String[] args){
		CommunityDetection graph = new CommunityDetection();
//		GraphLoader.loadGraph(graph, "data/betweenness_sample.txt");
//		GraphLoader.loadGraph(graph, "data/small_graph.txt");
//		GraphLoader.loadGraph(graph, "data/facebook_combined.txt");
//		GraphLoader.loadGraph(graph, "data/facebook_ucsd.txt");
		GraphLoader.loadGraph(graph, "data/facebook_1000.txt");

		System.out.println("Graph with " + graph.getNumVertices() + " vertices and " + graph.getNumEdges() + " edges");
		System.out.println("Finding communities...");
		int n = 5;
		long start = System.nanoTime();
		List<SocialGraph> communities = graph.findCommunities(n);
		long end = System.nanoTime();
		double time = (end-start)/Math.pow(10, 9);
		System.out.println("Split into " + n + " communities. Time: " + time);
		System.out.println("Final # communities: " + communities.size());
		for(SocialGraph scc: communities){
//			((CommunityDetection)scc).printGraph();
			System.out.println("Vertices: " + scc.getNumVertices() + " Edges: " + scc.getNumEdges());
			System.out.println("***");
		}

	}
}
