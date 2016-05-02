package capstone;

import java.util.HashMap;
import java.util.HashSet;
/**
 * 
 * @author praveens
 * Class that represents a node in the graph
 */
public class GraphNode {
	private int id;
	private HashSet<Integer> edges;
	private HashMap<Integer, GraphEdge> edgeObjects; // map between dest and its GraphEdge
	
	public GraphNode(int id){
		this.id = id;
		this.edges = new HashSet<Integer>();
		this.edgeObjects = new HashMap<Integer, GraphEdge>();
	}
	
	public void addDirectedEdge(Integer vertex){
		// HashSet will not add duplicate edges
		this.edges.add(vertex);
	}
	
	public HashSet<Integer> getEdges(){
		return new HashSet<Integer>(edges);
	}
	
	public int getID(){
		return this.id;
	}	
	
	public void addEdgeObject(int dest, GraphEdge e){
		this.edgeObjects.put(dest, e);
	}
	
	public HashMap<Integer, GraphEdge> getEdgeObjects(){
		return new HashMap<Integer, GraphEdge>(edgeObjects);
	}
	
	public void removeEdgeTo(int v){
		edges.remove(v);
		edgeObjects.remove(v);
	}
}
