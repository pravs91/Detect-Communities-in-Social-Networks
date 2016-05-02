package capstone;

import java.util.HashSet;
import java.util.List;
import java.util.Map;
/**
 * 
 * @author praveens
 * Class that represents a node in the graph
 */
public interface SocialGraph {
	/* Creates a vertex with the given num as ID*/
	public void addVertex(int num);
	
	/* Creates a directed edge from src to dest*/
	public void addEdge(int src,int dest);
	
	/* returns the number of vertices in the graph */
	public int getNumVertices();
	
	/* returns the number of edges in the graph */
	public int getNumEdges();
	
	/* Return the graph's connections in a readable format. 
     * The keys in this HashMap are the vertices in the graph.
     * The values are the nodes that are reachable via a directed
     * edge from the corresponding key. 
	 * The returned representation ignores edge weights and 
	 * multi-edges.  */
	public Map<Integer, HashSet<Integer>> exportGraph();
	
	/* Returns n communities from the graph, where a community is 
	 * a strongly connected component of the graph.
	 */
	public List<SocialGraph> findCommunities(int n);
}
