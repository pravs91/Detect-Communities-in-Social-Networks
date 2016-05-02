package capstone;

/**
 * Class that represents an edge in the graph
 * @author praveens
 *
 */
public class GraphEdge implements Comparable<GraphEdge> {
	private int src;
	private int dest;
	private double betweenness;
	
	public GraphEdge(int src,int dest){
		this.src = src;
		this.dest = dest;
		this.betweenness = 0;
	}
	
	public int getSrc(){
		return src;
	}
	public int getDest(){
		return dest;
	}
	public void addBetweenness(double val){
		this.betweenness += val;
	}
	
	public double getBetweenness(){
		return betweenness;
	}
	
	public void resetBetweenness(){
		this.betweenness = 0;
	}
	
	public int compareTo(GraphEdge edge){
		if(this.betweenness > edge.betweenness){
			return 1;
		} else if (this.betweenness < edge.betweenness){
			return -1;
		} else {
			return 0;
		}
	}
}
