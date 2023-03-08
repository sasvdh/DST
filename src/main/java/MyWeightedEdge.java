import org.jgrapht.graph.DefaultWeightedEdge;

class MyWeightedEdge extends DefaultWeightedEdge {
    @Override
    public String toString() {
        return Double.toString(getWeight());
    }
}