import gurobi.*;
import org.jgrapht.GraphPath;
import org.jgrapht.alg.shortestpath.AllDirectedPaths;
import org.jgrapht.alg.util.Pair;
import org.jgrapht.graph.DefaultDirectedWeightedGraph;

import java.io.File;
import java.util.*;

public class PathsMIP {
    GRBModel modelRMP;
    Map<Pair<String, String>, GRBVar> x;
    Map<List<String>, GRBVar> y;
    DefaultDirectedWeightedGraph<String, MyWeightedEdge> g;
    List<String> endVertices;
    List<GraphPath<String,MyWeightedEdge>> allPaths;
    HashMap<String, String> vertexNames;
    boolean ASML;
    long duration;
    public Map<Pair<String, String>, Double> solvePathsMIP(DefaultDirectedWeightedGraph<String, MyWeightedEdge> g,
                                                         List<String> endVertices, boolean ASML, HashMap<String, String> vertexNames) throws GRBException {
        long startTime = System.currentTimeMillis();
        this.g = g;
        this.endVertices = endVertices;
        //Compute all paths from the root to the end vertices
        int maxPathLength = g.vertexSet().size();
        AllDirectedPaths ADP = new AllDirectedPaths(g);
        Set<String> r = new HashSet<>(Collections.singleton("root"));
        //compute all paths; true means an edge can only be used once (so no cycles)
        this.allPaths = ADP.getAllPaths(r, new HashSet(endVertices), true, maxPathLength);
        this.vertexNames = vertexNames;
        this.ASML = ASML;

//        for (GraphPath<String, MyWeightedEdge> path : allPaths){
//            for (MyWeightedEdge e:path.getEdgeList()){
//                System.out.print("(" + g.getEdgeSource(e) + ", " + g.getEdgeTarget(e) + ") ");
//            }
//            System.out.print("\n");
//        }
        //Select a subset of variables: paths P' with arcs A' (by shortest path tree)
//        LetsGraph lg = new LetsGraph();
//        List rl = new ArrayList(Collections.singleton("root"));
//        DefaultDirectedWeightedGraph<String, MyWeightedEdge> shortestPathsForest = lg.shortestPathsForest(g, rl, endVertices);

        //Dit stond aan:
//        List<GraphPath<String, MyWeightedEdge>> initialPaths = new ArrayList<>();
//        for (String v : endVertices) {
//            initialPaths.add(DijkstraShortestPath.findPathBetween(g, "root", v));
//        }


        //Create environment:
        GRBEnv env = new GRBEnv(true);
        File logFile = new File("pmip.log");
        env.set("logFile", "pmip.log");
        env.set(GRB.IntParam.OutputFlag, 0);
        env.start();
        //Create empty model:
        this.modelRMP = new GRBModel(env);
        setupProblem();


        modelRMP.optimize();
//        System.out.println("Obj: " + modelRMP.get(GRB.DoubleAttr.ObjVal));


        //Save the results into a map
        Map<Pair<String, String>, Double> xresults = new HashMap<>();
        for(Pair<String, String> xe:x.keySet()) {
            xresults.put(xe, x.get(xe).get(GRB.DoubleAttr.X));
        }

        // Dispose of model and environment
        modelRMP.dispose();
        env.dispose();

        long endTime = System.currentTimeMillis();
        this.duration = (endTime - startTime);
//        System.out.println("Execution time solving paths MIP: " + duration + "ms");
        return xresults;
    }

    /**
     * Set up the problem
     */
    public void setupProblem() throws GRBException {
        int countvars = 0;
        int countcstrs = 0;
        //LP opstellen:
        //Create variables:
        this.x = new HashMap<>();
        for (MyWeightedEdge e:g.edgeSet()) {
            Pair<String, String> s = new Pair<>(g.getEdgeSource(e), g.getEdgeTarget(e));
            x.put(s, modelRMP.addVar(0.0, 1.0, g.getEdgeWeight(e), GRB.BINARY,
                    g.getEdgeSource(e) + "," + g.getEdgeTarget(e)));
            countvars += 1;
        }
        this.y = new HashMap<>();
        //Variable names for the paths are simply 0,1,2,3,... because other names quickly become too long
        int it = 0;
        for (GraphPath<String, MyWeightedEdge> path:allPaths) {
            String pathstring = Integer.toString(it);
            y.put(path.getVertexList(), modelRMP.addVar(0.0, 1.0, 0.0, GRB.CONTINUOUS, pathstring));
            countvars += 1;
            it += 1;
        }

        //Constraints aanmaken:
        //Each terminal has a path
        for (String t:endVertices) {
            GRBLinExpr expr = new GRBLinExpr();
            String cstr = "eachTerminalAPath(" + t + ")";
            for(GraphPath<String, MyWeightedEdge> path:allPaths) {
                if (path.getVertexList().get(path.getVertexList().size()-1).equals(t)) {
                    expr.addTerm(1.0, y.get(path.getVertexList()));
                }
            }
            modelRMP.addConstr(expr, GRB.GREATER_EQUAL, 1.0, cstr);
            countcstrs += 1;
        }
        //Path Edge Link
        int M = endVertices.size();
        for (MyWeightedEdge e:g.edgeSet()) {
            GRBLinExpr expr = new GRBLinExpr();
            String cstr = "pathEdgeLink(" + g.getEdgeSource(e) + "," + g.getEdgeTarget(e) + ")";
            expr.addTerm(M, x.get(new Pair<> (g.getEdgeSource(e), g.getEdgeTarget(e))));
            for(GraphPath<String, MyWeightedEdge> path:allPaths) {
                if(path.getEdgeList().contains(e)) {
                    expr.addTerm(-1.0, y.get(path.getVertexList()));
                }
            }
            modelRMP.addConstr(expr, GRB.GREATER_EQUAL, 0.0, cstr);
            countcstrs += 1;
        }
//        System.out.println("variables: " + countvars + ", constraints: " + countcstrs);
    }
}
