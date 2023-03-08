import gurobi.*;
import org.jgrapht.Graphs;
import org.jgrapht.alg.flow.PushRelabelMFImpl;
import org.jgrapht.alg.interfaces.MinimumSTCutAlgorithm;
import org.jgrapht.alg.util.Pair;
import org.jgrapht.graph.DefaultDirectedWeightedGraph;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Class for the directed cut linear program
 */
public class DirectedCutILP extends GRBCallback {
    boolean oldMethod = false;
    private int N;
    List<List<String>> allRightSubsets;
    List<List<String>> allRightSubsetsAlternative;
    List<String> endVertices;
    DefaultDirectedWeightedGraph<String, MyWeightedEdge> g;
    DefaultDirectedWeightedGraph<String, MyWeightedEdge> resultGraphDirectedCut;
    long timeMinimizing;
    long timeDictionary;
    long duration;
    boolean allSubsets = true; //true if traditional method, false if my own alternative method that takes only the paths from root to terminal
    long startTimeSubset = 0;
    long endTimeSubset = 0;
    private Map<Pair<String, String>, GRBVar> x;
    public GRBVar[] vars ;
    private List<MyWeightedEdge> arcList;
    private Map<Pair<String, String>, Double> solutionMap;
    FileWriter myWriter;
    public DirectedCutILP(DefaultDirectedWeightedGraph<String, MyWeightedEdge> g, List<String> endVertices){
        this.g = g;
        this.endVertices = endVertices;
    }


    /**
     * Solve the LP described in https://epubs.siam.org/doi/pdf/10.1137/15M1007185
     * Given an input graph and a list of end vertices, this returns the tree containing the solution that was found.
     */
    public Map<Pair<String, String>, Double> solveDirectedCutILP(DefaultDirectedWeightedGraph<String, MyWeightedEdge> g, List<String> endVertices) throws GRBException, IOException {
        long startTime = System.currentTimeMillis();
        List<String> array = new ArrayList<>(g.vertexSet());
        array.remove("root");
        this.g = g;
        this.endVertices = endVertices;
        this.resultGraphDirectedCut = new DefaultDirectedWeightedGraph<>(MyWeightedEdge.class);

        //LP opstellen:
        //Create environment:
        GRBEnv env = new GRBEnv(true);
        env.set("logFile", "blp.log");
        env.set(GRB.IntParam.OutputFlag, 0);
        env.set(GRB.DoubleParam.TimeLimit, 900.0);
        env.start();

        //Create empty model:
        GRBModel model = new GRBModel(env);

        //Create variables:
        Map<Pair<String, String>, GRBVar> x = new HashMap<>();
        for (MyWeightedEdge e:g.edgeSet()) {
            Pair<String, String> s = new Pair<>(g.getEdgeSource(e), g.getEdgeTarget(e));
            x.put(s, model.addVar(0.0, 1.0, g.getEdgeWeight(e), GRB.BINARY,
                    g.getEdgeSource(e) + "," + g.getEdgeTarget(e)));
        }

        //If we determine al subsets beforehand and set up constraints
        if(oldMethod) {
            setupConstraints(array, model, x);

            //Optimize model
            model.optimize();

            System.out.println("Obj: " + model.get(GRB.DoubleAttr.ObjVal));

            //Save the results into a map
            Map<Pair<String, String>, Double> xresults = new HashMap<>();
            for(Pair<String, String> xe:x.keySet()) {
                xresults.put(xe, x.get(xe).get(GRB.DoubleAttr.X));
            }

            // Dispose of model and environment
            model.dispose();
            env.dispose();

            long endTime = System.currentTimeMillis();
            long duration = (endTime - startTime);
            System.out.println("Dictionary time directed cup ILP: " + timeDictionary + "ms");
            System.out.println("Minimizing time directed cup ILP: " + timeMinimizing + "ms");
            if (oldMethod) {
                long durationPre = (endTimeSubset - startTimeSubset);
                System.out.println("Preprocessing time directed cup ILP: " + durationPre + "ms");
            }
            System.out.println("Execution time solving directed cup ILP: " + duration + "ms");
            return xresults;



        }else{ //////Method with callback//////////////////////////////////////////////////////////////////////////////////////////////////////
            //Optimize model
            myWriter = new FileWriter("cutsForDirectedCut.txt"); //we save the added cuts in a .txt file
            model.set(GRB.IntParam.LazyConstraints, 1); //set the lazy constraints parameter

            //Callback stuff
            DirectedCutILP BLP = new DirectedCutILP(g, endVertices);
            arcList = new ArrayList(g.edgeSet());
            BLP.vars = new GRBVar[x.size()]; //initialize vars
            for (int i = 0; i < x.size(); i++){ //set the vars of BLP to the GRBVars
                Pair<String, String> arc = new Pair<>(g.getEdgeSource(arcList.get(i)), g.getEdgeTarget(arcList.get(i)));
                BLP.vars[i] = x.get(arc);
            }
            model.setCallback(BLP);
            model.optimize();


            //In case of problems:
//        model.computeIIS();
//        model.write("model.ilp");

            myWriter.close();
//            System.out.println("Obj: " + model.get(GRB.DoubleAttr.ObjVal));

            //Save the results into a map
            Map<Pair<String, String>, Double> xresults = new HashMap<>();
            if (model.get(GRB.IntAttr.Status) == GRB.OPTIMAL) {//If an optimal solution was found
                for (Pair<String, String> xe : x.keySet()) {
                    xresults.put(xe, x.get(xe).get(GRB.DoubleAttr.X));
                }
            }else{//If no optimal solution was found
                for (Pair<String, String> xe : x.keySet()) {
                    xresults.put(xe, 0.0);
                }
            }

            // Dispose of model and environment
            model.dispose();
            env.dispose();

            long endTime = System.currentTimeMillis();
            this.duration = (endTime - startTime);
//            System.out.println("Dictionary time directed cup ILP: " + timeDictionary + "ms");
//            System.out.println("Minimizing time directed cup ILP: " + timeMinimizing + "ms");
//            System.out.println("Execution time solving directed cup ILP: " + duration + "ms");
            return xresults;
        }

    }



    /**
     * Callback method
     */
    protected void callback() {
            try {
                if (where == GRB.CB_MIPSOL) {
                    arcList = new ArrayList<>(g.edgeSet());
                    //Save the solution into the solutionMap for arc-solution, and grbvarMap for GRBVar-solution
                    double[] sol = getSolution(vars);
                    solutionMap = new HashMap<>(); //map with arc Pair<String,String> and the solution value
                    for (int i = 0; i < sol.length; i++){
                        Pair<String, String> arc = new Pair<>(g.getEdgeSource(arcList.get(i)), g.getEdgeTarget(arcList.get(i)));
                        solutionMap.put(arc, sol[i]);
                    }
                    Map<GRBVar, Double> grbvarMap = new HashMap<>(); //map with the Gurobi variable and the solution value
                    Map<Pair<String, String>, GRBVar> grbvarArc = new HashMap<>(); //map with the arc Pair<String, String> and the gurobi variable
                    for (GRBVar y:vars){
                        List<String> verts = Arrays.asList(y.get(GRB.StringAttr.VarName).split(","));
                        Pair<String, String> arc = new Pair<>(verts.get(0), verts.get(1));
                        grbvarMap.put(y, solutionMap.get(arc));
                        grbvarArc.put(arc, y);
                    }
                    //Create a graph containing the solution
                    DefaultDirectedWeightedGraph<String, MyWeightedEdge> gcopy = new DefaultDirectedWeightedGraph<>(MyWeightedEdge.class);
                    Graphs.addAllVertices(gcopy, g.vertexSet());
                    for (MyWeightedEdge a:g.edgeSet()){
                        String source = g.getEdgeSource(a);
                        String target = g.getEdgeTarget(a);
                        MyWeightedEdge anew = gcopy.addEdge(source, target);
                        gcopy.setEdgeWeight(anew, solutionMap.get(new Pair<>(source, target)));
                    }

                    //Compute the minimum of the minimum r-x cuts for terminals x
                    MinimumSTCutAlgorithm<String, MyWeightedEdge> MinCut = new PushRelabelMFImpl<>(gcopy);
                    double minCutValue = Double.POSITIVE_INFINITY;
                    Set<String> minCutSet = new HashSet<>();
                    for(String x:endVertices) {
                        double currentMinCutValue = MinCut.calculateMinCut("root", x);
                        if (currentMinCutValue < minCutValue){ //if the minimum cut for current terminal is better
                            minCutValue = currentMinCutValue;
                            minCutSet = MinCut.getSinkPartition();
                        }
                    }

                    if(minCutValue < 0.999) {//Check if this set induces a violated constraint+ "\n");
                        //add cut to restricted master problem
                        GRBLinExpr expr = new GRBLinExpr();
                        List<MyWeightedEdge> ingoingArcs = ingoingArcs(new ArrayList<>(minCutSet));
                        String cutStr = "";
                        for (MyWeightedEdge arc : ingoingArcs) { //for each arc going into S
                            expr.addTerm(1.0, grbvarArc.get(new Pair<>(gcopy.getEdgeSource(arc), gcopy.getEdgeTarget(arc))));
                            cutStr += " + (" + gcopy.getEdgeSource(arc) + "," + gcopy.getEdgeTarget(arc) + ")";
                        }
                        cutStr += " >= 1 \n";
                        addLazy(expr, GRB.GREATER_EQUAL, 1.0);
                        try {
                            FileWriter fw = new FileWriter("cutsForDirectedCut.txt", true);
                            BufferedWriter bw = new BufferedWriter(fw);
                            bw.write(cutStr);
                            bw.close();
                        } catch (IOException e) {
                            throw new RuntimeException(e);
                        }
                    }
//                    ///////////////////////////////////////////////////////////////////////////////////////////
//                    //Check if this solution contains a path from r to all terminals by computing the shortest paths
//                    DijkstraShortestPath<String, MyWeightedEdge> dijkstra = new DijkstraShortestPath(solGraph);
//                    ShortestPathAlgorithm.SingleSourcePaths<String, MyWeightedEdge> sptreeroot =  dijkstra.getPaths("root");
//                    for (String terminal:endVertices){
//                        if(sptreeroot.getPath(terminal) == null) { //if there is no path to terminal
//                            //Create a cut for each subset that contains at least one terminal but not the root
//                            List<String> usedVertices = new ArrayList<>();
//                            for (Pair<String, String> arc:solutionMap.keySet()) {
//                                if (solutionMap.get(arc) > 0.5) { //if y=1
//                                    usedVertices.add(arc.getFirst());
//                                    usedVertices.add(arc.getSecond());
//                                }
//                            }
//                            usedVertices.remove("root");
//                    ///////////////////////////////////////////////////////////////////////////////////////////

//                            ///////////////////////////////////////////////////////////////////////////////////////////
//                            //We add a subset cut
//                            List<List<String>> rightSubsets = allRightSubsets(usedVertices);
//                            for (List<String> S : rightSubsets) { //add a cut for each subset
//                                GRBLinExpr expr = new GRBLinExpr();
//                                List<MyWeightedEdge> ingoingArcs = ingoingArcs(S);
//                                String cutStr = "";
//                                for (MyWeightedEdge arc : ingoingArcs) { //for each arc going into S
//                                    expr.addTerm(1.0, grbvarArc.get(new Pair<>(g.getEdgeSource(arc), g.getEdgeTarget(arc))));
//                                    cutStr += " + (" + g.getEdgeSource(arc) + "," + g.getEdgeTarget(arc) + ")";
//                                }
//                                cutStr += " >= 1 \n";
//                                addLazy(expr, GRB.GREATER_EQUAL, 1.0);
//                                try {
//                                    FileWriter fw = new FileWriter("cutsForDirectedCut.txt", true);
//                                    BufferedWriter bw = new BufferedWriter(fw);
//                                    bw.write(cutStr);
//                                    bw.close();
//                                }catch(IOException e){
//                                    System.out.println("Error message : " + e.getMessage());
//                                    e.printStackTrace();
//                                }
//                            }
//                            ///////////////////////////////////////////////////////////////////////////////////////////



                            ///////////////////////////////////////////////////////////////////////////////////////////
//                            //we add a combinatorial cut
//                            String cutString = "";
//                            GRBLinExpr expr = new GRBLinExpr();
//                            for (GRBVar y : vars) { //add terms
//                                if (grbvarMap.get(y) > 0.5) { //if y=1, add term (1-y)
//                                    cutString += (" + (1 - y(" + y.get(GRB.StringAttr.VarName) + ") )");
//                                    expr.addConstant(1);
//                                    expr.addTerm(-1, y);
//                                } else { //if y=0, add term y
//                                    cutString += " + y(" + y.get(GRB.StringAttr.VarName) + ")";
//                                    expr.addTerm(1, y);
//                                }
//                            }
//                            System.out.print(expr.getConstant() + " + ");
//                            for(int i = 0; i < expr.size(); i++){
//                                System.out.print(expr.getCoeff(i) + " y(" + expr.getVar(i).get(GRB.StringAttr.VarName) + ") + ");
//                            }
//                            System.out.println(">= 1");
//                            cutString += ">= 1 \n \n";
//                            try {
//                                FileWriter fw = new FileWriter("cutsForDirectedCut.txt", true);
//                                BufferedWriter bw = new BufferedWriter(fw);
//                                bw.write(cutString);
//                                bw.close();
//                            }catch(IOException e){
//                                System.out.println("Error message : " + e.getMessage());
//                                e.printStackTrace();
//                            }
//                            addLazy(expr, GRB.GREATER_EQUAL, 1);
                            ///////////////////////////////////////////////////////////////////////////////////////////
//                        }
//                    }
                }
            } catch (GRBException e) {
                System.out.println("Error code : " + e.getErrorCode() + ". " + e.getMessage());
                e.printStackTrace();
            }
    }

    /**
     * This method adds the intitial constraints for the method that uses the callback
     */
    public void addConstraintsCallbackMethod(List<String> array, GRBModel model, Map<Pair<String, String>, GRBVar> x) throws GRBException {
        for (String terminal:endVertices) { //each terminal has an incoming arc
            GRBLinExpr expr = new GRBLinExpr();
            String cstr = "ingoingEdges[";
            for (MyWeightedEdge arc : g.incomingEdgesOf(terminal)) { //for each arc going into x
                expr.addTerm(1.0, x.get(new Pair<>(g.getEdgeSource(arc), g.getEdgeTarget(arc))));
                cstr += "(" + g.getEdgeSource(arc) + "," + g.getEdgeTarget(arc) + "), ";
            }
            cstr += "]";
            model.addConstr(expr, GRB.GREATER_EQUAL, 1.0, cstr);
        }
        GRBLinExpr expr = new GRBLinExpr(); //root has an outgoing arc
        String cstr = "ingoingEdges[";
        for (MyWeightedEdge arc : g.outgoingEdgesOf("root")) { //for each arc going into x
            expr.addTerm(1.0, x.get(new Pair<>(g.getEdgeSource(arc), g.getEdgeTarget(arc))));
            cstr += "(" + g.getEdgeSource(arc) + "," + g.getEdgeTarget(arc) + "), ";
        }
        cstr += "]";
        model.addConstr(expr, GRB.GREATER_EQUAL, 1.0, cstr);
    }

    /**
     * This method determines the subsets that contain the root and at least one terminal
     * And sets up the constraints
     */
    public void setupConstraints(List<String> array, GRBModel model, Map<Pair<String, String>, GRBVar> x) throws GRBException {
        //determine all subsets that don't contain the root but contain at least one terminal
        startTimeSubset = System.currentTimeMillis();
        List<List<String>> subsets;
        if (allSubsets) {
            subsets = allRightSubsets(array); //determine all subsets that don't contain the root
        } else {
            subsets = allRightSubsetsAlternative(); //determine all subsets that don't contain the root
        }
        endTimeSubset = System.currentTimeMillis();

        //Constraints aanmaken
        for (List<String> S : subsets) {
            GRBLinExpr expr = new GRBLinExpr();
            List<MyWeightedEdge> ingoingArcs = ingoingArcs(S);
            String cstr = "ingoingEdges[";
            for (MyWeightedEdge arc : ingoingArcs) { //for each edge going into S
                expr.addTerm(1.0, x.get(new Pair<>(g.getEdgeSource(arc), g.getEdgeTarget(arc))));
                cstr += "(" + g.getEdgeSource(arc) + "," + g.getEdgeTarget(arc) + "), ";
            }
            cstr += "]";
            model.addConstr(expr, GRB.GREATER_EQUAL, 1.0, cstr);
        }
    }



    /**
     * Edges with weight 0 may be in this graph even though they do not connect a root to a terminal.
     * This method is used to eliminate those edges and vertices in a recursive way.
     */
    public DefaultDirectedWeightedGraph<String, MyWeightedEdge> removeIf0(
                            DefaultDirectedWeightedGraph<String, MyWeightedEdge> g, String v,
                            DefaultDirectedWeightedGraph<String, MyWeightedEdge> resultGraphDirectedCut){
        this.resultGraphDirectedCut = resultGraphDirectedCut;
        removeIf0sub(g, v);
        return resultGraphDirectedCut;
    }

    /**
     * This method is used in the recursion in the removeIf0 method
     */
    public void removeIf0sub(DefaultDirectedWeightedGraph<String, MyWeightedEdge> g, String v){
        if (resultGraphDirectedCut.outDegreeOf(v) == 0){ //if the vertex has no outgoing edges
            if(!endVertices.contains(v)) { //and it is not an end vertex
                Set<MyWeightedEdge> incomingEdges = new HashSet<>(resultGraphDirectedCut.incomingEdgesOf(v));
                resultGraphDirectedCut.removeVertex(v);
                if (incomingEdges.size() > 0) { //if it has incoming edges
                    String source = g.getEdgeSource(incomingEdges.iterator().next());
                    removeIf0sub(g, source);
                }
            }
        }
    }

    /**
     * For an input graph and vertex set, this returns the list of arcs going into the vertex set.
     */
    public List<MyWeightedEdge> ingoingArcs(List<String> vertexSet){
        List<MyWeightedEdge> ingoingArcs = new ArrayList<>();
        for (String v: vertexSet){
            for (MyWeightedEdge arc:g.incomingEdgesOf(v)){
                if (!vertexSet.contains(g.getEdgeSource(arc))){ //if edge does not come from vertex that is in the vertex set
                    ingoingArcs.add(arc);
                }
            }
        }
        return ingoingArcs;
    }



    /**
     * For an input array, this returns the list of all subsets that doesn't contain the root and contains
     * at least one of the end vertices (note that the root was already excluded from the array).
     * I added a way to reduce this set to the minimal sets (see commented code) but doing this took way longer
     * than the time I win by reducing the model in this way, so this is not worth it.
     */
    public List<List<String>> allRightSubsets(List<String> array){
        N = array.size();
        allRightSubsets = new ArrayList<>();
        for(String x:endVertices) {
            array.remove(x);
            N = N-1;
            subset(0, new ArrayList<>(), array, x);
        }
//        ///////////////////Minimal subsets////////////////////////////////////////////////////
//        System.out.println("Total number of susbets: " + allRightSubsets.size());
//        //Minimal subsets
//        long startTime = System.currentTimeMillis();
//        //Create dictionary for vertex set and edges going into this set
//        Map<List<String>, List<MyWeightedEdge>> ingoing = new HashMap<>();
//        for (List<String> vertexSubset : allRightSubsets){
//            ingoing.put(vertexSubset, ingoingArcs(vertexSubset));
//        }
//        long endTimeDictionary = System.currentTimeMillis();
//        timeDictionary = endTimeDictionary - startTime;
//        //Remove the vertex sets that are not minimal
//        List<String>[] allRightSubsetsArray = allRightSubsets.toArray(new List[0]);
//        int len = allRightSubsetsArray.length;
//        for (int i = 0; i< len; i++){
//            List<String> Ev = allRightSubsetsArray[i];
//            if(Ev == null) continue;
//            List<MyWeightedEdge> E = ingoing.get(Ev);
//            for (int j = 0; j< len; j++){
//                List<String> Dv = allRightSubsetsArray[j];
//                if(Dv == null) continue;
//                List<MyWeightedEdge> D = ingoing.get(Dv);
//                if (D!=E){
//                    if (new HashSet<>(D).containsAll(E)){
//                        allRightSubsetsArray[j] = null;
//                    }else if (new HashSet<>(E).containsAll(D)){
//                        allRightSubsetsArray[i] = null;
//                    }
//                }
//            }
//        }
//        long endTime = System.currentTimeMillis();
//        timeMinimizing = endTime - startTime;
//
//        var output = Stream.of(allRightSubsetsArray).filter(Objects::nonNull).collect(Collectors.toList());
//        System.out.println("Total number of susbets: " + output.size());
//        return output;
//        /////////////////////////////////////////////////////////////////
        return allRightSubsets;
    }

    /**
     * This method is used in the allRightSubsets method. This is a backtracking algorithm that finds all subsets
     * of an array that contain at least one end vertex (note that the root is already excluded from the array).
     * Inspiration: <a href="https://www.topcoder.com/thrive/articles/print-all-subset-for-set-backtracking-and-bitmasking-approach?utm_source=thrive&utm_campaign=thrive-feed&utm_medium=rss-feed">...</a>
     */
    public void subset(int i, List<String> subset, List<String> array, String x){
        if (i==N){
            subset.add(x); //add the terminal
            allRightSubsets.add(subset);
            return;
        }
        //add element i of the array
        List<String> withi = new ArrayList<String>(subset);
        withi.add(array.get(i));
        subset(i+1, withi, array, x);
        //don't add element i of the array
        List<String> withouti = new ArrayList<String>(subset);
        subset(i+1, withouti, array, x);
    }

    /**
     * First, less efficient method
     * For an input array, this returns the list of all subsets that doesn't contain the root and contains
     * at least one of the end vertices (note that the root was already excluded from the array).
     * I added a way to reduce this set to the minimal sets (see commented code) but doing this took way longer
     * than the time I win by reducing the model in this way, so this is not worth it.
     */
    public List<List<String>> allRightSubsetsFirst(List<String> array){
        N = array.size();
        allRightSubsets = new ArrayList<>();
        subsetFirst(0, new ArrayList<>(), array);

        //Minimal subsets
//        long startTime = System.currentTimeMillis();
//        //Create dictionary for vertex set and edges going into this set
//        Map<List<String>, List<MyWeightedEdge>> ingoing = new HashMap<>();
//        for (List<String> vertexSubset : allRightSubsets){
//            ingoing.put(vertexSubset, ingoingArcs(vertexSubset));
//        }
//        //Remove the vertex sets that are not minimal
//        List<String>[] allRightSubsetsArray = allRightSubsets.toArray(new List[0]);
//        for (int i = 0; i< allRightSubsetsArray.length; i++){
//            List<String> Ev = allRightSubsetsArray[i];
//            if(Ev == null) continue;
//            List<MyWeightedEdge> E = ingoing.get(Ev);
//            for (int j = 0; j< allRightSubsetsArray.length; j++){
//                List<String> Dv = allRightSubsetsArray[j];
//                if(Dv == null) continue;
//                List<MyWeightedEdge> D = ingoing.get(Dv);
//                if (D!=E){
//                    if (new HashSet<>(D).containsAll(E)){
//                        allRightSubsetsArray[j] = null;
//                    }else if (new HashSet<>(E).containsAll(D)){
//                        allRightSubsetsArray[i] = null;
//                    }
//                }
//            }
//        }
//        long endTime = System.currentTimeMillis();
//        timeMinimizing = endTime - startTime;
//
//        var output = Stream.of(allRightSubsetsArray).filter(Objects::nonNull).collect(Collectors.toList());
//        return output;
        return allRightSubsets;
    }

    /**
     * First, less efficient method
     * This method is used in the allRightSubsets method. This is a backtracking algorithm that finds all subsets
     * of an array that contain at least one end vertex (note that the root is already excluded from the array).
     * Inspiration: <a href="https://www.topcoder.com/thrive/articles/print-all-subset-for-set-backtracking-and-bitmasking-approach?utm_source=thrive&utm_campaign=thrive-feed&utm_medium=rss-feed">...</a>
     */
    public void subsetFirst(int i, List<String> subset, List<String> array){
        if (i==N){
            if (!subset.stream().filter(endVertices::contains).collect(Collectors.toList()).isEmpty()) {
                //if the set contains an end vertex, we add this subset
                allRightSubsets.add(subset);
            }
            return;
        }
        //add element i of the array
        List<String> withi = new ArrayList<String>(subset);
        withi.add(array.get(i));
        subsetFirst(i+1, withi, array);
        //don't add element i of the array
        List<String> withouti = new ArrayList<String>(subset);
        subsetFirst(i+1, withouti, array);
    }





































    /**
     * Dit is illegaal
     * Dit is het alternatief: die neemt niet alle subsets, maar alle paden van root naar terminal
     * For an input array, this returns the list of all subsets that doesn't contain the root and contains
     * at least one of the end vertices (note that the root was already excluded from the array).
     * I added a way to reduce this set to the minimal sets (see commented code) but doing this took way longer
     * than the time I win by reducing the model in this way, so this is not worth it.
     */
    public List<List<String>> allRightSubsetsAlternative(){
        allRightSubsetsAlternative = new ArrayList<>();
        for (String x : endVertices) {
            ArrayList<String> startArray = new ArrayList<>();
            startArray.add(x);
            subsetAlternative(x, startArray);
        }
        return allRightSubsetsAlternative;
    }

    /**
     * Dit is illegaal
     * Dit is voor het alternatief: die dus niet alle subsets neemt maar alle paden van root naar terminal
     * This method is used in the allRightSubsets method. This is a backtracking algorithm that finds all subsets
     * of an array that contain at least one end vertex (note that the root is already excluded from the array).
     * Inspiration: <a href="https://www.topcoder.com/thrive/articles/print-all-subset-for-set-backtracking-and-bitmasking-approach?utm_source=thrive&utm_campaign=thrive-feed&utm_medium=rss-feed">...</a>
     */
    public void subsetAlternative(String i, List<String> subset){
        for (MyWeightedEdge e : g.incomingEdgesOf(i)){ //for each incoming edge == for each parent
            allRightSubsetsAlternative.add(subset);
            String parent = g.getEdgeSource(e);
            if (parent=="root"){
                return;
            }
            //add this parent to path
            List<String> withparent = new ArrayList<String>(subset);
            withparent.add(parent);
            subsetAlternative(parent, withparent);
        }
    }
}
