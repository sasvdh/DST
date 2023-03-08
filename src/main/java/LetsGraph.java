import com.mxgraph.layout.*;
import com.mxgraph.layout.hierarchical.mxHierarchicalLayout;
import com.opencsv.CSVWriter;
import gurobi.GRB;
import gurobi.GRBException;
import org.jgrapht.Graph;
import org.jgrapht.GraphPath;
import org.jgrapht.Graphs;
import org.jgrapht.alg.interfaces.ShortestPathAlgorithm;
import org.jgrapht.alg.shortestpath.AllDirectedPaths;
import org.jgrapht.alg.shortestpath.DijkstraShortestPath;
import org.jgrapht.alg.util.Pair;
import org.jgrapht.generate.GnmRandomGraphGenerator;
import org.jgrapht.generate.GraphGenerator;
import org.jgrapht.graph.DefaultDirectedWeightedGraph;
import org.jgrapht.ext.JGraphXAdapter;

import com.mxgraph.util.mxCellRenderer;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import org.jgrapht.util.SupplierUtil;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

import static java.lang.Double.POSITIVE_INFINITY;
import static java.lang.Double.min;

public class LetsGraph {

    private static HashMap<String, HashMap<Set<String>, DefaultDirectedWeightedGraph<String,MyWeightedEdge>>> optimalTreesMap;
    public static HashMap<String, String> vertexNames;
    private static DefaultDirectedWeightedGraph<String, MyWeightedEdge> g;
    private static List<String> valuesThisIteration;
    private static long startTime;




    public static <CSVParser> void main(String[] args) throws IOException, GRBException {
        //////////////Create graph//////////////////////////////////////////////////////////////////////////////////////
        //This is the input for creating the random graph
        int nmin = 8; //Smallest graph size
        int nmax = 8; //Largest graph size
        int nstep = 10; //Step size for graph size
        int numberOfGraphs = 1; //number of different graphs of the same size
        int maxWeight = 1; //maximum arc weight
        int minNumberOfStartVertices = 1;
        int maxNumberOfStartVertices = 1; //test for 1,...,10 start vertices
        int sstep = 5; //Step size for end vertices
        int minNumberOfEndVertices = 2;
        int maxNumberOfEndVertices = 2; //test for 1,...,10 end vertices
        int estep = 5; //Step size for end vertices
        int numberOfStartEndCombis = 1; //number of different random start/end combis for each number of start/end vertices (i.e. number of random seeds)
        boolean tooManyVariables = true;

        /////////////Initialise csv file////////////////////////////////////////////////////////////add to string: WithoutDWandDCILP
        CSVWriter writer = new CSVWriter(new FileWriter("AlgorithmTimesAppr,nmin=" + nmin + ",nmax=" + nmax +
                ",nstep=" + nstep + ",ngraphs=" + numberOfGraphs + ",maxw=" + maxWeight +
                ",mins=" + minNumberOfStartVertices + ",maxs=" + maxNumberOfStartVertices + ",sstep=" + sstep +
                ",mine=" + minNumberOfEndVertices + ",maxe=" + maxNumberOfEndVertices + ",estep=" + estep +
                ",ncombi=" + numberOfStartEndCombis + ".csv"));
        CSVWriter writerdistr = new CSVWriter(new FileWriter("random10-100" + ".csv"));
        List<String[]> values = new ArrayList<>();
        addHeader();
        values.add(valuesThisIteration.toArray(new String[valuesThisIteration.size()]));
        long startTime = System.currentTimeMillis();
        String progress = "|";
        if (tooManyVariables) {
            progress += " ".repeat(numberOfGraphs * ((nmax - nmin) / nstep + 1) * numberOfStartEndCombis);
        }else {
            progress += " ".repeat(numberOfGraphs * ((nmax - nmin) / nstep + 1) * ((maxNumberOfStartVertices - minNumberOfStartVertices)/sstep + 1) *
                    ((maxNumberOfEndVertices-minNumberOfEndVertices)/estep + 1) * numberOfStartEndCombis);
        }
        progress += "|";
        progress += 0;
        System.out.println("Progress:");
        System.out.print(progress);
        for (int seedn = 1; seedn <= numberOfGraphs; seedn++) {
            for (int n = nmin; n <= nmax; n = n + nstep) { //Graph size
                int q = n * 88 / 53; //number of arcs
                int seedWeights = (seedn - 1) * n + n; //Unique seed for each iteration
                //Deze ene weer aan:
//                g = createRandomGraph(n, q, seedn, seedWeights, maxWeight);
//                visualiseGraph(g, "UnprocessedGraph" + n + "," + seedn, false);
                //Hier uit://
//              ASML graph:
//                DirectedWeightedPseudograph<String, MyWeightedEdge> g0 = createASMLGraph();
//                preprocess(g0);
//                visualiseGraph(g, "ProcessedGraph", false); //laat niet de dubbele arcs zien
                //Hier aan//

//                AllDirectedPaths<String, MyWeightedEdge> allpaths = new AllDirectedPaths<>(g);
//                List<GraphPath<String, MyWeightedEdge>> allpathslist = allpaths.getAllPaths(g.vertexSet(), g.vertexSet(), true, g.vertexSet().size());
//                System.out.println("Number of paths for n=" + n + ": " + allpathslist.size());

                DirectedWeightedPseudograph<String, MyWeightedEdge> g0 = createOwnGraph();
                preprocess(g0);
//                visualiseGraph(g, "CustomGraph", false);
                List<String> startVertices = new ArrayList<>();
                List<String> endVertices = new ArrayList<>();
                startVertices.add("s1");
//                startVertices.add("s2");
                endVertices.add("x1");
                endVertices.add("x2");
//                endVertices.add("x3");

        //Analyse flow distribution
//                List<Integer> incoming = new ArrayList<>();
//                List<Integer> outgoing = new ArrayList<>();
////                List<Integer> flow = new ArrayList<>();
//                List<Integer> indistr = new ArrayList<>(Collections.nCopies(40, 0));
//                List<Integer> outdistr = new ArrayList<>(Collections.nCopies(40, 0));
////                List<Integer> flowdistr = new ArrayList<>(Collections.nCopies(40, 0));
//                for (String v:g.vertexSet()){
//                    incoming.add(g.inDegreeOf(v));
//                    indistr.set(g.inDegreeOf(v) +20, indistr.get(g.inDegreeOf(v)+20) + 1);
//                    outgoing.add(g.outDegreeOf(v));
//                    outdistr.set(g.outDegreeOf(v)+20, outdistr.get(g.outDegreeOf(v)+20) + 1);
////                    flow.add(g.inDegreeOf(v) - g.outDegreeOf(v));
////                    flowdistr.set(g.inDegreeOf(v) - g.outDegreeOf(v) +15, flowdistr.get(g.inDegreeOf(v) - g.outDegreeOf(v)+15) + 1);
//                }
//                List<String> indistrs = indistr.stream().map(String::valueOf).collect(Collectors.toList());
//                List<String> outdistrs = outdistr.stream().map(String::valueOf).collect(Collectors.toList());
////                List<String> flowdistrs = flowdistr.stream().map(String::valueOf).collect(Collectors.toList());
//                indistrs.add(Integer.toString(n));
//                outdistrs.add(Integer.toString(n));
////                flowdistrs.add(Integer.toString(n));
//                indistrs.add(Integer.toString(seedn));
//                outdistrs.add(Integer.toString(seedn));
////                flowdistrs.add(Integer.toString(seedn));
//                indistrs.add("in");
//                outdistrs.add("out");
////                flowdistrs.add("flow");
//                List<String[]> distr = new ArrayList<>();
//                distr.add(indistrs.toArray(new String[indistrs.size()]));
//                distr.add(outdistrs.toArray(new String[outdistrs.size()]));
////                distr.add(flowdistrs.toArray(new String[flowdistrs.size()]));
//                writerdistr.writeAll(distr, false);
//                writerdistr.flush();


//                List<Integer> startArray = new ArrayList<>();
//                List<Integer> endArray = new ArrayList<>();
//                for (int i = 1; i < 10; i++){
//                    startArray.add(i);
//                    endArray.add(i);
//                }
//                for (int i = 10; i <= 50; i=i+5){
//                    startArray.add(i);
//                    endArray.add(i);
//                }

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                for (int seedSE = 1; seedSE <= numberOfStartEndCombis; seedSE++) { //Number of graphs for each #start/#end pair
//                    for (int numberOfStartVertices:startArray){
//                        for(int numberOfEndVertices:endArray){
                    for (int numberOfStartVertices = minNumberOfStartVertices; numberOfStartVertices <= maxNumberOfStartVertices; numberOfStartVertices += sstep) { //Number of start vertices
                        for (int numberOfEndVertices = minNumberOfEndVertices; numberOfEndVertices <= maxNumberOfEndVertices; numberOfEndVertices += estep) { //Number of end vertices
                            if (numberOfStartVertices <= n && numberOfEndVertices <= n) {
//                                System.out.println("nt=" + numberOfEndVertices);
//                                System.out.println("\nseed = " + seedn + ", n = " + n + ", seedSE = " + seedSE + ", ns = "
//                                        + numberOfStartVertices + ", ne = " + numberOfEndVertices);
//                                //**Deze line aan
////                                DefaultDirectedWeightedGraph<String, MyWeightedEdge> g = (DefaultDirectedWeightedGraph<String, MyWeightedEdge>) g0.clone();
//                                //Randomly choose start vertices. Then compute allShortestPathsGraph for each start vertex and combine these into a graph.
//                                //Then, randomly choose end vertices from this graph.
//                                int seedStart = (seedSE - 1) * numberOfStartVertices + numberOfStartVertices + (seedSE - 1) * numberOfEndVertices + numberOfEndVertices; //Unique seed for each iteration
//                                List<String> startVertices = new ArrayList<>();
//                                List<String> endVertices = new ArrayList<>();
//                                List<String> vertexList = new ArrayList<>(g.vertexSet());
//                                Random randStart = new Random(seedStart);
//                                boolean possible = false;
//                                int it = 0;
//                                while (possible == false && it < 2*n) { //repeat until we have a set of start vertices from which enough end vertices can be reached
//                                    startVertices = new ArrayList<>();
//                                    endVertices = new ArrayList<>();
//                                    while (startVertices.size() < numberOfStartVertices) { //only add element to start vertices if it is not in there yet, to avoid double vertices
//                                        int randInt = randStart.nextInt(n);
//                                        if (!startVertices.contains(vertexList.get(randInt))) {
//                                            startVertices.add(vertexList.get(randInt));
//                                        }
//                                    }
//                                    Set<String> possibleEnds = new HashSet<>();
//                                    for (int i = 0; i < numberOfStartVertices; i++) { //for each start vertex, compute shortest paths tree
//                                        DijkstraShortestPath<String, MyWeightedEdge> dspg = new DijkstraShortestPath<>(g);
//                                        ShortestPathAlgorithm.SingleSourcePaths<String, MyWeightedEdge> paths = dspg.getPaths(startVertices.get(i));
//                                        for (String v : vertexList) {
//                                            if (paths.getPath(v) != null) {
//                                                possibleEnds.add(v);
//                                            }
//                                        }
//                                    }
//                                    List<String> possibleEndsList = new ArrayList<>(possibleEnds);
//                                    int seedEnd = maxNumberOfStartVertices * numberOfStartEndCombis + maxNumberOfEndVertices * numberOfStartEndCombis +
//                                            (seedSE - 1) * numberOfStartVertices + numberOfStartVertices + (seedSE - 1) * numberOfEndVertices + numberOfEndVertices; //Unique seed for each iteration
//                                    Random randEnd = new Random(seedEnd);
//                                    if (possibleEnds.size() >= numberOfEndVertices) { //if there are enough vertices reachable from the start vertices
//                                        possible = true;
//                                        while (endVertices.size() < numberOfEndVertices) { //only add element to end vertices if it is not in there yet, to avoid double vertices
//                                            int randInt = randEnd.nextInt(possibleEndsList.size());
//                                            if (!endVertices.contains(possibleEndsList.get(randInt))) {
//                                                endVertices.add(possibleEndsList.get(randInt));
//                                            }
//                                        }
//                                    }
//                                    it++;
//                                }
//                                if(!possible){
//                                    break;
//                                }
//                                //Tot hier

//                                System.out.println("Number of paths for n=" + n + ": " + allpathslist.size());

//                                System.out.println(g + ", " + startVertices + ", " +endVertices);
//                                System.out.print("SeedSE=" +seedSE + ",nstart=" + numberOfStartVertices + ",nend=" + numberOfEndVertices );

                                g = addRoot(startVertices, g);
                                visualiseGraph(g, "CustomGraph", false);
//                    visualiseGraph(g, "ProcessedGraph" + n + "," + seedSE + "," + numberOfStartVertices + "," + numberOfEndVertices, false);


                                ///////////////Run algorithms///////////////////////////////////////////////////////////////
                                initialiseValues(n, numberOfStartVertices, numberOfEndVertices);
//                    System.out.println("________________________________________________");
//                    System.out.println("\nPaths algorithm: ");
                                DefaultDirectedWeightedGraph<String, MyWeightedEdge> solutionPathsForest = shortestPathsForest(g, startVertices, endVertices);
                                valuesThisIteration.add(Double.toString(getGraphWeight(solutionPathsForest)));
//                                System.out.print("," + getGraphWeight(solutionPathsForest) + "\n");
                            visualiseGraph(solutionPathsForest, "PathsAlgorithm" + n + "," + seedSE + "," + numberOfStartVertices + "," + numberOfEndVertices, true);
//                    System.out.println("________________________________________________");
//        System.out.println("\nSingular bunch algorithm: ");
                                DefaultDirectedWeightedGraph<String, MyWeightedEdge> solutionBunch = shortestBunch(g, endVertices);
                                valuesThisIteration.add(Double.toString(getGraphWeight(solutionBunch)));
//                                System.out.print("," + getGraphWeight(solutionBunch) + "\n");
//        System.out.println(solutionBunch);
        visualiseGraph(solutionBunch, "SingularBunchAlgorithm" + seedSE, true);
//        System.out.println("________________________________________________");
//        System.out.println("\nBunches algorithm: ");
                                DefaultDirectedWeightedGraph<String, MyWeightedEdge> solutionBunchesForest = shortestBunchesForest(g, startVertices, endVertices);
                                valuesThisIteration.add(Double.toString(getGraphWeight(solutionBunchesForest)));
//                                System.out.print("," + getGraphWeight(solutionBunchesForest) + "\n");
//        System.out.println(solutionBunchesForest);
        visualiseGraph(solutionBunchesForest, "MultipleBunchesAlgorithm" + seedSE, true);
//        System.out.println("________________________________________________");
//        System.out.println("\nPrim variation algorithm: ");
                                DefaultDirectedWeightedGraph<String, MyWeightedEdge> solutionPrimVar = primVariation(g, endVertices);
                                valuesThisIteration.add(Double.toString(getGraphWeight(solutionPrimVar)));
//                                System.out.print("," + getGraphWeight(solutionPrimVar) + "\n");
//        System.out.println(solutionPrimVar);
        visualiseGraph(solutionPrimVar, "PrimvariationAlgorithm" + seedSE, true);
//        System.out.println("________________________________________________");
//          System.out.println("\nFlow ILP: ");
                                FlowILP filp = new FlowILP();
                                Map<Pair<String, String>, Double> resFlowILP = filp.solveFlowILP(g, endVertices);
                                DefaultDirectedWeightedGraph<String, MyWeightedEdge> resGraphFlowILP = lpSolToGraph(g, resFlowILP);
//                                DefaultDirectedWeightedGraph<String, MyWeightedEdge> resGraphFlowILP = extractCrucialArcs(resGraphFlowILP0, endVertices);
                                resGraphFlowILP.removeVertex("root");
                                valuesThisIteration.add(Long.toString(filp.duration));
                                valuesThisIteration.add(Double.toString(getGraphWeight(resGraphFlowILP)));
////                                System.out.print("floilp: " + filp.duration);
////                                System.out.print("," + getGraphWeight(resGraphFlowILP) + "\n");
        visualiseGraph(resGraphFlowILP, "FlowILP" + seedSE, true);
////        System.out.println("________________________________________________");
////        //Dit is de paths MIP
////        System.out.println("\nPaths MIP: ");
//                                PathsMIP pmip = new PathsMIP();
//                                Map<Pair<String, String>, Double> resPMIP = new HashMap<>();
//                                try {
//                                    resPMIP = pmip.solvePathsMIP(g, endVertices, false, vertexNames);
//                                } catch (OutOfMemoryError E){
//                                    System.out.println("OutOfMemoryError pathsMIP at n=" + g.vertexSet().size() + ", #s=" + startVertices.size() + ", #e=" + endVertices.size());
//                                }
//                                DefaultDirectedWeightedGraph<String, MyWeightedEdge> resGraphPathsMIP = lpSolToGraph(g, resPMIP);
////                                DefaultDirectedWeightedGraph<String, MyWeightedEdge> resGraphPathsMIP = extractCrucialArcs(resGraphPathsMIP0, endVertices);
//                                resGraphPathsMIP.removeVertex("root");
//                                valuesThisIteration.add(Long.toString(pmip.duration));
//                                valuesThisIteration.add(Double.toString(getGraphWeight(resGraphPathsMIP)));

//                                AllDirectedPaths<String, MyWeightedEdge> allpaths = new AllDirectedPaths<>(g);
//                                List<GraphPath<String, MyWeightedEdge>> allpathslist = allpaths.getAllPaths(new HashSet<>(Arrays.asList("root")), new HashSet<>(endVertices), true, g.vertexSet().size());
//                                valuesThisIteration.add(Integer.toString(allpathslist.size()));
////                                System.out.print("pmip: " + pmip.duration);
////                                System.out.print("," + getGraphWeight(resGraphPathsMIP) + "\n");
////        visualiseGraph(resGraphPathsMIP, "PathsMIP" + seedSE, true);
////        System.out.println("________________________________________________");
////        //Dit is de directed cut ILP
////        System.out.println("\nDirected cut ILP: ");
//                                DirectedCutILP dcilp = new DirectedCutILP(g, endVertices);
//                                Map<Pair<String, String>, Double> resDirectedCutILP = dcilp.solveDirectedCutILP(g, endVertices);
//                                DefaultDirectedWeightedGraph<String, MyWeightedEdge> resGraphDirectedCutILP = lpSolToGraph(g, resDirectedCutILP);
////                                DefaultDirectedWeightedGraph<String, MyWeightedEdge> resGraphDirectedCutILP = extractCrucialArcs(resGraphDirectedCutILP0, endVertices);
//                                resGraphDirectedCutILP.removeVertex("root");
//                                valuesThisIteration.add(Long.toString(dcilp.duration));
//                                valuesThisIteration.add(Double.toString(getGraphWeight(resGraphDirectedCutILP)));
//                                System.out.print("dcilp: " + dcilp.duration);
//                                System.out.print("," + getGraphWeight(resGraphDirectedCutILP) + "\n");
//        visualiseGraph(resGraphDirectedCutILP, "DirectedCutILP" + seedSE, true);
//        System.out.println("________________________________________________");
//        System.out.println("\nDreyfus-Wagner algorithm: ");
//                                DefaultDirectedWeightedGraph<String, MyWeightedEdge> solutionDreyfus = dreyfusWagner(g, startVertices, endVertices);
//                                valuesThisIteration.add(Double.toString(getGraphWeight(solutionDreyfus)));
//                                System.out.print("," + getGraphWeight(solutionDreyfus) + "\n");
//        System.out.println(solutionDreyfus);
//        visualiseGraph(solutionDreyfus, "DreyfusWagner" + seedSE, true);
//        System.out.println("________________________________________________");
                                values.add(valuesThisIteration.toArray(new String[valuesThisIteration.size()])); //add values of this iteration to list
                                writer.writeAll(values, false);
                                writer.flush();
                                values = new ArrayList<>();
                                if (!tooManyVariables) {
                                    //Print progress of everything
                                    long endTime = System.currentTimeMillis();
                                    long duration = endTime - startTime;
                                    progress = progress.replaceFirst(" ", "|");
                                    progress = progress.substring(0, numberOfGraphs * ((nmax - nmin) / nstep + 1) *
                                            ((maxNumberOfStartVertices - minNumberOfStartVertices)/sstep + 1) *
                                            ((maxNumberOfEndVertices-minNumberOfEndVertices)/estep + 1) * numberOfStartEndCombis + 2);
                                    progress += duration;
                                    System.out.flush();
                                    System.out.print("\r");
                                    System.out.print(progress);
                                }
                            }
                        }
                    }
                    if (tooManyVariables) {
                        //Print progress of everything
                        long endTime = System.currentTimeMillis();
                        long duration = endTime - startTime;
                        progress = progress.replaceFirst(" ", "|");
                        progress = progress.substring(0, numberOfGraphs * ((nmax - nmin) / nstep + 1) *
                                numberOfStartEndCombis + 2);
                        progress += duration;
                        System.out.flush();
                        System.out.print("\r");
                        System.out.print(progress);
                    }
                }
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            }
        }
        //////////////Write to CSV files//////////////////////////////////////////////////////////////////////////
        writer.close();
        writerdistr.flush();
    }

    /**
     * Initialise values for result array
     */
    private static void initialiseValues(int n, int numberOfStartVertices, int numberOfEndVertices){
        valuesThisIteration = new ArrayList<>();
        valuesThisIteration.add(Integer.toString(n));
        valuesThisIteration.add(Integer.toString(numberOfStartVertices));
        valuesThisIteration.add(Integer.toString(numberOfEndVertices));
    }


    /**
     * Add header to array
     */
    private static void addHeader(){
        valuesThisIteration = new ArrayList<>();
        valuesThisIteration.add("n");
        valuesThisIteration.add("#start");
        valuesThisIteration.add("#end");
        valuesThisIteration.add("time pa");
        valuesThisIteration.add("value pa");
        valuesThisIteration.add("time sba");
        valuesThisIteration.add("value sba");
        valuesThisIteration.add("time mba");
        valuesThisIteration.add("value mba");
        valuesThisIteration.add("time pva");
        valuesThisIteration.add("value pva");
        valuesThisIteration.add("time filp");
        valuesThisIteration.add("value filpa");
        valuesThisIteration.add("time pmip");
        valuesThisIteration.add("value pmip");
        valuesThisIteration.add("time dcilp");
        valuesThisIteration.add("value dcilp");
        valuesThisIteration.add("time dwa");
        valuesThisIteration.add("value dwa");
    }

    /**
     * NOTE: This is no longer necessary since there are no 0 weights arcs except for the ones from the root
     * Solutions may contain zero weight arcs that are not crucial
     * The crucial arcs are extracted by running Dijkstra's shortest paths algorithm
     */
    public static DefaultDirectedWeightedGraph<String, MyWeightedEdge> extractCrucialArcs(
            DefaultDirectedWeightedGraph<String, MyWeightedEdge> gsol, List<String> endVertices){
        DijkstraShortestPath<String, MyWeightedEdge> dsp = new DijkstraShortestPath<>(gsol);
        ShortestPathAlgorithm.SingleSourcePaths<String, MyWeightedEdge> paths = dsp.getPaths("root");
        DefaultDirectedWeightedGraph<String, MyWeightedEdge> gresult = new DefaultDirectedWeightedGraph<>(MyWeightedEdge.class);
        for (String e:endVertices){
            Graphs.addAllVertices(gresult, paths.getPath(e).getVertexList());
            Graphs.addAllEdges(gresult, gsol, paths.getPath(e).getEdgeList());
        }
        return gresult;
    }


    /**
     * Add a root and zero weight arcs from the root to the start vertices
     */
    public static DefaultDirectedWeightedGraph<String, MyWeightedEdge> addRoot(List<String> startVertices,
                                                                               DefaultDirectedWeightedGraph<String, MyWeightedEdge> g){
        //Add a root and zero weight arcs from the root to the start vertices
        if (g.vertexSet().contains("root")) {
            g.removeVertex("root");
        }
        g.addVertex("root");
        for (String s:startVertices){
            MyWeightedEdge e = g.addEdge("root", s);
            g.setEdgeWeight(e, 0);
        }
        return g;
    }

    /**
     * Preprocessing of the graph:
     * Remove self loops
     * Remove double arcs by choosing the minimum weight one
     */
    public static void preprocess(DirectedWeightedPseudograph<String, MyWeightedEdge> g0) throws IOException {
        g = new DefaultDirectedWeightedGraph<String, MyWeightedEdge>(MyWeightedEdge.class);
        //Remove double arcs and self loops
        for (String v:g0.vertexSet()){
            HashMap<String, Double> lowestWeight = new HashMap<>(); //dictionary with target vertex and arc weight
            HashMap<String, MyWeightedEdge> lowestWeightArcs = new HashMap<>(); //dictionary with target and lowest weight arc
            ArrayList<MyWeightedEdge> toRemove = new ArrayList<>();
            for (MyWeightedEdge a:g0.outgoingEdgesOf(v)){
                String target = g0.getEdgeTarget(a);
                Double curweight = g0.getEdgeWeight(a);
                if (target.equals(v)){ //a is a self loop
                    toRemove.add(a); //need to remove a
                } else if (!lowestWeight.containsKey(target)) { //target is not in dictionary yet
                    lowestWeight.put(target, curweight); //add target to dictionary
                    lowestWeightArcs.put(target, a);
                } else if (lowestWeight.get(target) < curweight) { //lowest weight so far < arc weight
                    toRemove.add(a); //need to remove arc
                }else{ //lowest weight so far >= arc weight
                    toRemove.add(lowestWeightArcs.get(target)); //need to remove current lowest weight arc
                    lowestWeight.put(target, curweight); //put our current arc in dictionary
                    lowestWeightArcs.put(target, a);
                }
            }
            g0.removeAllEdges(toRemove); //remove arcs (we don't do this earlier because we iterate over the arcs)
        }
        Graphs.addAllVertices(g, g0.vertexSet());
        Graphs.addAllEdges(g, g0, g0.edgeSet());
    }

    /**
     * Output a PNG file containing a visualisation of the graph.
     * Input is the graph you want to visualise, a String with what you want to name the png file,
     * and a boolean stating if the graph you want to visualise is a tree. If the graph is a tree,
     * we visualise it as a hierarchical structure.
     */
    public static void visualiseGraph(Graph<String, MyWeightedEdge> g,
                                      String name, boolean tree) throws IOException {

        JGraphXAdapter<String, MyWeightedEdge> graphAdapter =
                new JGraphXAdapter<>(g);
        mxIGraphLayout layout;
        if (!tree) {
            layout = new mxOrganicLayout(graphAdapter); //visualise in organic layout
        }else{
            layout = new mxHierarchicalLayout(graphAdapter); //visualise as hierarchy if we have a forest
        }
        layout.execute(graphAdapter.getDefaultParent());

        BufferedImage image =
                mxCellRenderer.createBufferedImage(graphAdapter, null, 2, Color.WHITE, true, null);
        File imgFile = new File(name + ".png");
        ImageIO.write(image, "PNG", imgFile);
    }

    public static DefaultDirectedWeightedGraph<String, MyWeightedEdge> lpSolToGraph(DefaultDirectedWeightedGraph<String, MyWeightedEdge> g,
            Map<Pair<String, String>, Double> y) {
        DefaultDirectedWeightedGraph<String, MyWeightedEdge> resultGraph = new DefaultDirectedWeightedGraph<>(MyWeightedEdge.class);
        for(Pair<String, String> edge:y.keySet()){
            if (y.get(edge) > 0.5){
                if (!resultGraph.containsVertex(edge.getFirst())){
                    resultGraph.addVertex(edge.getFirst());
                }
                if (!resultGraph.containsVertex(edge.getSecond())){
                    resultGraph.addVertex(edge.getSecond());
                }
                MyWeightedEdge e = resultGraph.addEdge(edge.getFirst(), edge.getSecond());
                resultGraph.setEdgeWeight(e, g.getEdgeWeight(g.getEdge(edge.getFirst(), edge.getSecond())));
            }
        }
        return resultGraph;
    }

    /**
     * Algorithm that is a variation of Prim's algorithm for finding a minimum weight spanning tree.
     * Starting with just the root, in every iteration we add a minimum weight edge that connects
     * a vertex in the tree to a vertex not in the tree. This is repeated until all terminals are in the tree.
     * This tree might contain unnecessary vertices. Therefore, the final tree that is returned is the set
     * of shortest path in the first solution tree from the root to each of the terminals.
     */
    public static DefaultDirectedWeightedGraph<String, MyWeightedEdge> primVariation(
            DefaultDirectedWeightedGraph<String, MyWeightedEdge> g,
            List<String> endVertices){
        long startTime = System.currentTimeMillis();
        DefaultDirectedWeightedGraph<String, MyWeightedEdge> firstSolution = new DefaultDirectedWeightedGraph<>(MyWeightedEdge.class);
        firstSolution.addVertex("root");
        boolean allTerminalsVisited = false;
        List<String> vertices = new ArrayList<>(g.vertexSet()); //convert set to list
        List<Boolean> visited = new ArrayList<>(Arrays.asList(new Boolean[g.vertexSet().size()])); //List of visisted vertices
        Collections.fill(visited, Boolean.FALSE);
        List<Double> minDistance = new ArrayList<>(Arrays.asList(new Double[g.vertexSet().size()])); //Initialise list of infinity elements
        Collections.fill(minDistance, POSITIVE_INFINITY);
        List<String> parent = new ArrayList<>(Arrays.asList(new String[g.vertexSet().size()])); //List of parent if vertex is to be added
        List<String> neigboursOfRoot = new ArrayList<>(g.outDegreeOf("root"));
        for (MyWeightedEdge a: g.outgoingEdgesOf("root")){
            neigboursOfRoot.add(g.getEdgeTarget(a));
        }
        for (int i = 0; i < vertices.size(); i++){ //Initialise distance and visited vectors
            String v = vertices.get(i);
            if (v.equals("root")){
                visited.set(i, true);
                minDistance.set(i, 0.0);
            }else if(neigboursOfRoot.contains(v)){
                minDistance.set(i, g.getEdgeWeight(g.getEdge("root", v)));
                parent.set(i, "root");
            }
        }
        while (!allTerminalsVisited) { //while not all terminals are in the tree
            int vi = minElement(minDistance, visited); //find minimum distance arc
            String v = vertices.get(vi);
            firstSolution.addVertex(v); //add vertex
            MyWeightedEdge e = firstSolution.addEdge(parent.get(vi), v); //add arc
            firstSolution.setEdgeWeight(e, g.getEdgeWeight(g.getEdge(parent.get(vi), v)));
            minDistance.set(vi, 0.0);
            visited.set(vi, true);
            for (MyWeightedEdge a: g.outgoingEdgesOf(v)){ //update distance matrix
                String u = g.getEdgeTarget(a);
                int ui = vertices.indexOf(u);
                double arcWeight = g.getEdgeWeight(a);
                if (!visited.get(ui) && arcWeight < minDistance.get(ui)){ //if closer than previously, set new distance
                    minDistance.set(ui, arcWeight);
                    parent.set(ui, v);
                }
            }
            //check if all terminals have been visited
            allTerminalsVisited = true;
            for (String ev: endVertices) {
                int i = vertices.indexOf(ev);
                if (!visited.get(i)){
                    allTerminalsVisited = false;
                }
            }
        }

        //We return the tree that is the combination of shortest paths from the root to each of the terminals
        //This way, all unnecessary vertices are removed
        DefaultDirectedWeightedGraph<String, MyWeightedEdge> solutionTree = new DefaultDirectedWeightedGraph<>(MyWeightedEdge.class);
        //Compute the shortest paths tree rooted at "root"
        DijkstraShortestPath<String, MyWeightedEdge> dijkstra = new DijkstraShortestPath<>(g);
        ShortestPathAlgorithm.SingleSourcePaths<String, MyWeightedEdge> sptreer =  dijkstra.getPaths("root");
        for (String x:endVertices){
//            GraphPath<String,MyWeightedEdge> sp = DijkstraShortestPath.findPathBetween(firstSolution, "root", x);
            GraphPath<String,MyWeightedEdge> sp = sptreer.getPath(x);
            Graphs.addAllVertices(solutionTree, sp.getVertexList());
            addEdgesToGraph(sp.getEdgeList(), g, solutionTree);
        }

        long endTime = System.currentTimeMillis();
        long duration = (endTime - startTime);
        valuesThisIteration.add(Long.toString(duration));
//        System.out.print("Prim variation algorithm: " + duration);
//        System.out.println("Execution time bunches algorithm: " + duration + "ms");
        return solutionTree;
    }

    /**
     * Algorithm that finds the minimum of a list given that same index element is false in another list
     */
    public static int minElement(List<Double> values, List<Boolean> conditions){
        int index = -1;
        int min = Integer.MAX_VALUE;
        for (int i = 0; i < values.size(); i++) {
            if (values.get(i) < min && !conditions.get(i)){
                index = i;
            }
        }
        return index;
    }


    /**
     * Exact algorithm by Dreyfus-Wagner (variation). (Zie word document voor precieze omschrijving)
     */
    public static DefaultDirectedWeightedGraph<String, MyWeightedEdge> dreyfusWagner(
            DefaultDirectedWeightedGraph<String, MyWeightedEdge> g,
            List<String> startVertices, List<String> endVertices){
        startTime = System.currentTimeMillis();
//        //add a root vertex and connect this to all start vertices with 0 weight edges
//        g.addVertex("root");
//        for (String s:startVertices){
//            MyWeightedEdge e = g.addEdge("root", s);
//            g.setEdgeWeight(e, 0);
//        }
        DefaultDirectedWeightedGraph<String, MyWeightedEdge> solutionTree = new DefaultDirectedWeightedGraph<>(MyWeightedEdge.class);
        try {
            optimalTreesMap = new HashMap<>(); //initiate optimalTreesMap

            //compute and save all optimal ways to connect a vertex to a set
            computeOptimalTrees(g, endVertices);

//            for(String v: optimalTreesMap.keySet()){
//                HashMap<Set<String>, DefaultDirectedWeightedGraph<String, MyWeightedEdge>> map = optimalTreesMap.get(v);
//                System.out.println("For vertex " + v + ":");
//                for (Set<String> set: map.keySet()){
//                    System.out.println("     " + set + " has " + map.get(set));
//                }
//            }

            if (optimalTreesMap.containsKey("root")) {
                if (optimalTreesMap.get("root").containsKey(new HashSet<>(endVertices))) {
                    solutionTree = optimalTreesMap.get("root").get(new HashSet<>(endVertices));
                    solutionTree.removeVertex("root");
                }
            }
        } catch (OutOfMemoryError E){
            System.out.println("OutOfMemoryError at n=" + g.vertexSet().size() + ", #s=" + startVertices.size() + ", #e=" + endVertices.size());
        }

        long endTime = System.currentTimeMillis();
        long duration = (endTime - startTime);
        valuesThisIteration.add(Long.toString(duration));
//        System.out.print("dreyfus: " + duration);
//        System.out.println("Execution time Dreyfus-Wagner algorithm: " + duration + "ms");
        return solutionTree;
    }

    /**
     * Computes the optimal ways to connect a vertex v to a set X of vertices, for all sets X of size at least 2.
     * Saves all this information in a hashmap.
     */
    public static void computeOptimalTrees(DefaultDirectedWeightedGraph<String, MyWeightedEdge> g, List<String> endVertices){
        //Maps (set of vertices, vertex) -> (vertex, partitions)
        HashMap<String, HashMap<Set<String>, Pair<String, Set<Set<String>>>>> optimalTrees = new HashMap<>();

        //For all subsets: for all vertices not in this set: compute the optimal way to connect the vertex to the set
        for (int i = 2; i <= endVertices.size(); i++){
            List<List<String>> subsetsOfSizei = fixedSizeSubsets(endVertices, i); //compute all subsets of size i
            for (List<String> X:subsetsOfSizei){ //for all subsets of terminals of size i
                for (String v:g.vertexSet()){
                    if (System.currentTimeMillis() - startTime > 900000.0){
//                        break;
                    }
//                    if (!X.contains(v)){ //for v in V\X
                        //we hebbben de solutionTree niet nodig maar de dingen die we nodig hebben worden opgeslagen in de optimalTree map
                        //Maar moet wel returnen voor recursie in lightestTreeXv
                        DefaultDirectedWeightedGraph<String, MyWeightedEdge> solutionTree = lightestTreeXv(g, v, new HashSet<>(X));
//                    }
                }
            }
        }
        if (endVertices.size() < 2){
            //If we only have one end vertex, we cannot compute bipartitions of the end vertices
            //Therefore, we simply compute the shortest path from the root to the end vertex
            DefaultDirectedWeightedGraph<String, MyWeightedEdge> solutionSingleEndvertex = new DefaultDirectedWeightedGraph<>(MyWeightedEdge.class);
            DijkstraShortestPath<String, MyWeightedEdge> dspg = new DijkstraShortestPath<>(g);
            GraphPath<String, MyWeightedEdge> sp = dspg.getPath("root", endVertices.get(0));
            Graphs.addAllVertices(solutionSingleEndvertex, sp.getVertexList());
            Graphs.addAllEdges(solutionSingleEndvertex, g, sp.getEdgeList());
            HashMap<Set<String>, DefaultDirectedWeightedGraph<String, MyWeightedEdge>> tempmap = new HashMap<>();
            tempmap.put(new HashSet<>(endVertices), solutionSingleEndvertex);
            optimalTreesMap.put("root", tempmap);
        }
    }

    /**
     * Computes the cheapest way to connect a vertex v to a tree containing the terminals X
     */
    public static DefaultDirectedWeightedGraph<String, MyWeightedEdge> lightestTreeXv(
            DefaultDirectedWeightedGraph<String, MyWeightedEdge> g,
            String v,
            Set<String> X){
        DefaultDirectedWeightedGraph<String, MyWeightedEdge> solutionTree = new DefaultDirectedWeightedGraph<>(MyWeightedEdge.class);
        double minValue = POSITIVE_INFINITY;
        Set<Set<Set<String>>> bipartitions = bipartitions(X);
        String wchoice = "";
        for (Set<Set<String>> bip:bipartitions) { //for all nontrivial bipartitions
            for (String w : g.vertexSet()) { //for all vertices
                DefaultDirectedWeightedGraph<String, MyWeightedEdge> treeA = null;
                DefaultDirectedWeightedGraph<String, MyWeightedEdge> treeB = null;
                int i = 0;
                for (Set<String> part : bip) {
                    if (i==0) {
                        if (part.size() == 1) {
                            String single = part.iterator().next();
                            GraphPath<String,MyWeightedEdge> sp = DijkstraShortestPath.findPathBetween(g, w, single);
                            if (sp != null) {
                                treeA = new DefaultDirectedWeightedGraph<>(MyWeightedEdge.class);
                                Graphs.addAllVertices(treeA, sp.getVertexList());
                                addEdgesToGraph(sp.getEdgeList(), g, treeA);
                            }
                        }else {
                            if (optimalTreesMap.containsKey(w)) {
                                treeA = optimalTreesMap.get(w).get(part);
                            }
                        }
                        i++;
                    }else{
                        if (part.size() == 1) {
                            String single = part.iterator().next();
                            GraphPath<String,MyWeightedEdge> sp = DijkstraShortestPath.findPathBetween(g, w, single);
                            if (sp != null) {
                                treeB = new DefaultDirectedWeightedGraph<>(MyWeightedEdge.class);
                                Graphs.addAllVertices(treeB, sp.getVertexList());
                                addEdgesToGraph(sp.getEdgeList(), g, treeB);
                            }
                        }else {
                            if(optimalTreesMap.containsKey(w)) {
                                treeB = optimalTreesMap.get(w).get(part);
                            }
                        }
                    }
                }
                GraphPath<String, MyWeightedEdge> pathvw = DijkstraShortestPath.findPathBetween(g, v, w);
                if(treeA != null && treeB != null && pathvw != null) {
                    double value = getGraphWeight(treeA) + getGraphWeight(treeB) + pathvw.getWeight();
                    if (value < minValue) { //if this value is lower than the current minimum, make solutionTree this tree
                        minValue = value;
                        DefaultDirectedWeightedGraph<String, MyWeightedEdge> tempTree = new DefaultDirectedWeightedGraph<>(MyWeightedEdge.class);
                        Graphs.addAllVertices(tempTree, treeA.vertexSet());
                        Graphs.addAllVertices(tempTree, treeB.vertexSet());
                        Graphs.addAllVertices(tempTree, pathvw.getVertexList());
                        Graphs.addAllEdges(tempTree, treeA, treeA.edgeSet());
                        Graphs.addAllEdges(tempTree, treeB, treeB.edgeSet());
                        Graphs.addAllEdges(tempTree, g, pathvw.getEdgeList());
                        solutionTree = (DefaultDirectedWeightedGraph<String, MyWeightedEdge>) tempTree.clone();
                        wchoice = w;
                    }
                }
            }
        }
//        System.out.println(solutionTree);
        if(solutionTree.edgeSet().size() > 0) {//if we have a solution
            if (!optimalTreesMap.containsKey(v)) { //if v is not in the map yet
                HashMap<Set<String>, DefaultDirectedWeightedGraph<String, MyWeightedEdge>> inner = new HashMap<>();
                inner.put(X, solutionTree);
                optimalTreesMap.put(v, inner);
//                System.out.println("T(" + X + " + " + v + ") = " + solutionTree + " via w = " + wchoice);
            } else {
                if (!optimalTreesMap.get(v).containsKey(X)) { //if X is not in the map yet
                    optimalTreesMap.get(v).put(X, solutionTree);
//                    System.out.println("T(" + X + " + " + v + ") = " + solutionTree + " via w = " + wchoice);
                }
            }
        }
        return solutionTree;
    }


    /**
     * Given an input list of strings, this method returns the list of nontrivial bipartitions of the
     * elements of the list.
     * Deze heb ik van https://stackoverflow.com/questions/19892115/algorithm-to-find-all-possible-splits-of-a-list
     */
    public static Set<Set<Set<String>>> bipartitions(Set<String> set) {
        //For this all to work we need a list as input, but for the hashmaps to work we need sets. So let's convert
        List<List<List<String>>> bipartitionsList = new ArrayList<>();
        List<String> list = new ArrayList<>(set);
        List<String> listA;
        List<String> listB;
        int n = list.size(); //number of elements in the list
        int m = (1 << n) - 1; // m=111111...11=2^n-1 in binary, with n ones (shifting 1 n bits to the left gives 10...0 with n zeros)
        for (int i = 0; i < m - 1; ++i) { //all numbers from 0 up to and including 2^n-3; so all binaries with up to n entries
            if (i > (m ^ i) || i == 0) continue;
            //m ^ i is the "complement" of i. If i > m ^ i, this was already considered (with 0 and 1 exchanged)
            // i == 0 to ensure we exclude the trivial partition of {empty set, list}
            listA = new ArrayList<>();
            listB = new ArrayList<>();
            for (int j = 0; j < n; ++j) {
                if (((1 << j) & i) > 0) { //if the jth entry of i is a 1
                    listA.add(list.get(j));
                } else { //if the jth entry of i is a 0
                    listB.add(list.get(j));
                }
            }
            //Add a list containing the two lists to bipartitions
            List<List<String>> dummy = new ArrayList<>();
            dummy.add(listA);
            dummy.add(listB);
            bipartitionsList.add(dummy);
        }

        //convert the list of lists of lists to a set of sets of sets
        Set<Set<Set<String>>> bipartitions = new HashSet<>();
        for (List<List<String>> twolist : bipartitionsList){
            Set<Set<String>> setset = new HashSet<>();
            bipartitions.add(setset);
            for (List<String> onelist : twolist){
                setset.add(new HashSet<>(onelist));
            }
        }
        return bipartitions;
    }

    /**
     * Given an input list of strings, and a size, this method returns the list of subsets of the list
     * of the given size.
     * Deze heb ik vgm ook van internet geplukt
     */
    private static void getSubsets(List<String> inputList, int k, int idx,
                                   List<String> current, List<List<String>> solution) {
        if (current.size() == k) { //if the desired size has been reached, we add the current list to the solution
            solution.add(new ArrayList<>(current));
            return;
        }
        for (int i = idx; i < inputList.size(); i++) { //for start number until end of the inputList
            String x = inputList.get(i);
            current.add(x);
            getSubsets(inputList, k, i + 1, current, solution);
            current.remove(current.size() - 1);
        }
    }

    /**
     * Given an input list of strings, and a size k, this method returns the list of subsets of the list
     * of the given size.
     */
    public static List<List<String>> fixedSizeSubsets(List<String> inputList, int k) {
        List<List<String>> res = new ArrayList<>();
        getSubsets(inputList, k, 0, new ArrayList<>(), res); //initialize recursive algorithm
        return res;
    }

    /**
     * Algorithm to find a tree, starting at root, with endVertices as leaves. Note that the input
     * for this algorithm is the tree that contains a virtual root. We compute a bunch starting at
     * the root, where each bunch contains all end vertices it should contain. The bunch
     * is computed by trying each vertex in the graph as "middle vertex". If there is no path from r to v,
     * or no path from v to e for some e, this v is not valid. If v is valid,
     * we check if the total weight of the tree for this v is better than that of the previously best v.
     * If so, the current v becomes the new best v. When we've tried all vertices v as middle vertices,
     * we set the current tree to the solutionForest. The graph solutionForest is returned.
     */
    public static DefaultDirectedWeightedGraph<String, MyWeightedEdge> shortestBunch(
            DefaultDirectedWeightedGraph<String, MyWeightedEdge> g,
            List<String> endVertices){
        long startTime = System.currentTimeMillis();
        DefaultDirectedWeightedGraph<String, MyWeightedEdge> solutionForest = new DefaultDirectedWeightedGraph<>(MyWeightedEdge.class);
//        Graphs.addAllVertices(solutionForest, g.vertexSet()); //add the vertices from g to the solutionForest

        //Compute the shortest paths tree rooted at "root"
        DijkstraShortestPath<String, MyWeightedEdge> dijkstra = new DijkstraShortestPath<>(g);
        ShortestPathAlgorithm.SingleSourcePaths<String, MyWeightedEdge> sptreer =  dijkstra.getPaths("root");

        //Compute the bunches
        DefaultDirectedWeightedGraph<String, MyWeightedEdge> dummyForest;
        double currentWeight = POSITIVE_INFINITY;

        //we try to create a bunch for all vertices as middle vertex
        for (String v:g.vertexSet()){
            dummyForest = new DefaultDirectedWeightedGraph<>(MyWeightedEdge.class);
            boolean validv = true; //boolean that states if the choice of v is valid (i.e. there are paths to every end vertex it needs to connect to)
//            GraphPath<String, MyWeightedEdge> spsv = DijkstraShortestPath.findPathBetween(g, "root", v);
            GraphPath<String, MyWeightedEdge> spsv = sptreer.getPath(v);

            if (spsv != null) { //if there is a path from s to v
                //we create a dummy graph where we compute a bunch for this v
//                DefaultDirectedWeightedGraph<String, MyWeightedEdge> dummyForest = new DefaultDirectedWeightedGraph<>(MyWeightedEdge.class);
                Graphs.addAllVertices(dummyForest, spsv.getVertexList());
                Graphs.addAllEdges(dummyForest, g, spsv.getEdgeList());

                //Compute the shortest paths tree rooted at v
                ShortestPathAlgorithm.SingleSourcePaths<String, MyWeightedEdge> sptreev =  dijkstra.getPaths(v);
                //for all end vertices, find a shortest path from v to e
                for (String e:endVertices){
                    GraphPath<String, MyWeightedEdge> spve = sptreev.getPath(e);
//                    GraphPath<String, MyWeightedEdge> spve = DijkstraShortestPath.findPathBetween(g, v, e);
                    if (spve == null) { //if no path exists from v to e, v is not the right vertex for this bunch
                        validv = false;
                        break;
                    }
                    Graphs.addAllVertices(dummyForest, spve.getVertexList());
                    addEdgesToGraph(spve.getEdgeList(), g, dummyForest);

                }

                //if there is a path from v to every end vertex e, and the weight of
                //this solution is better than the current best, we set the solution to the current dummy
                if (validv) {
                    double dummyWeight = getGraphWeight(dummyForest);
                    if (dummyWeight < currentWeight) {
//                        System.out.println("Intermediate vertex = " + v);
//                        tempSolutionForest = (DefaultDirectedWeightedGraph<String, MyWeightedEdge>) dummyForest.clone();
                        solutionForest = (DefaultDirectedWeightedGraph<String, MyWeightedEdge>) dummyForest.clone();
                        currentWeight = dummyWeight;
                    }
                }
            }
        }
//        //We add the temporary solution to the final solution
//        Graphs.addAllVertices(solutionForest, dummyForest.vertexSet());
//        Graphs.addAllEdges(solutionForest, dummyForest, dummyForest.edgeSet());
        long endTime = System.currentTimeMillis();
        long duration = (endTime - startTime);
        valuesThisIteration.add(Long.toString(duration));
//        System.out.print("Shortest bunch algorithm: " + duration);
//        System.out.println("Execution time bunches algorithm: " + duration + "ms");
        return solutionForest;
    }

    /**
     * Algorithm to find a set of trees, starting at startVertices, with endVertices as leaves.
     * We first compute which start vertex s each of the end vertices e is closest to. We then compute
     * bunches for each s, where each bunch contains all end vertices it should contain. The bunches
     * are computed by trying each vertex in the graph as "middle vertex". If there is no path from s to v,
     * or no path from v to e for some e s should be connected to, this v is not valid. If v is valid,
     * we check if the total weight of the tree for this v is better than that of the previously best v.
     * If so, the current v becomes the new best v. When we've tried all vertices v as middle vertices,
     * we add the current tree to the solutionForest. The graph solutionForest is returned.
     */
    public static DefaultDirectedWeightedGraph<String, MyWeightedEdge> shortestBunchesForest(
            DefaultDirectedWeightedGraph<String, MyWeightedEdge> g,
            List<String> startVertices, List<String> endVertices){
        long startTime = System.currentTimeMillis();
        DefaultDirectedWeightedGraph<String, MyWeightedEdge> solutionForest = new DefaultDirectedWeightedGraph<>(MyWeightedEdge.class);
//        Graphs.addAllVertices(solutionForest, g.vertexSet()); //add the vertices from g to the solutionForest
        DijkstraShortestPath<String, MyWeightedEdge> dijkstra = new DijkstraShortestPath<>(g);

//        Map<String, ShortestPathAlgorithm.SingleSourcePaths<String, MyWeightedEdge>> sptPerStartVertex = new HashMap<>(); //dictionary containing the shortest path tree for each start vertex
//        for (String s:startVertices) {
//            //Compute the shortest paths tree rooted at each start vertex
//            ShortestPathAlgorithm.SingleSourcePaths<String, MyWeightedEdge> sptrees = dijkstra.getPaths(s);
//            sptPerStartVertex.put(s, sptrees);
//        }

        Map<String, ShortestPathAlgorithm.SingleSourcePaths<String, MyWeightedEdge>> sptPerVertex = new HashMap<>(); //dictionary containing the shortest path tree for each start vertex
        for (String v:g.vertexSet()) {
            //Compute the shortest paths tree rooted at each start vertex
            ShortestPathAlgorithm.SingleSourcePaths<String, MyWeightedEdge> sptrees = dijkstra.getPaths(v);
            sptPerVertex.put(v, sptrees);
        }

        //Decide which start vertex each end vertex should be connected to
        Map<String, String> startEndCombi = new HashMap<>(); //dictionary stating which start vertex each end vertex should be connected to
        for (String e:endVertices){ //for each end vertex, see which start vertex is closest
            double bestWeight = POSITIVE_INFINITY;
            for (String s:startVertices){
                GraphPath<String, MyWeightedEdge> sp = sptPerVertex.get(s).getPath(e);
//                GraphPath<String, MyWeightedEdge> sp = DijkstraShortestPath.findPathBetween(g, s, e);
                if (sp != null && sp.getWeight() < bestWeight){
                    startEndCombi.put(e, s);
                    bestWeight = sp.getWeight();
                }
            }
            if (bestWeight == POSITIVE_INFINITY){
                System.out.println("!There is no path between any start vertex and " + e + ". Solution will not contain " + e + " as end vertex.");
            }
        }

        //Compute the bunches
        for (String s:startVertices){ //for each start vertex we compute a bunch
            DefaultDirectedWeightedGraph<String, MyWeightedEdge> tempSolutionForest = new DefaultDirectedWeightedGraph<>(MyWeightedEdge.class);
            double currentWeight = POSITIVE_INFINITY;

            //we try to create a bunch for all vertices as middle vertex
            for (String v:g.vertexSet()){
                boolean validv = true; //boolean that states if the choice of v is valid (i.e. there are paths to every end vertex it needs to connect to)
//                GraphPath<String, MyWeightedEdge> spsv = DijkstraShortestPath.findPathBetween(g, s, v);
                GraphPath<String, MyWeightedEdge> spsv = sptPerVertex.get(s).getPath(v);

                if (spsv != null) { //if there is a path from s to v
                    //we create a dummy graph where we compute a bunch for this v
                    DefaultDirectedWeightedGraph<String, MyWeightedEdge> dummyForest = new DefaultDirectedWeightedGraph<>(MyWeightedEdge.class);
                    Graphs.addAllVertices(dummyForest, spsv.getVertexList());
                    Graphs.addAllEdges(dummyForest, g, spsv.getEdgeList());
//                    ShortestPathAlgorithm.SingleSourcePaths<String, MyWeightedEdge> sptreev = dijkstra.getPaths(s);

                    //for all end vertices that s should be connected to, find a shortest path from v to e
                    for (String e:endVertices){
                        if (startEndCombi.get(e).equals(s)){

                            GraphPath<String, MyWeightedEdge> spve = sptPerVertex.get(v).getPath(e);//sptreev.getPath(e);
//                            GraphPath<String, MyWeightedEdge> spve = DijkstraShortestPath.findPathBetween(g, v, e);
                            if (spve == null) { //if no path exists from v to e, v is not the right vertex for this bunch
                                validv = false;
                                break;
                            }
                            Graphs.addAllVertices(dummyForest, spve.getVertexList());
                            addEdgesToGraph(spve.getEdgeList(), g, dummyForest);
                        }
                    }

                    //if there is a path from v to every end vertex e that s needs to connect with, and the weight of
                    //this solution is better than the current best, we set the temporary solution to the current dummy
                    if (validv) {
                        double dummyWeight = getGraphWeight(dummyForest);
                        if (dummyWeight < currentWeight) {
                            tempSolutionForest = (DefaultDirectedWeightedGraph<String, MyWeightedEdge>) dummyForest.clone();
                            currentWeight = dummyWeight;
                        }
                    }
                }
            }
            //For each start vertex, we add the temporary solution to the final solution
            Graphs.addAllVertices(solutionForest, tempSolutionForest.vertexSet());
            Graphs.addAllEdges(solutionForest, tempSolutionForest, tempSolutionForest.edgeSet());
        }
        long endTime = System.currentTimeMillis();
        long duration = (endTime - startTime);
        valuesThisIteration.add(Long.toString(duration));
//        System.out.print("Shortest bunches algorithm: " + duration);
//        System.out.println("Execution time bunches algorithm: " + duration + "ms");
        return solutionForest;
    }

    /**
     * Algorithm to find a set of trees, starting at startVertices, with endVertices as leaves.
     * We add a vertex x to the graph. x is connected to all startVertices through zero-weight edges.
     * We find the shortest path from x to each of the endVertices. The edges of the shortest paths
     * (except for the one from x to the next one)  are each added to an edgeless graph on the vertices
     * of g (without x). They will form a set of trees rooted at the start vertices, containing the end
     * vertices. This graph (solutionForest) is returned.
     */
    public static DefaultDirectedWeightedGraph<String, MyWeightedEdge> shortestPathsForest(
            DefaultDirectedWeightedGraph<String, MyWeightedEdge> g,
            List<String> startVertices, List<String> endVertices){
        long startTime = System.currentTimeMillis();
        DefaultDirectedWeightedGraph<String, MyWeightedEdge> solutionForest = new DefaultDirectedWeightedGraph<>(MyWeightedEdge.class);

//        g.addVertex("x"); //add a vertex x to the graph
//        for(String s:startVertices){ //add a zero weight edge from x to each of the start vertices
//            MyWeightedEdge ee1 = g.addEdge("x", s);
//            g.setEdgeWeight(ee1, 0);
//        }

        //Compute the shortest paths tree rooted at "x"
        DijkstraShortestPath<String, MyWeightedEdge> dijkstra = new DijkstraShortestPath<>(g);
        ShortestPathAlgorithm.SingleSourcePaths<String, MyWeightedEdge> sptree =  dijkstra.getPaths("root");

        //Find a shortest path from x to each of the endVertices
        for(String e:endVertices){
            try {
                GraphPath<String, MyWeightedEdge> sp = sptree.getPath(e);
                List<MyWeightedEdge> spEdgeList = new ArrayList<>(sp.getEdgeList());
                spEdgeList.remove(0); //remove the edge coming from vertex x

                List<String> spVertexList = new ArrayList<>(sp.getVertexList());
                spVertexList.remove(0); //remove vertex x
                Graphs.addAllVertices(solutionForest, spVertexList); //add the vertices to the graph

                addEdgesToGraph(spEdgeList, g, solutionForest); //add the edges to the graph
            }
            catch(Exception ex){
                System.out.println("!There is no path between x and " + e + ". Solution does not contain " + e + " as end vertex.");
            }
        }
//        g.removeVertex("x");
        long endTime = System.currentTimeMillis();
        long duration = (endTime - startTime);
        valuesThisIteration.add(Long.toString(duration));
//        System.out.print("Paths algorithm: " + duration);
//        System.out.println("Execution time paths algorithm: " + duration + "ms");
        return solutionForest;
    }

    /**
     * Creates a random graph with n nodes. There is an edge between two vertices with probability p=probFactor/(n-1).
     * This means the expected number of edges will be n*(n-1)*p = n*probFactor. Furthermore, the random seed has to
     * be given as input.
     */
    public static DefaultDirectedWeightedGraph<String, MyWeightedEdge> createRandomGraph(
            int n, int probFactor, int seedGraph, int seedWeight, int maxWeight) {
        GraphGenerator<String, MyWeightedEdge, String> gen = new GnmRandomGraphGenerator<>(n, probFactor, seedGraph, false, false);
        DefaultDirectedWeightedGraph<String, MyWeightedEdge> g = new DefaultDirectedWeightedGraph<>(
                SupplierUtil.createStringSupplier(), SupplierUtil.createSupplier(MyWeightedEdge.class));
        gen.generateGraph(g);
        //assign random weights to the edges
        Random r = new Random(seedWeight);
        for (MyWeightedEdge e:g.edgeSet()){
            double randomValue = 0;
            while (randomValue == 0) {
                randomValue = maxWeight * r.nextDouble(); //random double that is not 0
            }
            g.setEdgeWeight(e, randomValue);
//            g.setEdgeWeight(e,1.0);
        }
        return g;
    }

    /**
     * Create a custom graph. Edit this to create whichever graph you want.
     */
    public static DirectedWeightedPseudograph<String, MyWeightedEdge> createASMLGraph(){
        //Create graph
        DirectedWeightedPseudograph<String, MyWeightedEdge> g = new DirectedWeightedPseudograph<>(MyWeightedEdge.class);

        ////////////////////////////WAFER-UP CONTEXT///////////////////////////////////////////////////////
        //Add vertices: Wafer-up Context
        g.addVertex("Flow");
        g.addVertex("LISJob");
        g.addVertex("Document");
        g.addVertex("PhysicalLayer");
        g.addVertex("ProcessStep");
        g.addVertex("Equipment");
        g.addVertex("WaferLot");
        g.addVertex("MeasureProcessJob");
        g.addVertex("ExposureProcessJob");
        g.addVertex("Wafer");
        g.addVertex("WaferMeasureProcessJob");
        g.addVertex("WaferExposureProcessJob");
        g.addVertex("Chuck");
        g.addVertex("WaferProcessJobId");
        g.addVertex("ProcessJobId");

        //Add edges: Wafer-Up Context
        g.addEdge("ProcessStep", "PhysicalLayer");
        g.addEdge("MeasureProcessJob", "ExposureProcessJob");
        g.addEdge("ExposureProcessJob", "MeasureProcessJob");
        g.addEdge("WaferLot", "Wafer");
        g.addEdge("MeasureProcessJob", "WaferMeasureProcessJob");
        g.addEdge("ExposureProcessJob", "WaferExposureProcessJob");
        g.addEdge("WaferMeasureProcessJob", "WaferExposureProcessJob");
        g.addEdge("WaferExposureProcessJob", "WaferMeasureProcessJob");
        g.addEdge("WaferExposureProcessJob", "Chuck");
        g.addEdge("WaferProcessJobId", "ProcessJobId");
        g.addEdge("MeasureProcessJob", "LISJob");
        g.addEdge("ExposureProcessJob", "LISJob");
        g.addEdge("MeasureProcessJob", "Flow");
        g.addEdge("ExposureProcessJob", "Flow");
        g.addEdge("MeasureProcessJob", "Document");
        g.addEdge("MeasureProcessJob", "Document");
        g.addEdge("ExposureProcessJob", "Document");
        g.addEdge("ExposureProcessJob", "Document");
        g.addEdge("MeasureProcessJob", "Equipment");
        g.addEdge("ExposureProcessJob", "Equipment");
        g.addEdge("MeasureProcessJob", "ProcessStep");
        g.addEdge("ExposureProcessJob", "ProcessStep");
        g.addEdge("MeasureProcessJob", "WaferLot");
        g.addEdge("ExposureProcessJob", "WaferLot");
        g.addEdge("WaferMeasureProcessJob", "Wafer");
        g.addEdge("WaferExposureProcessJob", "Wafer");
        g.addEdge("WaferMeasureProcessJob", "WaferProcessJobId");
        g.addEdge("WaferExposureProcessJob", "WaferProcessJobId");
        g.addEdge("WaferMeasureProcessJob", "ProcessJobId");
        g.addEdge("WaferExposureProcessJob", "ProcessJobId");


        ////////////////////////////WA-CONTEXT-EXTENSION///////////////////////////////////////////////////////
        //Add vertices: WA Context Extension
        g.addVertex("WaferLayout");
        g.addVertex("Die");
        g.addVertex("DieLayoutGrid");
        g.addVertex("Image");
        g.addVertex("Reticle");
        g.addVertex("MeasureLogicalWafer");
        g.addVertex("ExposureLogicalWafer");
        g.addVertex("ExposedField");
        g.addVertex("Field");
        g.addVertex("IntraFieldPoint");
        g.addVertex("MarkMeasurement");
        g.addVertex("CoefficientType");
        g.addVertex("AlignmentPhase");
        g.addVertex("AlignmentExposureProcessJob");
        g.addVertex("InterFieldPoint");

        //Add edges: WA Context Extension
        g.addEdge("ExposureProcessJob", "Reticle");
        g.addEdge("WaferExposureProcessJob", "MeasureLogicalWafer");
        g.addEdge("WaferExposureProcessJob", "ExposureLogicalWafer");
        g.addEdge("Reticle", "Image");
        g.addEdge("DieLayoutGrid", "Image");
        g.addEdge("WaferLayout", "DieLayoutGrid");
        g.addEdge("MeasureProcessJob", "WaferLayout");
        g.addEdge("ExposureProcessJob", "WaferLayout");
        g.addEdge("WaferLayout", "Die");
        g.addEdge("WaferLayout", "Field");
        g.addEdge("AlignmentExposureProcessJob", "ExposureProcessJob");
        g.addEdge("MeasureLogicalWafer", "MarkMeasurement");
        g.addEdge("ExposureLogicalWafer", "ExposedField");
        g.addEdge("ExposedField", "Image");
        g.addEdge("ExposedField", "Field");
        g.addEdge("IntraFieldPoint", "Field");
        g.addEdge("MarkMeasurement", "AlignmentPhase");
        g.addEdge("MarkMeasurement", "IntraFieldPoint");
        g.addEdge("AlignmentPhase", "CoefficientType");
        g.addEdge("MarkMeasurement", "InterFieldPoint");


        ////////////////////////////SUB-WAFER CONTEXT///////////////////////////////////////////////////////
        //Add vertices: Sub-Wafer Context
        g.addVertex("WaferComputeJob");
        g.addVertex("WaferPoint");

        //Add edges: Sub-Wafer Context
        g.addEdge("WaferPoint", "IntraFieldPoint");
        g.addEdge("WaferPoint", "InterFieldPoint");
        g.addEdge("WaferComputeJob", "WaferPoint");
        g.addEdge("WaferMeasureProcessJob", "WaferPoint");
        g.addEdge("MeasureLogicalWafer", "WaferPoint");


        ////////////////////////////OV CONTEXT EXTENSION///////////////////////////////////////////////////////
        //Add vertices: OV Context Extension
        g.addVertex("OverlayMeasureProcessJob");
        g.addVertex("YSWaferMeasurement");

        //Add edges: OV Context Extension
        g.addEdge("OverlayMeasureProcessJob", "MeasureProcessJob");
        g.addEdge("WaferMeasureProcessJob", "YSWaferMeasurement");
        g.addEdge("YSWaferMeasurement", "IntraFieldPoint");


        ////////////////////////////EP CONTEXT EXTENSION///////////////////////////////////////////////////////
        //Add vertices: EP Context Extension
        g.addVertex("epMeasureProcessJob");
        g.addVertex("Cutline");
        g.addVertex("PointLayout");
        g.addVertex("GridDieLayout");
        g.addVertex("FieldLayout");
        g.addVertex("ScatteredDieLayout");

        //Add edges: EP Context Extension
        g.addEdge("epMeasureProcessJob", "Cutline");
        g.addEdge("epMeasureProcessJob", "MeasureProcessJob");
        g.addEdge("Cutline", "PhysicalLayer");
        g.addEdge("Cutline", "PhysicalLayer");
        g.addEdge("Cutline", "IntraFieldPoint");
        g.addEdge("WaferComputeJob", "WaferMeasureProcessJob");
        g.addEdge("WaferComputeJob", "WaferExposureProcessJob");
        g.addEdge("WaferComputeJob", "PointLayout");
        g.addEdge("PointLayout", "IntraFieldPoint");
        g.addEdge("PointLayout", "InterFieldPoint");
        g.addEdge("MeasureLogicalWafer", "PointLayout");
        g.addEdge("WaferMeasureProcessJob", "PointLayout");
        g.addEdge("GridDieLayout", "Image");
        g.addEdge("WaferLayout", "GridDieLayout");
        g.addEdge("WaferLayout", "FieldLayout");
        g.addEdge("FieldLayout", "Field");
        g.addEdge("WaferLayout", "ScatteredDieLayout");
        g.addEdge("ScatteredDieLayout", "Die");


        ////////////////////////////YSMETRO CONTEXT EXTENSION///////////////////////////////////////////////////////
        //Add vertices: YSmetro Context Extension
        g.addVertex("YSmetroWaferMeasureProcessJob");
        g.addVertex("YSmetroMeasureProcessJob");
        g.addVertex("YSmetroSite");
        g.addVertex("YSmetroSetupFlow");
        g.addVertex("YSmetroModeledTargetContext");
        g.addVertex("YSmetroMeasurementContext");
        g.addVertex("YSmetroSubTarget");
        g.addVertex("YSmetroAcquirementContext");
        g.addVertex("YSmetroTarget");
        g.addVertex("YSmetroTargetToMeasure");
        g.addVertex("YSmetroMeasurementCluster");
        g.addVertex("YSmetroCluster");
        g.addVertex("YSmetroChannel");

        //Add edges: YSmetro Context Extension
        g.addEdge("WaferMeasureProcessJob", "YSmetroTarget");
        g.addEdge("YSmetroWaferMeasureProcessJob", "WaferMeasureProcessJob");
        g.addEdge("YSmetroMeasureProcessJob", "MeasureProcessJob");
        g.addEdge("YSmetroSetupFlow", "Flow");
        g.addEdge("YSmetroSetupFlow", "YSmetroModeledTargetContext");
        g.addEdge("YSmetroSubTarget", "YSmetroSite");
        g.addEdge("YSmetroSubTarget", "YSmetroAcquirementContext");
        g.addEdge("YSmetroTarget", "YSmetroSubTarget");
        g.addEdge("YSmetroTarget", "WaferPoint");
        g.addEdge("YSmetroTargetToMeasure", "YSmetroMeasurementContext");
        g.addEdge("YSmetroTarget", "YSmetroTargetToMeasure");
        g.addEdge("YSmetroAcquirementContext", "YSmetroChannel");
        g.addEdge("YSmetroTarget", "YSmetroCluster");
        g.addEdge("YSmetroTarget", "IntraFieldPoint");
        g.addEdge("YSmetroTargetToMeasure", "YSmetroMeasurementCluster");


        vertexNames = new HashMap<String, String>();
        for (String v:g.vertexSet()) {
            String vnew = "";
            for(int i = 0; i < v.length(); i++){
                if(Character.isUpperCase(v.charAt(i))){
                    vnew += v.charAt(i);
                }
            }
            vertexNames.put(v, vnew);
        }

//        //Check if there are duplicates in the names:
//        List<String> vertexList = new ArrayList<>(g.vertexSet());
//        for (int i = 0; i < vertexNames.size(); i++){
//            for (int j = 0; j < vertexNames.size(); j++){
//                if (i != j){
//                    if (vertexList.get(i) == vertexList.get(j)){
//                        System.out.println(vertexList.get(i) + " and " + vertexList.get(j) + " are both named " + vertexNames.get(vertexList.get(i)));
//                    }
//                }
//            }
//        }

//        //Weights are all 1:
//        for(MyWeightedEdge e:g.edgeSet()) {
//            g.setEdgeWeight(e, 1);
//        }
//        Weights are random between 0 and 1:
        Random r = new Random(1);
        for (MyWeightedEdge e:g.edgeSet()){
//            int randomValue = r.nextInt(maxWeight - 1) + 1; //random int between 1 and maxweight
            double randomValue = 0;
            while (randomValue == 0) {
                randomValue = 1 * r.nextDouble(); //random double that is not 0
            }
            g.setEdgeWeight(e, randomValue);
        }

//        //Analyse the arc distribution
//        System.out.println("Number of vertices: " + g.vertexSet().size());
//        System.out.println("Number of arcs: " + g.edgeSet().size());
////        System.out.println(g.vertexSet());
//        List<Integer> indegrees = new ArrayList<>();
//        List<Integer> outdegrees = new ArrayList<>();
//        List<Integer> flow = new ArrayList<>();
//        int selfloops = 0;
//        for (String v:g.vertexSet()) {
//            for (MyWeightedEdge a:g.incomingEdgesOf(v)) {
//                if (g.getEdgeSource(a) == v) {
//                    selfloops += 1;
//                }
//            }
//            indegrees.add(g.inDegreeOf(v));
//            outdegrees.add(g.outDegreeOf(v));
//            flow.add(g.inDegreeOf(v)-g.outDegreeOf(v));
//        }
//        System.out.println("Number of selfloops: " + selfloops);
//        System.out.println("Incoming arc distribution: " + indegrees);
//        System.out.println("Outgoing arc distribution: " + outdegrees);
//        System.out.println("Flow distribution: " + flow);


        return g;
    }

    /**
     * Create a custom graph. Edit this to create whichever graph you want.
     */
    public static DirectedWeightedPseudograph<String, MyWeightedEdge> createOwnGraph(){
        //Create graph
        DirectedWeightedPseudograph<String, MyWeightedEdge> g = new DirectedWeightedPseudograph<>(MyWeightedEdge.class);
        g.addVertex("s1");
        g.addVertex("v1");
        g.addVertex("v2");
        g.addVertex("v3");
        g.addVertex("v4");
        g.addVertex("v5");
        g.addVertex("v6");
        g.addVertex("x1");
        g.addVertex("x2");
        MyWeightedEdge e1 = g.addEdge("s1", "v1");
        g.setEdgeWeight(e1, 3);
        MyWeightedEdge e2 = g.addEdge("s1", "v2");
        g.setEdgeWeight(e2, 1);
        MyWeightedEdge e3 = g.addEdge("v1", "v3");
        g.setEdgeWeight(e3, 6);
        MyWeightedEdge e4 = g.addEdge("v1", "x1");
        g.setEdgeWeight(e4, 2);
        MyWeightedEdge e5 = g.addEdge("v1", "v4");
        g.setEdgeWeight(e5, 5);
        MyWeightedEdge e6 = g.addEdge("v1", "v2");
        g.setEdgeWeight(e6, 4);
        MyWeightedEdge e7 = g.addEdge("v2", "v5");
        g.setEdgeWeight(e7, 2);
        MyWeightedEdge e8 = g.addEdge("v4", "x2");
        g.setEdgeWeight(e8, 2);
        MyWeightedEdge e9 = g.addEdge("v6", "x2");
        g.setEdgeWeight(e9, 5);
        MyWeightedEdge e10 = g.addEdge("v4", "v2");
        g.setEdgeWeight(e10, 1);
        MyWeightedEdge e11 = g.addEdge("x1", "v6");
        g.setEdgeWeight(e11, 3);
        MyWeightedEdge e12 = g.addEdge("v2", "v4");
        g.setEdgeWeight(e12, 4);
//        MyWeightedEdge e13 = g.addEdge("v5", "x3");
//        g.setEdgeWeight(e13, 2);
//        MyWeightedEdge e14 = g.addEdge("v6", "x1");
//        g.setEdgeWeight(e14, 1);
//        MyWeightedEdge e15 = g.addEdge("x2", "x3");
//        g.setEdgeWeight(e15, 2);



//        //Voorbeeld die ik bij paths MIP wil hebben
//        g.addVertex("s1");
//        g.addVertex("v1");
//        g.addVertex("v2");
//        g.addVertex("v3");
//        g.addVertex("x1");
//        g.addVertex("x2");
//        // add edges
//        MyWeightedEdge e1 = g.addEdge("s1", "v1");
//        g.setEdgeWeight(e1, 1);
//        MyWeightedEdge e2 = g.addEdge("s1", "v2");
//        g.setEdgeWeight(e2, 1);
//        MyWeightedEdge e8 = g.addEdge("v1", "v3");
//        g.setEdgeWeight(e8, 1);
//        MyWeightedEdge e3 = g.addEdge("v2", "v3");
//        g.setEdgeWeight(e3, 1);
//        MyWeightedEdge e4 = g.addEdge("v3", "x1");
//        g.setEdgeWeight(e4, 1);
//        MyWeightedEdge e5 = g.addEdge("v3", "x2");
//        g.setEdgeWeight(e5, 1);


//        //Voorbeeld van de plaatjes
//        g.addVertex("s1");
//        g.addVertex("s2");
//        g.addVertex("v1");
//        g.addVertex("v2");
//        g.addVertex("v3");
//        g.addVertex("v4");
//        g.addVertex("v5");
//        g.addVertex("v6");
//        g.addVertex("v7");
//        g.addVertex("v8");
//        g.addVertex("x1");
//        g.addVertex("x2");
//        g.addVertex("x3");
//        // add edges
//        MyWeightedEdge e000 = g.addEdge("s1", "s1");
//        g.setEdgeWeight(e000, 150);
//        MyWeightedEdge e00 = g.addEdge("s1", "v1");
//        g.setEdgeWeight(e00, 100);
//        MyWeightedEdge e0 = g.addEdge("s1", "v1");
//        g.setEdgeWeight(e0, 99);
//        MyWeightedEdge e1 = g.addEdge("s1", "v1");
//        g.setEdgeWeight(e1, 3);
//        MyWeightedEdge e2 = g.addEdge("s1", "v2");
//        g.setEdgeWeight(e2, 1);
//        MyWeightedEdge e8 = g.addEdge("v8", "s2");
//        g.setEdgeWeight(e8, 4);
//        MyWeightedEdge e3 = g.addEdge("s2", "v6");
//        g.setEdgeWeight(e3, 2);
//        MyWeightedEdge e4 = g.addEdge("s2", "v7");
//        g.setEdgeWeight(e4, 3);
//        MyWeightedEdge e5 = g.addEdge("v4", "v1");
//        g.setEdgeWeight(e5, 6);
//        MyWeightedEdge e6 = g.addEdge("v1", "x1");
//        g.setEdgeWeight(e6, 2);
//        MyWeightedEdge e7 = g.addEdge("v3", "v1");
//        g.setEdgeWeight(e7, 5);
//        MyWeightedEdge e12 = g.addEdge("v2", "v3");
//        g.setEdgeWeight(e12, 4);
//        MyWeightedEdge e13 = g.addEdge("v3", "x2");
//        g.setEdgeWeight(e13, 2);
//        MyWeightedEdge e14 = g.addEdge("v6", "v3");
//        g.setEdgeWeight(e14, 3);
//        MyWeightedEdge e15 = g.addEdge("v6", "x3");
//        g.setEdgeWeight(e15, 4);
//        MyWeightedEdge e16 = g.addEdge("v7", "x3");
//        g.setEdgeWeight(e16, 5);
//        MyWeightedEdge e17 = g.addEdge("v5", "x1");
//        g.setEdgeWeight(e17, 5);


//        //Ik weet niet meer wat dit was, simpel voorbeeld
////        g.addVertex("r");
//        g.addVertex("s1");
//        g.addVertex("v1");
//        g.addVertex("v2");
//        g.addVertex("v3");
//        g.addVertex("x1");
//        g.addVertex("x2");
//        g.addVertex("x3");
//        // add edges
////        MyWeightedEdge e1 = g.addEdge("r", "s1");
////        g.setEdgeWeight(e1, 0);
//        MyWeightedEdge e2 = g.addEdge("s1", "v1");
//        g.setEdgeWeight(e2, 5);
//        MyWeightedEdge e8 = g.addEdge("s1", "v2");
//        g.setEdgeWeight(e8, 6);
//        MyWeightedEdge e3 = g.addEdge("v1", "x1");
//        g.setEdgeWeight(e3, 2);
//        MyWeightedEdge e4 = g.addEdge("v1", "v2");
//        g.setEdgeWeight(e4, 2);
//        MyWeightedEdge e5 = g.addEdge("v2", "x3");
//        g.setEdgeWeight(e5, 6);
//        MyWeightedEdge e6 = g.addEdge("v2", "v3");
//        g.setEdgeWeight(e6, 5);
//        MyWeightedEdge e7 = g.addEdge("v3", "x2");
//        g.setEdgeWeight(e7, 2);
//        MyWeightedEdge e12 = g.addEdge("v3", "x3");
//        g.setEdgeWeight(e12, 2);
////        MyWeightedEdge e9 = g.addEdge("5", "9");
////        g.setEdgeWeight(e9, 1);
////        MyWeightedEdge e10 = g.addEdge("5", "10");
////        g.setEdgeWeight(e10, 2);
////        MyWeightedEdge e11 = g.addEdge("6", "10");
////        g.setEdgeWeight(e11, 2);

        return g;
    }

    /**
     * Adds a list of edges to a graph
     */
    public static void addEdgesToGraph(List<MyWeightedEdge> edges,
                                       DefaultDirectedWeightedGraph<String, MyWeightedEdge> sourceGraph,
                                       DefaultDirectedWeightedGraph<String, MyWeightedEdge> destinationGraph){
        for(MyWeightedEdge d:edges){ //add the edges from the shortest path to the solutionTree
            if (!destinationGraph.containsEdge(sourceGraph.getEdgeSource(d), sourceGraph.getEdgeTarget(d))) { //only add edge if an edge doesn't already exist
                MyWeightedEdge d1 = destinationGraph.addEdge(sourceGraph.getEdgeSource(d), sourceGraph.getEdgeTarget(d));
                destinationGraph.setEdgeWeight(d1, sourceGraph.getEdgeWeight(d));
            }
        }
    }

    /**
     * Returns the total weight of a graph
     */
    public static double getGraphWeight(DefaultDirectedWeightedGraph<String, MyWeightedEdge> graph){
        double weight = 0;
        for(MyWeightedEdge e:graph.edgeSet()){
            weight += graph.getEdgeWeight(e);
        }
        return weight;
    }

}