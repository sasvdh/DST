import gurobi.*;
import org.jgrapht.alg.util.Pair;
import org.jgrapht.alg.util.Triple;
import org.jgrapht.graph.DefaultDirectedWeightedGraph;

import java.io.File;
import java.util.*;

public class FlowILP {
    List<String> endVertices;
    DefaultDirectedWeightedGraph<String, MyWeightedEdge> resultGraphFlow;
    long duration;

    /**
     * Solve the ILP described in https://arxiv.org/pdf/1111.5473.pdf
     * Given an input graph and a list of end vertices, this returns the tree containing the solution that was found.
     */
    public Map<Pair<String, String>, Double> solveFlowILP(DefaultDirectedWeightedGraph<String, MyWeightedEdge> g, List<String> endVertices) throws GRBException {
        long startTime = System.currentTimeMillis();
        this.endVertices = endVertices;
        int countvars = 0;
        int countcstrs = 0;
        //ILP opstellen:
        //Create environment:
        GRBEnv env = new GRBEnv(true);
        File logFile = new File("filp.log");
        env.set("logFile", "filp.log");
        env.set(GRB.IntParam.OutputFlag, 0);
        env.start();
        //Create empty model:
        GRBModel model = new GRBModel(env);
        //Create variables:
        Map<Pair<String, String>, GRBVar> y = new HashMap<>();
        for (MyWeightedEdge e:g.edgeSet()) {
            Pair<String, String> s = new Pair<>(g.getEdgeSource(e), g.getEdgeTarget(e));
            y.put(s, model.addVar(0.0, 1.0, g.getEdgeWeight(e), GRB.BINARY,
                    g.getEdgeSource(e) + "," + g.getEdgeTarget(e)));
            countvars += 1;
        }
        Map<Triple<String, String, String>, GRBVar> f = new HashMap<>();
        for (MyWeightedEdge e:g.edgeSet()) {
            for (String t : endVertices) {
                Triple<String, String, String> s = new Triple<>(t, g.getEdgeSource(e), g.getEdgeTarget(e));
                f.put(s, model.addVar(0.0, 1.0, 0.0, GRB.BINARY,
                        t + "," + g.getEdgeSource(e) + "," + g.getEdgeTarget(e)));
                countvars += 1;
            }
        }

        //Constraints aanmaken:
//        //flowpreservation
//        for (String s:endVertices) {
//            for(String v:g.vertexSet()) {
//                GRBLinExpr expr = new GRBLinExpr();
//                String cstr = "edgeTerminalLink(" + s + "," + v + ")";
//                for(MyWeightedEdge e:g.outgoingEdgesOf(v)) {
//                    expr.addTerm(1.0, f.get(new Triple<>(s, g.getEdgeSource(e), g.getEdgeTarget(e))));
//                }
//                for(MyWeightedEdge e:g.incomingEdgesOf(v)){
//                    expr.addTerm(-1.0, f.get(new Triple<>(s, g.getEdgeSource(e), g.getEdgeTarget(e))));
//                }
//                if (expr.size() > 0) {
//                    if (v.equals("root")) { //v is the root
//                        model.addConstr(expr, GRB.EQUAL, 1.0, cstr); // sum_outgoing fse - sum_ingoing fse = 1
//                    } else if (v.equals(s)) { //v is this end vertex
//                        model.addConstr(expr, GRB.EQUAL, -1.0, cstr); // sum_outgoing fse - sum_ingoing fse = -1
//                    } else {
//                        model.addConstr(expr, GRB.EQUAL, 0.0, cstr); // sum_outgoing fse - sum_ingoing fse = 0
//                    }
//                }
//            }
//        }
        //flowpreservation
        for (String s:endVertices) {
            for(String v:g.vertexSet()) {
                if (v.equals("root")) { //v is the root
                    //Constraint incoming edges
                    GRBLinExpr exprin = new GRBLinExpr();
                    String cstrin = "flowPreservationrin(" + s + "," + v + ")";
                    for(MyWeightedEdge e:g.incomingEdgesOf(v)){
                        exprin.addTerm(1.0, f.get(new Triple<>(s, g.getEdgeSource(e), g.getEdgeTarget(e))));
                    }
                    if (exprin.size() > 0) {
                        model.addConstr(exprin, GRB.EQUAL, 0.0, cstrin); // sum_incoming fsa = 0
                    }
                    //Constraint outgoing edges
                    GRBLinExpr exprout = new GRBLinExpr();
                    String cstrout = "flowPreservationrout(" + s + "," + v + ")";
                    for(MyWeightedEdge e:g.outgoingEdgesOf(v)){
                        exprout.addTerm(1.0, f.get(new Triple<>(s, g.getEdgeSource(e), g.getEdgeTarget(e))));
                    }
                    if (exprout.size() > 0) {
                        model.addConstr(exprout, GRB.EQUAL, 1.0, cstrout); // sum_outgoing fsa = 1
                        countcstrs += 1;
                    }
                } else if (v.equals(s)) { //v is this end vertex
                    //Constraint incoming edges
                    GRBLinExpr exprin = new GRBLinExpr();
                    String cstrin = "flowPreservationsin(" + s + "," + v + ")";
                    for(MyWeightedEdge e:g.incomingEdgesOf(v)){
                        exprin.addTerm(1.0, f.get(new Triple<>(s, g.getEdgeSource(e), g.getEdgeTarget(e))));
                    }
                    if (exprin.size() > 0) {
                        model.addConstr(exprin, GRB.EQUAL, 1.0, cstrin); // sum_incoming fsa = 1
                        countcstrs += 1;
                    }
                    //Constraint outgoing edges
                    GRBLinExpr exprout = new GRBLinExpr();
                    String cstrout = "flowPreservationsout(" + s + "," + v + ")";
                    for(MyWeightedEdge e:g.outgoingEdgesOf(v)){
                        exprout.addTerm(1.0, f.get(new Triple<>(s, g.getEdgeSource(e), g.getEdgeTarget(e))));
                    }
                    if (exprout.size() > 0) {
                        model.addConstr(exprout, GRB.EQUAL, 0.0, cstrout); // sum_outgoing fsa = 0
                        countcstrs += 1;
                    }
                } else { //not r or s
                    //flowin = flowout
                    GRBLinExpr expr = new GRBLinExpr();
                    String cstr = "flowPreservationNonrs(" + s + "," + v + ")";
                    for(MyWeightedEdge e:g.outgoingEdgesOf(v)) {
                        expr.addTerm(1.0, f.get(new Triple<>(s, g.getEdgeSource(e), g.getEdgeTarget(e))));
                    }
                    for(MyWeightedEdge e:g.incomingEdgesOf(v)){
                        expr.addTerm(-1.0, f.get(new Triple<>(s, g.getEdgeSource(e), g.getEdgeTarget(e))));
                    }
                    if (expr.size() > 0) {
                        model.addConstr(expr, GRB.EQUAL, 0.0, cstr); // sum_outgoing fsa - sum_incoming fsa = 0
                        countcstrs += 1;
                    }
                    //flowin <= 1
                    GRBLinExpr expr1 = new GRBLinExpr();
                    String cstr1 = "flowPreservationNonrs(" + s + "," + v + ")";
                    for(MyWeightedEdge e:g.incomingEdgesOf(v)){
                        expr1.addTerm(1.0, f.get(new Triple<>(s, g.getEdgeSource(e), g.getEdgeTarget(e))));
                    }
                    if (expr1.size() > 0) {
                        model.addConstr(expr1, GRB.LESS_EQUAL, 1.0, cstr1); // sum_outgoing fsa - sum_incoming fsa = 0
                        countcstrs += 1;
                    }
                }
            }
        }
        //edgeTerminalLink f_{s,e} <= y_e
        for (String s:endVertices) {
            for(MyWeightedEdge e: g.edgeSet()) {
                GRBLinExpr expr = new GRBLinExpr();
                String cstr = "edgeTerminalLink(" + s + "," + g.getEdgeSource(e) + "," + g.getEdgeTarget(e) + ")";
                expr.addTerm(1.0, f.get(new Triple<>(s, g.getEdgeSource(e), g.getEdgeTarget(e))));
                expr.addTerm(-1.0, y.get(new Pair<>(g.getEdgeSource(e), g.getEdgeTarget(e))));
                model.addConstr(expr, GRB.LESS_EQUAL, 0.0, cstr); // fse - ye <= 0
                countcstrs += 1;
            }
        }
        //deze constraint haal ik weg want het paper (Lasarre Hierarchy) had hem alleen voor analysis purposes
        //en het wordt er niet sneller op, alleen langzamer zelfs vgm
//        //treeterminal sum_{incoming v} y_a <= 1
//        for (String v:g.vertexSet()) {
//            GRBLinExpr expr = new GRBLinExpr();
//            String cstr = "edgeTerminalLink(" + v + ")";
//            for (MyWeightedEdge e:g.incomingEdgesOf(v)) {
//                expr.addTerm(1.0, y.get(new Pair<>(g.getEdgeSource(e), g.getEdgeTarget(e))));
//            }
//            model.addConstr(expr, GRB.LESS_EQUAL, 1.0, cstr); // sum_incoming ye <= 1
//        }
        //yzero y_e <= sum_{terminals s} f_{s,e}
        //We may get a solution that contains a 0 weight edge that is not used in the solution,
        //because f_se <= y_e does not ensure that y_e = 0 if all f_se are 0
        //Therefore, I add an extra constraint
        for (MyWeightedEdge e:g.edgeSet()) {
            GRBLinExpr expr = new GRBLinExpr();
            String cstr = "yzero(" + g.getEdgeSource(e) + "," + g.getEdgeTarget(e) + ")";
            expr.addTerm(1.0, y.get(new Pair<>(g.getEdgeSource(e), g.getEdgeTarget(e))));
            for (String s:endVertices) {
                expr.addTerm(-1.0, f.get(new Triple<>(s, g.getEdgeSource(e), g.getEdgeTarget(e))));
            }
            model.addConstr(expr, GRB.LESS_EQUAL, 0.0, cstr); // y_e - sum_terminals fse <= 0
            countcstrs += 1;
        }
//        System.out.println("variables: " + countvars + ", constraints: " + countcstrs);

        //Optimize model
        model.optimize();

//        model.computeIIS();
//        model.write("model.ilp");

//        System.out.println("Obj: " + model.get(GRB.DoubleAttr.ObjVal));

        //Save the results into a map
        Map<Pair<String, String>, Double> yresults = new HashMap<>();
        for(Pair<String, String> ye:y.keySet()) {
            yresults.put(ye, y.get(ye).get(GRB.DoubleAttr.X));
        }

        // Dispose of model and environment
        model.dispose();
        env.dispose();

        long endTime = System.currentTimeMillis();
        this.duration = (endTime - startTime);
//        System.out.println("Execution time solving flow ILP: " + duration + "ms");

        return yresults;
    }
}
