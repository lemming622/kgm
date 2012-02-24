package kgm.ml.dt;
import kgm.ml.*;
import kgm.image.*;
import kgm.utility.*;
import java.awt.*;
import java.util.*;

/**
 * A boosted binary decision tree.
 * @author Kenton McHenry
 */
public class BoostedBinaryDecisionTree
{
	private Vector<BinaryDecisionTree> trees = new Vector<BinaryDecisionTree>();
	private double h0;
	public boolean VERBOSE = false;
	public boolean PLOT = false;
	
	/**
	 * Train a boosted decision tree.
	 * @param data the training data
	 * @param validation_data a validation set
	 * @param levels the levels for individual decision trees
	 * @param iterations the number of boosting iterations
	 * @param iteration_sample the data sample size to use at each iteration
	 * @param bias the data bias
	 */
	public void build(Data data, Data validation_data, int levels, int iterations, double iteration_sample, double bias)
	{
		Vector<double[]> X = data.X;
		Vector<Double> Y = data.Y;
		Vector<Double> W = data.W;
		Vector<Double> B = new Vector<Double>();
		Vector<Double> C = new Vector<Double>();
		Vector<Double> V = new Vector<Double>();
		Vector<Double> Y1;
		int positives, negatives;
		double positive_weight = 0, negative_weight = 0;
		double positive_bias = 1, negative_bias = 1;
		double accuracy, precision, recall;
		double validation_accuracy = 0, validation_precision = 0, validation_recall = 0;
		PlotPanel plot = null;

		if(PLOT){
			plot = PlotPanel.framedPlotPanel(600, 300, "Plot");
			plot.setAxis(1, iterations, 0, 100);
			plot.setDrawPoints(true);
		}
		
		//Initialize
		positives = 0;
		negatives = 0;
		
		for(int i=0; i<Y.size(); i++){
			if(Y.get(i) > 0){
				positives++;
			}else if(Y.get(i) < 0){
				negatives++;
			}
		}
		
		h0 = Math.log(((double)positives)/(Y.size()-positives));
		if(positives > 0) positive_weight = 0.5 * 1.0/positives;
		if(negatives > 0) negative_weight = 0.5 * 1.0/negatives;
		if(bias > 0) positive_bias = bias;
		if(bias < 0) negative_bias = -bias;
				
		for(int i=0; i<Y.size(); i++){
			B.add(1.0);
			C.add(h0);
			
			if(Y.get(i) > 0){
				W.set(i, positive_weight);
				B.set(i, positive_bias);
			}else if(Y.get(i) < 0){
				W.set(i, negative_weight);
				B.set(i, negative_bias);
			}
		}
		
		Data.timesEquals(W, B);
		Data.divideEquals(W, Data.sum(W));

		//Iteratively boost
		Data training_data = data;
		BinaryDecisionTree tree;
		
		for(int i=0; i<iterations; i++){
			if(iteration_sample < 1){		//Use only a sample of the training data
				training_data = data.split(iteration_sample, false);
				Data.divideEquals(training_data.W, Data.sum(training_data.W));
			}
						
			//tree = BinaryDecisionTree.stump(training_data);
			tree = new BinaryDecisionTree(); tree.build(training_data, levels);
			
			if(tree == null){
				System.out.println("Invalid tree!");
				break;
			}
						
			trees.add(tree);
			Data.plusEquals(C, tree.classify(X, true));
			
			//Use validation set to prevent over-fitting
			if(validation_data != null){
				Y1 = tree.classify(validation_data.X, false);
				validation_accuracy = 100*validation_data.accuracy(Y1);
				
				if(!V.isEmpty() && validation_accuracy < V.lastElement()){
					break;
				}else{
					V.add(validation_accuracy);
				}
			}
			
			//View progress
			if(VERBOSE || PLOT){
				Y1 = Data.signs(C);
				accuracy = 100*data.accuracy(Y1);
				precision = 100*data.precision(Y1);
				recall = 100*data.recall(Y1);
				
				if(validation_data != null){
					Y1 = tree.classify(validation_data.X, false);
					validation_accuracy = 100*validation_data.accuracy(Y1);
					validation_precision = 100*validation_data.precision(Y1);
					validation_recall = 100*validation_data.recall(Y1);
				}

		    if(VERBOSE){
		    	System.out.print((i+1) + ": " + Utility.round(accuracy, 2) + " (" + Utility.round(recall, 2) + ")");
		    	if(validation_data != null) System.out.print(" [validation: " + Utility.round(validation_accuracy, 2) + " (" + Utility.round(validation_recall, 2) + ")]");
		    	System.out.println();
		    }
		    
		    if(PLOT){
		    	plot.add(0, i+1, accuracy, Color.blue);
		    	plot.add(1, i+1, precision, new Color(0x8888ff));
		    	plot.add(2, i+1, recall, new Color(0xccccff));
		    	
		    	if(validation_data != null){
			    	plot.add(3, i+1, validation_accuracy, Color.red);
			    	plot.add(4, i+1, validation_precision, new Color(0xff8888));
			    	plot.add(5, i+1, validation_recall, new Color(0xffcccc));
		    	}
		    }		    
			}
			
			//Rebuild weights
			for(int j=0; j<W.size(); j++){
				if(C.get(j).isNaN()){
					//System.out.println("Found an invalid confidence: " + C.get(j));
					W.set(j, 0.0);
				}else{
					W.set(j, 1/(Math.exp(Y.get(j)*C.get(j))+1));
				}
				
				if(W.get(j).isNaN() || W.get(j).isInfinite()){
					System.out.println("Found an invalid weight: " + W.get(j) + "(" + C.get(j) + ")");
				}
			}
			
			Data.divideEquals(W, Data.sum(W));
			Data.timesEquals(W, B);
			Data.divideEquals(W, Data.sum(W));
		}
	}
	
	/**
	 * Train a boosted decision tree.
	 * @param data the training data
	 * @param validation_data a validation set
	 * @param levels the levels for individual decision trees
	 * @param iterations the number of boosting iterations
	 */
	public void build(Data data, Data validation_data, int levels, int iterations)
	{
		build(data, validation_data, levels, iterations, 1, 0);
	}
	
	/**
	 * Train a boosted decision tree.
	 * @param data the training data
	 * @param levels the levels for individual decision trees
	 * @param iterations the number of boosting iterations
	 */
	public void build(Data data, int levels, int iterations)
	{
		build(data, null, levels, iterations, 1, 0);
	}
	
	/**
	 * Train a boosted decision tree.
	 * @param data the training data
	 * @param validation_data a validation set
	 */
	public void build(Data data, Data validation_data)
	{
		build(data, validation_data, 2, 10, 1, 0);
	}
	
	/**
	 * Train a boosted decision tree.
	 * @param data the training data
	 */
	public void build(Data data)
	{
		build(data, null, 2, 10, 1, 0);
	}
	
	/**
	 * Classify the given set of feature vectors.
	 * @param X the feature vectors to classify
	 * @return the feature vector classifications
	 */
	public Vector<Double> classify(Vector<double[]> X)
	{
		Vector<Double> Y = new Vector<Double>();
		
		//Initialize
		for(int i=0; i<X.size(); i++){
			Y.add(h0);
		}
		
		//Accumulate confidences from each tree
		for(int i=0; i<trees.size(); i++){
			Data.plusEquals(Y, trees.get(i).classify(X, true));
		}
		
		//Extract signs
		Y = Data.signs(Y);
		
		return Y;
	}
	
	/**
	 * A simple main for debug purposes.
	 * @param args the command line arguments
	 */
	public static void main(String args[])
	{
		//Data data = new Data("data/debug.txt"); //data.print();
		Data data = new Data("data/diabetes/diabetes.txt"); //data.print();
		//Data data = new Data("data/cod-rna/cod-rna.txt.gz"); //data.print();
		Data training_data, validation_data = null;
		Vector<Double> Y;
		BoostedBinaryDecisionTree tree = new BoostedBinaryDecisionTree(); tree.VERBOSE = true; tree.PLOT = true;
				
		training_data = data.split(0.5);
		training_data.scale("tmp/scale.txt", true);
		validation_data = training_data.split(0.5);
    tree.build(training_data, validation_data, 2, 10, 0.5, 0.0);

    System.out.println("\n[Training data]");
    Y = tree.classify(training_data.X);
    Data.printClassDistribution(Y);
    System.out.println("Accuracy: " + Utility.round(100*training_data.accuracy(Y), 2));
       
    System.out.println("\n[Test data]");
    data.scale("tmp/scale.txt");
    Y = tree.classify(data.X);
    Data.printClassDistribution(Y);
    System.out.println("Accuracy: " + Utility.round(100*data.accuracy(Y), 2));
	}
}