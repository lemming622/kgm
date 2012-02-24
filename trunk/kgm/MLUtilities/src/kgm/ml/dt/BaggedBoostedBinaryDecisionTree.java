package kgm.ml.dt;
import kgm.ml.*;
import kgm.utility.*;
import java.util.*;

/**
 * A bagged boosted binary decision tree.
 * @author Kenton McHenry
 */
public class BaggedBoostedBinaryDecisionTree
{
	private Vector<BoostedBinaryDecisionTree> trees = new Vector<BoostedBinaryDecisionTree>();
	
	/**
	 * Train a bagged boosted decision tree.
	 * @param data the training data
	 * @param validation_data a validation set
	 * @param bags the number of bags to use
	 * @param bag_sample the data sample size to use for each bag
	 */
	public void build(Data data, Data validation_data, int bags, double bag_sample)
	{
		BoostedBinaryDecisionTree tree;
		Data bag_data;
		Vector<Double> Y1;
		
		trees.clear();
		
		for(int i=0; i<bags; i++){			
			bag_data = data.split(bag_sample, false);
			tree = new BoostedBinaryDecisionTree();
			tree.build(bag_data);
			
			if(tree != null){
				trees.add(tree);
				
				Y1 = classify(data.X);
	    	System.out.print((i+1) + ": " + Utility.round(data.accuracy(Y1), 2) + " (" + Utility.round(data.recall(Y1), 2) + ")");
	    	
	    	if(validation_data != null){
					Y1 = classify(validation_data.X);
	    		System.out.print(" [validation: " + Utility.round(validation_data.accuracy(Y1), 2) + " (" + Utility.round(validation_data.recall(Y1), 2) + ")]");
	    	}
	    	
	    	System.out.println();
			}
		}
	}
	
	/**
	 * Train a bagged boosted decision tree.
	 * @param data the training data
	 * @param bags the number of bags to use
	 * @param bag_sample the data sample size to use for each bag
	 */
	public void build(Data data, int bags, double bag_sample)
	{
		build(data, null, bags, bag_sample);
	}
	
	/**
	 * Train a bagged boosted decision tree.
	 * @param data the training data
	 * @param bags the number of bags to use
	 */
	public void build(Data data, int bags)
	{
		build(data, null, bags, 0.5);
	}
	
	/**
	 * Train a bagged boosted decision tree.
	 * @param data the training data
	 */
	public void build(Data data)
	{
		build(data, null, 5, 0.5);
	}
	
	/**
	 * Classify the given set of feature vectors.
	 * @param X the feature vectors to classify
	 * @return the feature vector classifications
	 */
	public Vector<Double> classify(Vector<double[]> X)
	{
		Vector<Double> Y = trees.get(0).classify(X);
		
		for(int i=1; i<trees.size(); i++){
			Data.plusEquals(Y, trees.get(i).classify(X));
		}
		
		Y = Data.signs(Y);
		
		return Y;
	}
	
	/**
	 * A simple main for debug purposes.
	 * @param args the command line arguments
	 */
	public static void main(String args[])
	{		
		BaggedBoostedBinaryDecisionTree tree = new BaggedBoostedBinaryDecisionTree();
		//Data data = new Data("data/debug.txt"); //data.print();
		Data data = new Data("data/diabetes/diabetes.txt"); //data.print();
		//Data data = new Data("data/cod-rna/cod-rna.txt.gz"); //data.print();
		Data training_data, validation_data = null;
		Vector<Double> Y;
				
		training_data = data.split(0.5);
		training_data.scale("tmp/scale.txt", true);
		validation_data = training_data.split(0.5);
    tree.build(training_data, validation_data, 5, 0.5);

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