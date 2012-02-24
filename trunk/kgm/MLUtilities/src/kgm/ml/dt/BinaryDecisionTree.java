package kgm.ml.dt;
import kgm.ml.*;
import kgm.utility.*;
import edu.uci.ics.jung.graph.*;
import edu.uci.ics.jung.algorithms.layout.*;
import edu.uci.ics.jung.visualization.*;
import edu.uci.ics.jung.visualization.control.*;
import edu.uci.ics.jung.visualization.decorators.*;
import edu.uci.ics.jung.visualization.picking.*;
import edu.uci.ics.jung.visualization.renderers.Renderer;
import org.apache.commons.collections15.functors.*;
import java.awt.*;
import java.awt.geom.*;
import javax.swing.*;
import java.util.*;
/**
 * A binary decision tree.
 * @author Kenton McHenry
 */
public class BinaryDecisionTree
{
	private Node root = null;
	
	public BinaryDecisionTree() {}
	
	/**
	 * Class constructor.
	 * @param root the root node
	 */
	public BinaryDecisionTree(Node root)
	{
		this.root = root;
	}
	
	/**
	 * Create a one level decision tree.
	 * @param data the training data
	 * @return the created decision tree
	 */
	public static BinaryDecisionTree stump(Data data)
	{
		BinaryDecisionTree tree = null;
		Vector<double[]> X = data.X;
		Vector<Double> Y = data.Y;
		Vector<Double> W = data.W;
		int n = X.size();
		int d = X.get(0).length;
		boolean GOOD_DATA = false;
		
		//Check if the given data is adequate for creating a tree
		if(n > 0){
			for(int i=1; i<n; i++){
				if(!Y.get(i).equals(Y.get(0))){
					GOOD_DATA = true;
					break;
				}
			}
		}
		
		//Build the tree
		if(GOOD_DATA){
			Vector<Pair<Double,Integer>> Xorder = new Vector<Pair<Double,Integer>>();
			Vector<Double> Xsort = new Vector<Double>();
			Vector<Double> Ysort = new Vector<Double>();
			Vector<Double> Wsort = new Vector<Double>();
			Vector<double[]> thresholds = new Vector<double[]>();
			Vector<double[]> scores = new Vector<double[]>();
			Vector<int[]> labels = new Vector<int[]>();		//Labels for features greater than threshold
			double Wsum = 0;
			double epsilon = 1e-10;
			int index, at;
			double threshold = 0;
			double matches, score;
			int label = 0;
			
			//Initialize
			for(int i=0; i<n; i++){
				Wsum += W.get(i);
			}
						
			for(int i=0; i<=n; i++){
				scores.add(new double[d]);
				thresholds.add(new double[d]);
				labels.add(new int[d]);
			}
			
			//Build feature scores
			for(int j=0; j<d; j++){
				Xorder.clear();
				Xsort.clear();
				Ysort.clear();
				Wsort.clear();
				
				//Get jth feature and sort based on it
				for(int i=0; i<n; i++){
					Xorder.add(new Pair<Double,Integer>(X.get(i)[j], i));
				}
				
				Collections.sort(Xorder);
				
				for(int i=0; i<n; i++){
					index = Xorder.get(i).second;
					Xsort.add(X.get(index)[j]);
					Ysort.add(Y.get(index));
					Wsort.add(W.get(index));
				}
				
				//Initialize matches with first possible threshold
				threshold = Xsort.get(0) - epsilon;
				matches = 0;
				
				for(int i=0; i<n; i++){		//Everything is greater than threshold and classified as +1
					if(Ysort.get(i) > 0) matches += Wsort.get(i);
				}
				
				score = matches / Wsum;
				label = 1;
								
				if(score < 0.5){					//Flip label (and thus the score)
					score = 1.0 - score;
					label = -1;
				}
				
				thresholds.get(0)[j] = threshold;
				scores.get(0)[j] = score;
				labels.get(0)[j] = label;
				
				//Evaluate other thresholds
				at = 0;
				
				for(int i=0; i<n; i++){
					//Set threshold
					if(i < (n-1)){
						threshold = Xsort.get(i) + (Xsort.get(i+1)-Xsort.get(i))/2;
					}else{
						threshold = Xsort.get(i) + epsilon;
					}
										
					//Update matches
					while(at < n && Xsort.get(at) < threshold){
						if(Ysort.get(at) < 0){
							matches += Wsort.get(at);
						}else{
							matches -= Wsort.get(at);
						}
						
						at++;
					}
					
					score = matches / Wsum;
					label = 1;
					
					if(score < 0.5){				//Flip label (and thus the score)
						score = 1.0 - score;
						label = -1;
					}
					
					//System.out.println(j + "," + i + ": " + threshold + " => " + score);
					
					thresholds.get(i+1)[j] = threshold;
					scores.get(i+1)[j] = score;
					labels.get(i+1)[j] = label;
				}
			}
			
			//Find the maximum score
			double max_score = -Double.MAX_VALUE;
			int feature = -1;
			
			for(int j=0; j<d; j++){
				for(int i=0; i<=n; i++){
					if(scores.get(i)[j] >= max_score){
						max_score = scores.get(i)[j];
						feature = j;
						threshold = thresholds.get(i)[j];
						label = labels.get(i)[j];
					}
				}
			}
			
			//Debug
			if(false){				
				for(int j=0; j<d; j++){
					System.out.print(j + ": ");
					
					for(int i=0; i<=n; i++){
						System.out.print(scores.get(i)[j] + " ");
					}
					
					System.out.println();
				}
				
				System.out.println();
				System.out.println("Max score: " + max_score);
				System.out.println("Feature: " + feature);
				System.out.println("Threshold: " + threshold);
			}
			
			//Failsafe
			if(feature < 0){
				System.out.println("Invalid stump!");
				return null;
			}
			
			//Get left/right weight sums (greater than threshold will always go right)
			double Wsum_left = 0;
			double Wsum_right = 0;
			
			for(int i=0; i<n; i++){
				if(X.get(i)[feature] >= threshold){
					Wsum_right += W.get(i);
				}else{
					Wsum_left += W.get(i);
				}
			}
			
			//Get observed positive labels for each node and compute node confidence
			double root_positives = 0;
			double left_positives = 0;
			double right_positives = 0;
			double root_confidence, left_confidence, right_confidence;
			
			for(int i=0; i<n; i++){
				if(Y.get(i) > 0){
					root_positives += W.get(i);
				}
				
				if(X.get(i)[feature] >= threshold){
					if(Y.get(i) > 0){
						right_positives += W.get(i);
					}
				}else{
					if(Y.get(i) > 0){
						left_positives += W.get(i);
					}
				}
			}
			
			if(false){		//Debug, helpful when all weights are 1
				System.out.println("left positives: " + left_positives + " of " + Wsum_left);
				System.out.println("right positives: " + right_positives + " of " + Wsum_right);
				System.out.println();
			}
			
			root_positives /= Wsum;
			left_positives /= Wsum_left;
			right_positives /= Wsum_right;

			root_confidence = (Math.log(root_positives) - Math.log(1.0-root_positives)) * 0.5;
			left_confidence = (Math.log(left_positives) - Math.log(1.0-left_positives)) * 0.5;
			right_confidence = (Math.log(right_positives) - Math.log(1.0-right_positives)) * 0.5;

			//Build the decision tree
			Node left = new Node(-label, left_confidence);
			Node right = new Node(label, right_confidence);
			Node root = new Node(left, right, feature, threshold, root_confidence);
			tree = new BinaryDecisionTree(root);
		}
		
		return tree;
	}
	
	/**
	 * Grow a leaf by training another stump at this level.
	 * @param leaf a leaf in a binary decision tree
	 * @param data the training data
	 * @param levels the number of levels to grow
	 */
	public static void growLeaf(Node leaf, Data data, int levels)
	{
		BinaryDecisionTree tree;
		Node leaf_parent = null;
		Data data_left = new Data();
		Data data_right = new Data();
		
		//Build a new level
		tree = stump(data);
		
		if(tree != null){
			leaf_parent = leaf.parent;
			leaf.set(tree.root);
			leaf.parent = leaf_parent;
			
			//Split data
			for(int i=0; i<data.X.size(); i++){			
				if(data.X.get(i)[leaf.feature] >= leaf.threshold){
					data_right.X.add(data.X.get(i));
					data_right.Y.add(data.Y.get(i));
					data_right.W.add(data.W.get(i));
				}else{
					data_left.X.add(data.X.get(i));
					data_left.Y.add(data.Y.get(i));
					data_left.W.add(data.W.get(i));			
				}
			}
			
			//Grow further
			levels--;
			
			if(levels > 0){
				growLeaf(leaf.left, data_left, levels);
				growLeaf(leaf.right, data_right, levels);
			}
		}
	}
	
	/**
	 * Build a multi-level binary decision tree based on the given data.
	 * @param data the training data
	 * @param levels the number of levels in the tree
	 */
	public void build(Data data, int levels)
	{
		root = new Node();
		growLeaf(root, data, levels);
	}
	
	/**
	 * Classify the given feature vector.
	 * @param x the feature vector to classify
	 * @param USE_CONFIDENCE true if confidences should be used rather than labels
	 * @return the classification
	 */
	public double classify(double[] x, boolean USE_CONFIDENCE)
	{
		double y = 0;
		Node node = root;
		
		while(node.feature >= 0){
			if(x[node.feature] >= node.threshold){
				node = node.right;
			}else{
				node = node.left;
			}
		}
		
		if(USE_CONFIDENCE){
			y = node.confidence;
		}else{
			y = node.label;
		}
				
		return y;
	}
	
	/**
	 * Classify the given feature vector.
	 * @param x the feature vector to classify
	 * @return the classification
	 */
	public double classify(double[] x)
	{
		return classify(x, false);
	}
	
	/**
	 * Classify the given data.
	 * @param X the feature vectors to classify
	 * @param USE_CONFIDENCE true if confidences should be used rather than labels
	 * @return the feature vector classifications
	 */
	public Vector<Double> classify(Vector<double[]> X, boolean USE_CONFIDENCE)
	{
		Vector<Double> Y = new Vector<Double>();
		Node node;
		
		for(int i=0; i<X.size(); i++){
			Y.add(classify(X.get(i), USE_CONFIDENCE));
		}
		
		return Y;
	}
	
	/**
	 * Classify the given data.
	 * @param X the feature vectors to classify
	 * @return the feature vector classifications
	 */
	public Vector<Double> classify(Vector<double[]> X)
	{
		return classify(X, false);
	}
	
	/**
	 * Print the tree under the given node.
	 * @param node the root node
	 * @param count a count indicating what number the current node is in the tree
	 * @return the number of nodes already displayed
	 */
	public static int print(Node node, Integer count)
	{
		System.out.println("[Node-" + count + "]");
		
		if(node.feature >= 0){
			System.out.println("feature: " + node.feature);
			System.out.println("threshold: " + node.threshold);
			//System.out.println("confidence: " + node.confidence);
		}else{
			System.out.println("label: " + node.label);
			System.out.println("confidence: " + node.confidence);
		}
		
		System.out.println();
		
		count++;
		
		if(node.left != null) count = print(node.left, count);
		if(node.right != null) count = print(node.right, count);
		
		return count;
	}

	/**
	 * Print the stored tree.
	 */
	public void print()
	{
		print(root, 0);
	}
	
	/**
	 * Show the stored tree.
	 */
	public void show()
	{
		DelegateTree<Node,Integer> tree = new DelegateTree<Node,Integer>();
  	Stack<Node> stack = new Stack<Node>();
  	Node node;
  	int edge_count = 0;

  	tree.addVertex(root);
  	stack.push(root);
  	
  	while(!stack.isEmpty()){
  		node = stack.pop();
  		
  		if(node.left != null){
  			tree.addChild(edge_count++, node, node.left);
  			stack.push(node.left);
  		}
  		
  		if(node.right != null){
  			tree.addChild(edge_count++, node, node.right);
  			stack.push(node.right);
  		}
  	}
  	
  	TreeLayout<Node,Integer> layout = new TreeLayout(tree); 
  	DefaultModalGraphMouse gm = new DefaultModalGraphMouse();
  	gm.setMode(ModalGraphMouse.Mode.TRANSFORMING);
  	
  	VisualizationViewer<Node,Integer> vv = new VisualizationViewer<Node,Integer>(layout);
  	vv.setPreferredSize(new Dimension(600,600));
  	vv.setBackground(Color.white);
  	vv.setGraphMouse(gm);
  	vv.getRenderContext().setVertexShapeTransformer(new ConstantTransformer(new Rectangle2D.Float(-30,-6,60,12)));
  	vv.getRenderContext().setVertexFillPaintTransformer(new PickableVertexPaintTransformer<Node>(new MultiPickedState<Node>(), Color.white, Color.blue));
  	vv.getRenderContext().setVertexLabelTransformer(new ToStringLabeller());  	
  	vv.getRenderer().getVertexLabelRenderer().setPosition(Renderer.VertexLabel.Position.CNTR);
  	//vv.getRenderContext().setEdgeLabelTransformer(new ToStringLabeller());
  	vv.getRenderContext().setEdgeShapeTransformer(new EdgeShape.Line<Node,Integer>());
  	
  	JFrame frame = new JFrame("Tree");
  	frame.setSize(600, 600);
  	frame.add(vv);
  	frame.setVisible(true);
	}
	
	/**
	 * A main for debug purposes.
	 * @param args command line arguments
	 */
	public static void main(String[] args)
	{
		BinaryDecisionTree tree = new BinaryDecisionTree();
		//Data data = new Data("data/debug.txt"); //data.print();
		Data data = new Data("data/diabetes/diabetes.txt"); //data.print();
		//Data data = new Data("data/cod-rna/cod-rna.txt.gz"); //data.print();
		Data training_data;
		Vector<Double> Y;
		
		training_data = data.split(0.5);
		//Data.divideEquals(training_data.W, Data.sum(training_data.W));
		training_data.scale("tmp/scale.txt", true);
    //tree = stump(training_data);
    tree.build(training_data, 2);
    //tree.print();
    tree.show();
    
    System.out.println("[Training data]");
    Y = tree.classify(training_data.X);
    Data.printClassDistribution(Y);
    System.out.println("Accuracy: " + (100*training_data.accuracy(Y)));
       
    System.out.println("\n[Test data]");
    data.scale("tmp/scale.txt");
    Y = tree.classify(data.X);
    Data.printClassDistribution(Y);
    System.out.println("Accuracy: " + (100*data.accuracy(Y)));
	}
}