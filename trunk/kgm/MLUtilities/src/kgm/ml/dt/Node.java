package kgm.ml.dt;
import kgm.utility.*;

/**
 * A decision tree node.
 * @author Kenton McHenry
 */
public class Node
{
	public Node parent = null;
	public Node left = null;
	public Node right = null;
	
	public int feature = -1;
	public double threshold = 0;
	public int label = 0;
	public double confidence = 0;
	
	private double max_confidence = 10;
	
	public Node() {}
	
	/**
	 * Class constructor.
	 * @param label the label this node assigns
	 * @param confidence the confidence of this nodes label
	 */
	public Node(int label, double confidence)
	{
		this.label = label;
		
		if(!Double.isNaN(confidence)){
			if(Double.isInfinite(confidence)){
				if(confidence < 0){
					this.confidence = -max_confidence;
				}else{
					this.confidence = max_confidence;
				}
			}else{
				this.confidence = confidence;
			}
		}
	}

	/**
	 * Class constructor.
	 * @param left the left child node
	 * @param right the right child node
	 * @param feature the index of the feature to use
	 * @param threshold the threshold to use on the feature
	 * @param confidence the confidence of this nodes label
	 */
	public Node(Node left, Node right, int feature, double threshold, double confidence)
	{
		this.left = left;
		this.right = right;
		left.parent = this;
		right.parent = this;
		
		this.feature = feature;
		this.threshold = threshold;
		
		if(!Double.isNaN(confidence)){
			if(Double.isInfinite(confidence)){
				if(confidence < 0){
					this.confidence = -max_confidence;
				}else{
					this.confidence = max_confidence;
				}
			}else{
				this.confidence = confidence;
			}
		}	
	}
	
	/**
	 * Return a string representation of this node.
	 * @return the string representation of this node
	 */
	public String toString()
	{
		if(feature >= 0){
			return feature + " >= " + Utility.round(threshold,2);
		}else{
			if(Double.isInfinite(confidence)){
				if(confidence < 0){
					return "-Inf";
				}else{
					return "Inf";
				}
			}else{
				return Utility.round(confidence,2);
			}
		}
	}
	
	/**
	 * Set this node's values according to the given node.
	 * @param node the node to copy
	 */
	public void set(Node node)
	{
		left = node.left;
		right = node.right;
		parent = node.parent;
		
		feature = node.feature;
		threshold = node.threshold;
		label = node.label;
		confidence = node.confidence;
	}
}