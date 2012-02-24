package kgm.ml.svm;
import kgm.ml.*;
import libsvm.*;
import java.util.*;

/**
 * A Support Vector Machine classifier.
 * @author Kenton McHenry
 */
public class SVMClassifier
{
  /**
	 * Train an SVM classifier.
	 * @param data the training data
	 * @return the trained model
	 */
	public static svm_model train(Data data)
  {
    //Setup libsvm parameters
    svm_parameter param = new svm_parameter();
    param.svm_type = svm_parameter.C_SVC;
    param.kernel_type = svm_parameter.RBF;
    param.degree = 3;
    param.gamma = 0.5;  //1/k
    param.coef0 = 0;
    param.nu = 0.5;
    param.cache_size = 100;
    param.C = 1;
    param.eps = 1e-3;
    param.p = 0.1;
    param.shrinking = 1;
    param.probability = 0;
    param.nr_weight = 0;
    param.weight_label = new int[0];
    param.weight = new double[0];
    
    //Convert training data to libsvm's structure
    Vector<double[]> X = data.X;
    Vector<Double> Y = data.Y;
    Vector<svm_node[]> X_nodes = v2n(X);
    
    //Setup libsvm problem
    svm_problem prob = new svm_problem();
    prob.l = Y.size();
    prob.x = new svm_node[prob.l][];
    prob.y = new double[prob.l];    
    for(int i=0;i<prob.l;i++) prob.x[i] = X_nodes.get(i);
    for(int i=0;i<prob.l;i++) prob.y[i] = Y.get(i);
    
    //Train SVM
    svm_model model = svm.svm_train(prob,param);
    
    return model;
  }
  
  /**
	 * Train an SVM classifier.
	 * @param filename the file to store the model
	 * @param data the training data
	 * @return the trained model
	 */
	public static svm_model train(String filename, Data data)
	{
		svm_model model = train(data);
		
	  try{
	    svm.svm_save_model(filename, model);
	  }catch(Exception e) {e.printStackTrace();}
	  
	  return model;
	}

	/**
	 * Classify the given data using an SVM model.
	 * @param model the SVM model
	 * @param X the test data
	 * @return the classifications of the test data
	 */
	public static Vector<Double> predict(svm_model model, Vector<double[]> X)
	{
	  Vector<Double> Y = new Vector<Double>();
	  double y;
	  
	  //Convert training data to libsvm's structure
	  Vector<svm_node[]> X_nodes = v2n(X);
	  
	  //Test each example
	  for(int i=0; i<X_nodes.size(); i++){
	    y = svm.svm_predict(model, X_nodes.get(i));
	    Y.add(y);
	  }
	  
	  return Y;
	}

	/**
   * Classify the given data using an SVM model.
   * @param filename the file of the SVM model
   * @param X the test data
   * @return the classifications of the test data
   */
  public static Vector<Double> predict(String filename, Vector<double[]> X)
  {
  	Vector<Double> Y = null;
    svm_model model;
    
    try{
      model = svm.svm_load_model(filename);
      Y = predict(model, X);
    }catch(Exception e) {e.printStackTrace();}
    
    return Y;
  }
  
  /**
   * Convert a vector of N-D points into a vector svm_node's
   * @param X the data points
   * @return the vector of svm_nodes's
   */
  private static Vector<svm_node[]> v2n(Vector<double[]> X)
  {
    Vector<svm_node[]> X_nodes = new Vector<svm_node[]>();
    svm_node[] nodes;
    int d = X.get(0).length;    //Assume all points are of the same dimension
    
    for(int i=0; i<X.size(); i++){
      nodes = new svm_node[d];
      
      for(int j=0; j<X.get(i).length; j++){
        nodes[j] = new svm_node();
        nodes[j].index = j + 1;
        nodes[j].value = X.get(i)[j];
      }
      
      X_nodes.add(nodes);
    }
    
    return X_nodes;
  }
  
  /**
   * A simple main for debugging.
   */
  public static void main(String args[])
  {
		Data data = new Data("data/diabetes/diabetes.txt"); //data.print();
		Data training_data;
		Vector<Double> Y;

		training_data = data.split(0.2);
		training_data.scale("tmp/scale.txt", true);
    train("tmp/model.txt", training_data);
    
    System.out.println("\n[Training data]");
    Y = predict("tmp/model.txt", training_data.X);
    Data.printClassDistribution(Y);
    System.out.println("Accuracy: " + (100*training_data.accuracy(Y)));
    
    System.out.println("\n[Test data]");
    data.scale("tmp/scale.txt");
    Y = predict("tmp/model.txt", data.X);
    Data.printClassDistribution(Y);
    System.out.println("Accuracy: " + (100*data.accuracy(Y)));
  }
}