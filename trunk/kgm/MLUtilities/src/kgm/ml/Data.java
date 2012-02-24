package kgm.ml;
import kgm.utility.*;
import java.io.*;
import java.util.*;

/**
 * A class to hold example data.
 * @author Kenton McHenry
 */
public class Data
{
	public Vector<double[]> X = new Vector<double[]>();
	public Vector<Double> Y = new Vector<Double>();
	public Vector<Double> W = new Vector<Double>();
	
	public Data() {}
	
	/**
	 * Class constructor.
	 * @param filename the data file to load
	 */
	public Data(String filename)
	{
		load(filename);
	}
	
	/**
	 * Class constructor.
	 * @param X the feature vectors to store
	 * @param Y the feature vector labels
	 * @param W the feature vector weights
	 */
	public Data(Vector<double[]> X, Vector<Double> Y, Vector<Double> W)
	{
		this.X = X;
		this.Y = Y;
		this.W = W;
	}
	
	/**
	 * Clear current data.
	 */
	public void clear()
	{
		X.clear();
		Y.clear();
		W.clear();
	}
	
	/**
	 * Get the number of examples.
	 * @return the number of examples
	 */
	public int size()
	{
		return X.size();
	}

	/**
	 * Load a data file.
	 * @param filename the name of the data file
	 */
	public void load(String filename)
	{
		Scanner ins, scanner;
		String line, tmp;
		double[] x;
		double xi, y;
		int i, tmpi;
		int features = 0;
		boolean FIRST_EXAMPLE = true;
			
		clear();
		
		//Read the file in once to determine the feature vector length
		ins = new Scanner(Utility.getInputStream(filename));
		
		while(ins.hasNextLine()){
			line = ins.nextLine();
			scanner = new Scanner(line);
			
			y = scanner.nextDouble();
			
			while(scanner.hasNext()){
				tmp = scanner.next();
				tmpi = tmp.indexOf(':');
				
				if(tmpi >= 0){						//Sparse representation
					i = Integer.valueOf(tmp.substring(0, tmpi));
					if(i > features) features = i;
				}else if(FIRST_EXAMPLE){	//Dense representation
					features++;
				}
			}
			
			FIRST_EXAMPLE = false;
		}
		
		ins.close();

		//Read in the files contents
		ins = new Scanner(Utility.getInputStream(filename));
		
		while(ins.hasNextLine()){
			line = ins.nextLine();
			scanner = new Scanner(line);
			
			y = scanner.nextDouble();
			x = new double[features];
			i = 0;
			
			while(scanner.hasNext()){
				tmp = scanner.next();
				tmpi = tmp.indexOf(':');
				
				if(tmpi >= 0){		//Sparse representation
					i = Integer.valueOf(tmp.substring(0, tmpi)) - 1;
					xi = Double.valueOf(tmp.substring(tmpi+1));
					x[i] = xi;
				}else{						//Dense representation
					xi = Double.valueOf(tmp);
					x[i++] = xi;
				}
			}
			
			add(x, y);
		}
		
		ins.close();
	}

	/**
	 * Add example to data.
	 * @param x the example feature vector to add
	 * @param y the example class
	 * @param w the example weight
	 */
	public void add(double[] x, double y, double w)
	{
		X.add(x);
		Y.add(y);
		W.add(w);
	}

	/**
	 * Add example to data.
	 * @param x the example feature vector to add
	 * @param y the example class
	 */
	public void add(double[] x, double y)
	{
		add(x, y, 1);
	}
	
	/**
	 * Set data based on another data instance.
	 * @param data the data to use
	 */
	public void set(Data data)
	{
		X = data.X;
		Y = data.Y;
		W = data.W;
	}

	/**
	 * Scale data between 0 and 1.
	 * @return the min/max values for each feature
	 */
	public Pair<double[],double[]> scale()
	{
		int features = X.get(0).length;
		double[] min = new double[features];
		double[] max = new double[features];
		
		//Initialize feature min/max values
		for(int i=0; i<features; i++){
			min[i] = Double.MAX_VALUE;
			max[i] = -Double.MAX_VALUE;
		}
		
		//Find min/max value for each feature
		for(int i=0; i<X.size(); i++){
			for(int j=0; j<X.get(i).length; j++){
				if(X.get(i)[j] < min[j]) min[j] = X.get(i)[j];
				if(X.get(i)[j] > max[j]) max[j] = X.get(i)[j];
			}
		}
		
		//Scale the data
		scale(min, max);
		
		return new Pair<double[],double[]>(min,max);
	}
	
	/**
	 * Scale data between 0 and 1 based on min/max values.
	 * @param min the minimum values for each feature
	 * @param max the maximum values for each feature
	 */
	public void scale(double[] min, double[] max)
	{
		for(int i=0; i<X.size(); i++){
			for(int j=0; j<X.get(i).length; j++){
				X.get(i)[j] = (X.get(i)[j]-min[j]) / max[j];
			}
		}
	}
	
	/**
	 * Scale data between 0 and 1 based on min/max values.
	 * @param filename a file holding these min/max values
	 * @param SAVE true if the min/max values should be saved to this file
	 * @return the min/max values for each feature
	 */
	public Pair<double[],double[]> scale(String filename, boolean SAVE)
	{
		Pair<double[],double[]> minmax = null;
		double[] min;
		double[] max;
		
		if(SAVE){
			minmax = scale();
			min = minmax.first;
			max = minmax.second;
			
	    try{
	      FileWriter outs = new FileWriter(filename);
	      
	      for(int i=0; i<min.length; i++){
	        outs.write(min[i] + " ");
	      }
	      
	      outs.write("\n");
	      
	      for(int i=0; i<max.length; i++){
	        outs.write(max[i] + " ");
	      }
	      
	      outs.close();
	    }catch(Exception e) {e.printStackTrace();}
		}else{
			Vector<Double> tmpv = new Vector<Double>();
			String line;
			
			try{
				Scanner ins = new Scanner(new File(filename));
				Scanner scanner;
				
				//Read in min values
				line = ins.nextLine();
				scanner = new Scanner(line);
				tmpv.clear();
				
				while(scanner.hasNext()){
					tmpv.add(scanner.nextDouble());
				}
				
				min = new double[tmpv.size()];
				
				for(int i=0; i<tmpv.size(); i++){
					min[i] = tmpv.get(i);
				}
				
				//Read in max values
				line = ins.nextLine();
				scanner = new Scanner(line);
				tmpv.clear();
				
				while(scanner.hasNext()){
					tmpv.add(scanner.nextDouble());
				}
				
				max = new double[tmpv.size()];
				
				for(int i=0; i<tmpv.size(); i++){
					max[i] = tmpv.get(i);
				}
				
				ins.close();
			
				//Scale the data
				scale(min, max);
			}catch(Exception e) {e.printStackTrace();}
		}
		
		return minmax;
	}
	
	/**
	 * Scale data between 0 and 1 based on min/max values.
	 * @param filename a file holding these min/max values
	 * @return the min/max values for each feature
	 */
	public Pair<double[],double[]> scale(String filename)
	{
		return scale(filename, false);
	}

	/**
	 * Split off a portion of the data set.
	 * @param p the percentage of the data set to split off.
	 * @param REMOVE true if the split off portion should also be removed
	 * @return the split off portion of the data set
	 */
	public Data split(double p, boolean REMOVE)
	{
		Data data1 = new Data();	
		Data data2 = new Data();
		Random random = new Random();
		Vector<Pair<Double,Integer>> order = new Vector<Pair<Double,Integer>>();
		int index;
		int threshold = (int)Math.round(p*Y.size());
		
		//Assign random values to each index
		for(int i=0; i<Y.size(); i++){
			order.add(new Pair<Double,Integer>(random.nextDouble(), i));
		}
		
		Collections.sort(order);
		
		//Split data based on random values
		for(int i=0; i<order.size(); i++){
			index = order.get(i).second;
			
			if(i < threshold){
				data1.add(X.get(index), Y.get(index), W.get(index));
			}else{
				data2.add(X.get(index), Y.get(index), W.get(index));
			}
		}
		
		//Erase split off portion
		if(REMOVE){
			set(data2);
		}
		
		return data1;
	}
	
	/**
	 * Split off a portion of the data set.
	 * @param p the percentage of the data set to split off.
	 * @return the split off portion of the data set
	 */
	public Data split(double p)
	{
		return split(p, true);
	}
	
	/**
	 * Calculate the accuracy of the given classified results.
	 * @param Y1 obtained classifications for this data
	 * @return the accuracy
	 */
	public double accuracy(Vector<Double> Y1)
	{
		double result = 0;
		
		if(Y.size() == Y1.size()){			
			for(int i=0; i<Y.size(); i++){
				if(Y.get(i).equals(Y1.get(i))){
					result++;
				}
			}
			
			result /= Y.size();
		}
		
		return result;
	}
	
	/**
	 * Calculate the precision of the classified results given the ground truth.
	 * @param Y1 the classified labels
	 * @return the precision
	 */
	public double precision(Vector<Double> Y1)
	{
		double tp = 0, fp = 0, tn = 0, fn = 0;
		double tmpd = 0;
		
		if(Y.size() == Y1.size()){
			for(int i=0; i<Y.size(); i++){
				if(Y.get(i) > 0){
					if(Y1.get(i) > 0){
						tp++;
					}else if(Y1.get(i) < 0){
						fn++;
					}
				}else if(Y.get(i) < 0){
					if(Y1.get(i) > 0){
						fp++;
					}else if(Y1.get(i) < 0){
						tn++;
					}
				}
			}
			
			tmpd = tp / (tp + fp);
		}
		
		return tmpd;
	}

	/**
	 * Calculate the recall of the classified results given the ground truth.
	 * @param Y1 the classified labels
	 * @return the recall
	 */
	public double recall(Vector<Double> Y1)
	{
		double tp = 0, fp = 0, tn = 0, fn = 0;
		double tmpd = 0;
		
		if(Y.size() == Y1.size()){
			for(int i=0; i<Y.size(); i++){
				if(Y.get(i) > 0){
					if(Y1.get(i) > 0){
						tp++;
					}else if(Y1.get(i) < 0){
						fn++;
					}
				}else if(Y.get(i) < 0){
					if(Y1.get(i) > 0){
						fp++;
					}else if(Y1.get(i) < 0){
						tn++;
					}
				}
			}
			
			tmpd = tp / (tp + fn);
		}
		
		return tmpd;
	}

	/**
	 * Display the stored data.
	 */
	public void print()
	{
		for(int i=0; i<Y.size(); i++){
			System.out.print(Y.get(i) + ":");
			
			for(int j=0; j<X.get(i).length; j++){
				System.out.print(" " + X.get(i)[j]);
			}
			
			System.out.println();
		}
	}
	
	/**
	 * Print classifications.
	 * @param Y the classifications to print
	 */
	public static void print(Vector<Double> Y)
	{
		for(int i=0; i<Y.size(); i++){
			System.out.print(Y.get(i) + " ");
		}
		
		System.out.println();
	}
	
	/**
	 * Print binary results.
	 * @param Y the classifications to use
	 */
	public static void printClassDistribution(Vector<Double> Y)
	{
		double positives = 0;
		double negatives = 0;
		double neutrals = 0;
		
		for(int i=0; i<Y.size(); i++){
			if(Y.get(i) > 0) positives++;
			if(Y.get(i) < 0) negatives++;
			if(Y.get(i) == 0) neutrals++;
		}
		
		positives = 100 * (positives / Y.size());
		negatives = 100 * (negatives / Y.size());
		neutrals = 100 * (neutrals / Y.size());
		
		System.out.println("Size: " + Y.size());
		System.out.println("Positives: " + Utility.round(positives, 2));
		System.out.println("Negatives: " + Utility.round(negatives, 2));
		//System.out.println("Neutrals: " + Utility.round(neutrals, 20));
	}
	
	/**
	 * Extract the signs of the given doubles.
	 * @param Y a vector of doubles
	 * @return the signs of the given doubles
	 */
	public static Vector<Double> signs(Vector<Double> Y)
	{
		Vector<Double> Ysigns = new Vector<Double>();
		
		for(int i=0; i<Y.size(); i++){
			if(Y.get(i) > 0){
				Ysigns.add(1.0);
			}else if(Y.get(i) < 0){
				Ysigns.add(-1.0);
			}else{
				Ysigns.add(0.0);
			}
		}
		
		return Ysigns;
	}
	
	/**
	 * Sum up the values in a vector of doubles.
	 * @param v the vector to sum up
	 * @return the sum of the vectors values
	 */
	public static double sum(Vector<Double> v)
	{
		double tmpd = 0;
		
		for(int i=0; i<v.size(); i++){
			tmpd += v.get(i);
		}
		
		return tmpd;
	}
	
	/**
	 * Get the magnitued of a vector of doubles.
	 * @param v the vector
	 * @return the magnitude
	 */
	public static double magnitude(Vector<Double> v)
	{
		double tmpd = 0;
		
		for(int i=0; i<v.size(); i++){
			tmpd += v.get(i)*v.get(i);
		}
		
		tmpd = Math.sqrt(tmpd);
		
		return tmpd;
	}
	
	/**
	 * Add the values of the second vector to the values of the first.
	 * @param v1 the first vector
	 * @param v2 the second vector
	 */
	public static void plusEquals(Vector<Double> v1, Vector<Double> v2)
	{		
		if(v1.size() == v2.size()){			
			for(int i=0; i<v1.size(); i++){
				v1.set(i, v1.get(i)+v2.get(i));
			}
		}
	}
	
	/**
	 * Multiply the values of the first vector by the values of the second.
	 * @param v1 the first vector
	 * @param v2 the second vector
	 */
	public static void timesEquals(Vector<Double> v1, Vector<Double> v2)
	{		
		if(v1.size() == v2.size()){			
			for(int i=0; i<v1.size(); i++){
				v1.set(i, v1.get(i)*v2.get(i));
			}
		}
	}

	/**
	 * Divide the values of a vector of doubles.
	 * @param v the vector of doubles
	 * @param d the divisor
	 */
	public static void divideEquals(Vector<Double> v, double d)
	{
		for(int i=0; i<v.size(); i++){
			v.set(i, v.get(i)/d);
		}
	}
}