package kgm.ml;
import kgm.matrix.*;
import java.util.*;

/**
 * Cluster a vector of N-dimensional points into K groups.
 * @author Kenton McHenry
 */
public class KMeans
{
  /**
   * Cluster the given points.
   * @param points a vector N-dimensional points
   * @param k the number of clusters
   * @param iterations the number of iterations to re-cluster
   * @return a vector of vectors containing the point indices within the groups
   */
  public static Vector<Vector<Integer>> cluster(Vector<double[]> points, int k, int iterations)
  {    
  	Vector<double[]> cluster_centers = new Vector<double[]>();
  	Vector<Vector<double[]>> cluster_points = null;
  	Vector<Vector<Integer>> cluster_indices = null;
    Vector<double[]> minmax;
    int d;
    int minc;
    double mind, tmpd;
    
    if(points != null && points.size() > 0){
    	cluster_points = new Vector<Vector<double[]>>();
    	cluster_indices = new Vector<Vector<Integer>>();
      d = points.get(0).length;
      
      //Set random initial cluster centers
      minmax = MatrixUtility.minmax(points);
      
      for(int i=0; i<k; i++){
      	cluster_centers.add(MatrixUtility.random(d, minmax));
        cluster_points.add(new Vector<double[]>());
        cluster_indices.add(new Vector<Integer>());
      }
      
      //Cluster points
      for(int it=1; it<=iterations; it++){
        //System.out.println("Iteration: " + it);
        
        //Clear previous points
        for(int i=0; i<k; i++){
          cluster_points.get(i).clear();
          cluster_indices.get(i).clear();
        }
        
        //Assign points to nearest cluster
        for(int i=0; i<points.size(); i++){
          minc = 0;
          mind = Double.MAX_VALUE;
          
          for(int j=0; j<cluster_centers.size(); j++){
            tmpd = MatrixUtility.ssd(points.get(i), cluster_centers.get(j));
            
            if(tmpd < mind){
              minc = j;
              mind = tmpd;
            }
          }
          
          cluster_points.get(minc).add(points.get(i));
          cluster_indices.get(minc).add(i);
        }
        
        //Calculate new cluster centers
        for(int i=0; i<k; i++){
          cluster_centers.set(i, MatrixUtility.mean(cluster_points.get(i)));
          
          if(cluster_centers.get(i) == null){
            cluster_centers.set(i, MatrixUtility.random(d, minmax));
          }
        }
      }      
    }
    
    return cluster_indices;
  }
  
  /**
   * Get the cluster points.
   * @param points a vector of N-dimensional points
   * @param cluster_indices the point indices within each cluster
   * @return the points within each cluster
   */
  public static Vector<Vector<double[]>> getClusterPoints(Vector<double[]> points, Vector<Vector<Integer>> cluster_indices)
  {
  	Vector<Vector<double[]>> cluster_points = new Vector<Vector<double[]>>();
  	
  	for(int i=0; i<cluster_indices.size(); i++){
  		cluster_points.add(new Vector<double[]>());
  		
  		for(int j=0; j<cluster_indices.get(i).size(); j++){
  			cluster_points.lastElement().add(points.get(cluster_indices.get(i).get(j)));
  		}
  	}
  	
  	return cluster_points;
  }
  
  /**
   * Get cluster centers.
   * @param cluster_points the groups of points returned by clustering
   * @return the centers of each cluster
   */
  public static Vector<double[]> getClusterCenters(Vector<Vector<double[]>> cluster_points)
  {
  	Vector<double[]> cluster_centers = new Vector<double[]>();
    int k = cluster_points.size();
    int d = 0;
    double[] mean;
    
    //Determine d (in case some clusters have no points!)
    for(int i=0; i<cluster_points.size(); i++){
      if(cluster_points.get(i)!=null && cluster_points.get(i).size()>0){
        if(cluster_points.get(i).get(0).length > d){
          d = cluster_points.get(i).get(0).length;
        }
      }
    }
    
    //Set centers    
    for(int i=0; i<k; i++){
      mean = MatrixUtility.mean(cluster_points.get(i));
      
      if(mean == null){
        //System.out.println("Found null!");
        mean = MatrixUtility.vector(d, 0);
      }
      
      cluster_centers.add(mean);
    }
      
    return cluster_centers;
  }
    
  /**
   * Get cluster extremes.
   * @param cluster_points the groups of points returned by clustering
   * @return the extremes in each dimension of each cluster
   */
  public static Vector<Vector<double[]>> getClusterExtremes(Vector<Vector<double[]>> cluster_points)
  {
  	Vector<Vector<double[]>> cluster_extremes = new Vector<Vector<double[]>>();
    int k = cluster_points.size();

    for(int i=0; i<k; i++){
     cluster_extremes.add(MatrixUtility.minmax(cluster_points.get(i))); 
    }
    
    return cluster_extremes;
  }
  
  /**
   * Get cluster sizes.
   * @param cluster_points the groups of points returned by clustering
   * @return the number of points in each cluster
   */
  public static Vector<Integer> getClusterSizes(Vector<Vector<double[]>> cluster_points)
  {
  	Vector<Integer> cluster_sizes = new Vector<Integer>();
    int k = cluster_points.size();

    for(int i=0; i<k; i++){
      cluster_sizes.add(cluster_points.get(i).size());
    }
    
    return cluster_sizes;
  }
}