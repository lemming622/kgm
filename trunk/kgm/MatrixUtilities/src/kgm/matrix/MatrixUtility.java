package kgm.matrix;
import java.io.*;
import java.util.*;

/**
 * Utility functions used to manipulate matrices (many are matlab like).
 *  @author Kenton McHenry
 */
public class MatrixUtility
{
	/**
	 * An interface for an arbitrary distance function.
	 */
	public static interface Distance{
		public double evaluate(double[] x0, double[] x1);
	}
	
  /**
   * Convert a 1D array into a string so it can be printed to the standard output.
   *  @param arr the 1D array to convert
   *  @return the string representation of the array
   */
  public static String toString(double[] arr)
  {
  	if(arr == null) return "null";
  	
    String tmp = "";
    
    for(int i=0; i<arr.length; i++){
      tmp += arr[i] + " ";
    }
    
    return tmp;
  }
  
  /**
   * Convert a 2D array into a string so it can be printed to the standard output.
   *  @param mat the 2D array to convert
   *  @return the string representation of the array
   */
  public static String toString(double[][] mat)
  {
  	if(mat == null) return "null";
  	
    String tmp = "";
    
    for(int j=0; j<mat.length; j++){
      for(int i=0; i<mat[j].length; i++){
        tmp += mat[j][i] + " ";
      }
      
      tmp += "\n";
    }
    
    return tmp;
  }
  
  /**
   * Display the contents of a 1D array.
   *  @param arr the array to display
   */
  public static void print(double[] arr)
  {
    for(int i=0; i<arr.length; i++){
      System.out.print(arr[i] + " ");
    }
  }

	/**
   * Print a 1D array to standard output.
   *  @param arr the 1D array to print
   */
  public static void println(double[] arr)
  {
  	System.out.println(toString(arr));
  }
  
  /**
   * Print a 2D array to standard output.
   *  @param mat the 2D array to print
   */
  public static void println(double[][] mat)
  {
  	System.out.println(toString(mat));
  }
  
  /**
   * Collapse a 2D array into a 1D array (uses row major order).
   *  @param mat the 2D array
   *  @return the resulting 1D array
   */
  public static double[] to1D(double [][] mat)
  {
    int m = mat.length;
    int n = mat[0].length;
    double[] arr = new double[m*n];
    
    for(int j=0; j<m; j++){   
      for(int i=0; i<n; i++){
        arr[j*n+i] = mat[j][i];
      }
    }
    
    return arr;
  }
  
  /**
   * Convert a 1D array into a 2D array (uses row major order).
   *  @param m the number of rows the 2D version should have
   *  @param n the number of columns
   *  @param arr the 1D array
   *  @return the resulting 2D array
   */
  public static double[][] to2D(int m, int n, double[] arr)
  {
    double[][] mat = new double[m][n];
    
    for(int j=0; j<m; j++){    
      for(int i=0; i<n; i++){
        mat[j][i] = arr[j*n+i];
      }
    }
    
    return mat;
  }
  
  /**
   * Convert a 1D array into a 2D array (uses row major order).
   *  @param mat the 2D array to fill
   *  @param arr the 1D array
   */
  public static void to2D(double[][] mat, double[] arr)
  {
    int m = mat.length;
    int n = mat[0].length;
    
    for(int j=0; j<m; j++){    
      for(int i=0; i<n; i++){
        mat[j][i] = arr[j*n+i];
      }
    }
  }
  
  /**
   * Convert a vector of vectors into a 2D array.
   *  @param v the vector of vectors
   *  @return the resulting 2D array
   */
  public static double[][] to2D(Vector<Vector<Double>> v)
  {
    double[][] mat = new double[v.size()][v.get(0).size()];
    
    for(int j=0; j<v.size(); j++){
      for(int i=0; i<v.get(j).size(); i++){
        mat[j][i] = v.get(j).get(i);
      }
    }
    
    return mat;
  }
  
  /**
   * Copy values from one matrix over to the another without allocating new memory.
   *  @param T the target matrix
   *  @param S the source matrix
   */
  public static void set(double[][] T, double[][] S)
  {
  	for(int j=0; j<T.length; j++){
  		for(int i=0; i<T[j].length; i++){
  			T[j][i] = S[j][i];
  		}
  	}
  }

	/**
   * Copy an array.
   *  @param a the array to copy
   *  @return the duplicate array
   */
  public static double[] copy(double[] a)
  {
    int n = a.length;
    double[] b = new double[n];
    
    for(int i=0; i<n; i++){
      b[i] = a[i];
    }
    
    return b;
  }
  
	/**
   * Copy a matrix.
   *  @param A the matrix to copy
   *  @return the duplicate matrix
   */
  public static double[][] copy(double[][] A)
  {
    int m = A.length;
    int n = A[0].length;
    double[][] B = new double[m][n];
    
    for(int j=0; j<m; j++){
      for(int i=0; i<n; i++){
        B[j][i] = A[j][i];
      }
    }
    
    return B;
  }

	/**
   * Create a new array of the specified length with the given initial value.
   *  @param d the size of the array
   *  @param v0 the initial value
   *  @return the created array
   */
  public static double[] vector(int d, double v0)
  {
    double[] tmpv = new double[d];
    
    for(int i=0; i<d; i++){
      tmpv[i] = v0;
    }
    
    return tmpv;
  }

	/**
   * Get a sub-array given based on a set of indices.
   * @param a the array to retrieve elements from
   * @param indices the desired indices
   * @return the sub-array
   */
  public static double[] vector(double[] a, double[] indices)
  {
  	double[] b = new double[indices.length];
  	
  	for(int i=0; i<indices.length; i++){
  		b[i] = a[(int)indices[i]];
  	}
  	
  	return b;
  }

	/**
	 * Create an array from a vector of values.
	 * @param v a vector of values
	 * @return the created array
	 */
	public static double[] vector(Vector v)
	{
		double[] a = new double[v.size()];
		
		for(int i=0; i<v.size(); i++){
			if(v.get(i) instanceof Double){
				a[i] = (Double)v.get(i);
			}else if(v.get(i) instanceof Integer){
				a[i] = (Integer)v.get(i);
			}
		}
		
		return a;
	}

	/**
   * Create an mxn matrix filled with the given value.
   * @param m the number of rows
   * @param n the number of columns
   * @param a the value to fill the matrix with
   * @return a matrix filled with the given value
   */
  public static double[][] matrix(int m, int n, double a)
  {  	
  	double[][] A = new double[m][n];
  	
  	for(int j=0; j<m; j++){
  		for(int i=0; i<n; i++){
  			A[j][i] = a;
  		}
  	}
  	
  	return A;
  }

	/**
   * Create an mxn matrix with all zero values.
   * @param m the number of rows
   * @param n the number of columns
   * @return the matrix of zeros
   */
  public static double[][] zeros(int m, int n)
  {  	
  	return new double[m][n];
  }
  
  /**
   * Create an mxn matrix with all one values.
   * @param m the number of rows
   * @param n the number of columns
   * @return the matrix of ones
   */
  public static double[][] ones(int m, int n)
  {  	
  	double[][] A = new double[m][n];
  	
  	for(int j=0; j<m; j++){
  		for(int i=0; i<n; i++){
  			A[j][i] = 1;
  		}
  	}
  	
  	return A;
  }
  
  /**
   * Create an nxn identity matrix.
   *  @param n the number of rows/columns
   *  @return the resulting identity matrix (as a 2D array)
   */
  public static double[][] eye(int n)
  {
    double[][] M = new double[n][n];
    
    for(int i=0; i<n; i++){
      for(int j=0; j<n; j++){
        if(i == j){
          M[j][i] = 1;
        }else{
          M[j][i] = 0;
        }
      }
    }
    
    return M;
  }

	/**
   * Create an array of the specified size with random values.
   *  @param d the size of the array
   *  @param minmax the range from which the random values should lay
   *  @return the resulting array
   */
  public static double[] random(int d, Vector<double[]> minmax)
  {
    double[] tmpv = vector(d, 0);
    double[] min = minmax.get(0);
    double[] max = minmax.get(1);
    double a, b;
    
    for(int i=0; i<d; i++){
      a = max[i] - min[i];
      b = min[i];
      tmpv[i] = a*Math.random() + b;
    }
    
    return tmpv;
  }
  
	/**
   * Resize an array using interpolation.
   * @param x the image
   * @param l_new the desired array length
   */
  public static double[] resizeAndInterpolate(double[] x, int l_new)
  {    
  	double[] x_new = new double[l_new];
  	double scale = ((double)(l_new-1)) / ((double)(x.length-1));
  	double i_old;
    int i_old_whole;
    double i_old_remainder;
    double A, B;
  
    //Pad so interpolation is always correct
    double[] x_padded = new double[x.length+1];
    
    for(int i=0; i<x.length; i++){
      x_padded[i] = x[i];
    }
    
    x_padded[x.length] = x[x.length-1];
    
    //Derive resized array  
    for(int i=0; i<l_new; i++){
      i_old = ((double)i) / scale;
      i_old_whole = (int)Math.round(i_old);
      i_old_remainder = i_old - (double)i_old_whole;

      A = x_padded[i_old_whole];
      B = x_padded[i_old_whole+1];

      x_new[i] = (1.0-i_old_remainder)*A + i_old_remainder*B;
    }
  
    return x_new;
  }
  
  /**
   * Subscripted reference.
   * @param a an array
   * @param indices the indices
   * @return an array containing the subscripted elements
   */
  public static double[] subsref(double[] a, double[] indices)
  {
  	double[] b = new double[indices.length];
  	
		for(int i=0; i<indices.length; i++){
			b[i] = a[(int)indices[i]];
		}
  	
  	return b;
  }
  
  /**
   * Retrieve a submatrix from a larger matrix.
   * @param A a matrix
   * @param j0 the starting vertical index
   * @param j1 the ending vertical index
   * @param i0 the starting horizontal index
   * @param i1 the ending horizontal index
   * @return the submatrix
   */
  public static double[][] subsref(double[][] A, int j0, int j1, int i0, int i1)
  {
  	int m = j1 - j0 + 1;
  	int n = i1 - i0 + 1;
  	double[][] B = new double[m][n];
  	
  	for(int j=0; j<m; j++){
  		for(int i=0; i<n; i++){
  			B[j][i] = A[j+j0][i+i0];
  		}
  	}
  	
  	return B;
  }
  
  /**
	 * Subscripted assignment.
	 * @param a an array
	 * @param indices the indices of a to assign
	 * @param value the value to assign
	 * @return a copy of a with value assign according to the given indices
	 */
	public static double[] subsasgn(double[] a, double[] indices, double value)
	{
		double[] b = copy(a);
		
		for(int i=0; i<indices.length; i++){
			b[(int)indices[i]] = value;
		}
		
		return b;
	}

	/**
   * Subscripted assignment.
   * @param a an array
   * @param indices the indices of a to assign
   * @param b the array whose values will be assigned to a
   * @return a copy of a with assignments from b
   */
  public static double[] subsasgn(double[] a, double[] indices, double[] b)
  {
  	double[] c = copy(a);
  	
  	for(int i=0; i<indices.length; i++){
  		c[(int)indices[i]] = b[i];
  	}
  	
  	return c;
  }
  
  /**
   * Assign a matrix as a submatrix to the given matrix.
   * @param A a matrix
   * @param j0 the starting vertical location of the submatrix
   * @param i0 the starting horizontal location of the submatrix
   * @param B the submatrix
   * @return the combined matrix
   */
  public static double[][] subsasgn(double[][] A, int j0, int i0, double[][] B)
  {
  	double[][] C = copy(A);
  	int m = B.length;
  	int n = B[0].length;
  	
  	for(int j=0; j<m; j++){
  		for(int i=0; i<n; i++){
  			C[j+j0][i+i0] = B[j][i];
  		}
  	}
  	
  	return C;
  }
  
  /**
   * Find indices of non-zero elements.
   * @param a the array to search for non-zero elements
   * @return an array of indices of non-zero elements
   */
  public static double[] find(double[] a)
  {
  	Vector<Integer> tmpv = new Vector<Integer>();
  	double[] indices = null;
  	int n = a.length;
  	
  	for(int i=0; i<n; i++){
  		if(a[i] > 0){
  			tmpv.add(i);
  		}
  	}
  	
  	indices = new double[tmpv.size()];
  	
  	for(int i=0; i<tmpv.size(); i++){
  	  indices[i] = tmpv.get(i);	
  	}
  	
  	return indices;
  }
  
  /**
   * Retrieve the unique values within a given array.
   * @param a an array
   * @return the unique values within the given array 
   */
  public static double[] unique(double[] a)
  {
  	TreeSet<Double> set = new TreeSet<Double>();
  	Iterator<Double> itr;
  	double[] b = null;
  	int at;
  	
  	for(int i=0; i<a.length; i++){
  	  if(!set.contains(a[i])) set.add(a[i]);	
  	}
  	
  	b = new double[set.size()];
  	itr = set.iterator();
  	at = 0;
  	
  	while(itr.hasNext()){
  		b[at] = itr.next();
  		at++;
  	}
  	
  	return b;
  }
  
  /**
   * True for set members.
   * @param a an array
   * @param b a set
   * @return an array with a 1 at indices of a that were members of b
   */
  public static double[] ismember(double[] a, double[] b)
  {
  	int n = a.length;
  	double[] c = new double[n];
  	TreeSet<Double> set = new TreeSet<Double>();
  	
  	for(int i=0; i<b.length; i++){
  		set.add(b[i]);
  	}
  	
  	for(int i=0; i<a.length; i++){
  		if(set.contains(a[i])){
  			c[i] = 1;
  		}
  	}
  	
  	return c;
  }
  
  /**
   * Convert a non-homogeneous 3x3 transformation to a 4x4 homogeneous transformation.
   *  @param A the non-homogeneous rotation
   *  @return the homogeneous transformation
   */
  public static double[][] homogeneousRotation(double[][] A)
  {
  	double[][] B = eye(4);
  	
  	for(int j=0; j<3; j++){
  		for(int i=0; i<3; i++){
  			B[j][i] = A[j][i];
  		}
  	}
  	
  	return B;
  }
  
	/**
   * Construct a 2D scale matrix.
   * @param sx the scale in the x direction
   * @param sy the scale in the y direction
   * @return the scale matrix
   */
  public static double[][] scale(double sx, double sy)
  {
  	double[][] S = eye(3);
  	
  	S[0][0] = sx;
  	S[1][1] = sy;
  	
  	return S;
  }
  
	/**
   * Construct a 2D translation matrix.
   * @param tx the translation in the x direction
   * @param ty the translation in the y direction
   * @return the translation matrix
   */
  public static double[][] translate(double tx, double ty)
  {
  	double[][] T = eye(3);
  	
  	T[0][2] = tx;
  	T[1][2] = ty;
  	
  	return T;
  }
  
	/**
   * Construct a 2D rotation matrix.
   * @param theta the angle in degrees to rotate
   * @return the rotation matrix
   */
  public static double[][] rotate(double theta)
  {
  	double[][] R = eye(3);
  	
  	theta = Math.PI * theta / 180.0;
  
  	R[0][0] = Math.cos(theta);
  	R[0][1] = -Math.sin(theta);
  	R[1][0] = Math.sin(theta);
  	R[1][1] = Math.cos(theta);
  	
  	return R;
  }
  
	/**
   * Construct a rotation matrix about the x-axis.
   *  @param theta the angle in degrees around the x-axis
   *  @return the rotation matrix about the x-axis
   */
  public static double[][] rotateX(double theta)
  {
  	double[][] R = eye(4);
  	
  	theta = Math.PI * theta / 180.0;
  
  	R[1][1] = Math.cos(theta);
  	R[1][2] = -Math.sin(theta);
  	R[2][1] = Math.sin(theta);
  	R[2][2] = Math.cos(theta);
  	
  	return R;
  }

	/**
   * Construct a rotation matrix about the y-axis.
   *  @param theta the angle in degrees around the y-axis
   *  @return the rotation matrix about the y-axis
   */
  public static double[][] rotateY(double theta)
  {
  	double[][] R = eye(4);
  	
  	theta = Math.PI * theta / 180.0;
  
  	R[0][0] = Math.cos(theta);
  	R[0][2] = Math.sin(theta);
  	R[2][0] = -Math.sin(theta);
  	R[2][2] = Math.cos(theta);
  	
  	return R;
  }

	/**
   * Construct a rotation matrix about the z-axis.
   *  @param theta the angle in degrees around the z-axis
   *  @return the rotation matrix about the z-axis
   */
  public static double[][] rotateZ(double theta)
  {
  	double[][] R = eye(4);
  	
  	theta = Math.PI * theta / 180.0;
  	
  	R[0][0] = Math.cos(theta);
  	R[0][1] = -Math.sin(theta);
  	R[1][0] = Math.sin(theta);
  	R[1][1] = Math.cos(theta);
  	
  	return R;
  }

	/**
   * Construct a rotation matrix given an angle about a normalized axis.
   *  @param theta the angle in degress around the axis
   *  @param x the x component of the axis
   *  @param y the y component of the axis
   *  @param z the z component of the axis
   *  @return the rotation matrix
   */
  public static double[][] rotate(double theta, double x, double y, double z)
  {
  	double[][] R = eye(4);
  	double c, s;
  	
  	theta = Math.PI * theta / 180.0;
    c = Math.cos(theta);
    s = Math.sin(theta);
  	
  	R[0][0] = x*x+(1-x*x)*c;
  	R[0][1] = x*y*(1-c)-z*s;
  	R[0][2] = x*z*(1-c)+y*s;
  	R[1][0] = x*y*(1-c)+z*s;
  	R[1][1] = y*y+(1-y*y)*c;
  	R[1][2] = y*z*(1-c)-x*s;
  	R[2][0] = x*z*(1-c)-y*s;
  	R[2][1] = y*z*(1-c)+x*s;
  	R[2][2] = z*z+(1-z*z)*c;
  	
  	return R;
  }

	/**
   * Rotate about Z, Y, then X axis given the given angles in degrees.
   *  @param rx the rotation in degrees around the x-axis
   *  @param ry the rotation in degrees around the y-axis
   *  @param rz the rotation in degrees around the z-axis
   *  @return the overall rotation matrix
   */
  public static double[][] rotateXYZ(double rx, double ry, double rz)
  {
  	/*
  	double[][] Rx = MatrixUtils.rotateX(rx);
  	double[][] Ry = MatrixUtils.rotateY(ry);
  	double[][] Rz = MatrixUtils.rotateZ(rz);
  	*/
  	  	  	
  	double[][] Rx = MatrixUtility.rotate(rx, 1, 0, 0);
  	double[][] Ry = MatrixUtility.rotate(ry, 0, 1, 0);
  	double[][] Rz = MatrixUtility.rotate(rz, 0, 0, 1);
  	
  	double[][] R = MatrixUtility.mtimes(MatrixUtility.mtimes(Rx, Ry), Rz);
  	
  	return R;
  }

	/**
   * Construct a 4x4 rotation matrix that rotates the coordinate space represented by the vectors (ib, jb, kb)
   * to the coordinate space represented by the vectors (ia, ja, ka).
   *  @param ia the x axis in coordinate space a
   *  @param ja the y axis in coordinate space a
   *  @param ka the z axis in coordinate space a
   *  @param ib the x axis in coordinate space b
   *  @param jb the y axis in coordinate space b
   *  @param kb the z axis in coordinate space b
   *  @return the resulting rotation matrix
   */
  public static double[][] rotate(double[] ia, double[] ja, double[] ka, double[] ib, double[] jb, double[] kb)
  {
    double[][] M = new double[4][4];
    
    M[0][0] = dot(ia,ib);
    M[0][1] = dot(ja,ib);
    M[0][2] = dot(ka,ib);
    M[0][3] = 0;
    
    M[1][0] = dot(ia,jb);
    M[1][1] = dot(ja,jb);
    M[1][2] = dot(ka,jb);
    M[1][3] = 0;
    
    M[2][0] = dot(ia,kb);
    M[2][1] = dot(ja,kb);
    M[2][2] = dot(ka,kb);
    M[2][3] = 0;
    
    M[3][0] = 0;
    M[3][1] = 0;
    M[3][2] = 0;
    M[3][3] = 1;
    
    return M;
  }

	/**
   * Construct a 4x4 rotation matrix that rotates the coordinate space represented by the vectors (ib, jb, kb)
   * to the coordinate space represented by the vectors (ia, ja, ka), where ia = (1, 0, 0), 
   * ja = (0, 1, 0), and ka = (0, 0, 1)
   *  @param ib the x axis in coordinate space a
   *  @param jb the y axis in coordinate space a
   *  @param kb the z axis in coordinate space a
   *  @return the resulting rotation matrix
   */
  public static double[][] rotate(double[] ib, double[] jb, double[] kb)
  {
    double[][] M = new double[4][4];
    
    M[0][0] = ib[0];
    M[0][1] = ib[1];
    M[0][2] = ib[2];
    M[0][3] = 0;
    
    M[1][0] = jb[0];
    M[1][1] = jb[1];
    M[1][2] = jb[2];
    M[1][3] = 0;
    
    M[2][0] = kb[0];
    M[2][1] = kb[1];
    M[2][2] = kb[2];
    M[2][3] = 0;
    
    M[3][0] = 0;
    M[3][1] = 0;
    M[3][2] = 0;
    M[3][3] = 1;
    
    return M;
  }
  
  /**
   * Construct a 4x4 rotation matrix from the given forward and up directions of a camera
   * @param forward the forward looking vector
   * @param up an upward pointing vector
   * @return a rotation orienting a canonical camera to the given directions
   */
  public static double[][] rotateCamera(double[] forward, double[] up)
  {
  	double[] right;
  	
  	forward = divide(forward, norm(forward));
  	up = divide(up, norm(up));
  	right = cross(forward, up);
  	right = divide(right, norm(right));
  	
  	return rotateCamera(forward, up, right);
  }

  /**
   * Construct a 4x4 rotation matrix from the given forward, up, and side directions of a camera
   * @param forward the forward looking vector
   * @param up an upward pointing vector
   * @param right the side direction vector
   * @return a rotation orienting a canonical camera to the given directions
   */
  public static double[][] rotateCamera(double[] forward, double[] up, double[] right)
  {
  	double[][] M = new double[4][4];
  	double[] tmp;
  	
  	forward = divide(forward, norm(forward));
  	up = divide(up, norm(up));
  	right = divide(right, norm(right));
  	
  	tmp = cross(right, forward);
  	tmp = divide(tmp, norm(tmp));
  	
    M[0][0] = right[0];
    M[0][1] = right[1];
    M[0][2] = right[2];
    M[0][3] = 0;
    
    M[1][0] = tmp[0];
    M[1][1] = tmp[1];
    M[1][2] = tmp[2];
    M[1][3] = 0;
    
    M[2][0] = -forward[0];
    M[2][1] = -forward[1];
    M[2][2] = -forward[2];
    M[2][3] = 0;
    
    M[3][0] = 0;
    M[3][1] = 0;
    M[3][2] = 0;
    M[3][3] = 1;
  	
  	return M;
  }
  
	/**
   * Transpose a matrix.
   *  @param A the matrix to transpose
   *  @return the transposed matrix
   */
  public static double[][] transpose(double[][] A)
  {
  	if(A == null) return null;
  	
  	int m = A.length;
  	int n = A[0].length;
  	double[][] B = new double[n][m];
  	
  	for(int j=0; j<m; j++){
  		for(int i=0; i<n; i++){
  			B[i][j] = A[j][i];
  		}
  	}
  	
  	return B;
  }
  
  /**
   * Create a skew symmetric matrix for the given 3-vector that can be used for computing a cross product with
   * another vector.
   * @param a a 3-vector
   * @return the vectors skew matrix
   */
  public static double[][] skew(double[] a)
  {
  	double[][] B = null;
  	
  	if(a.length == 3){
  		B = new double[3][3];
	  	B[0][1] = -a[2];
	  	B[0][2] = a[1];
	  	B[1][0] = a[2];
	  	B[1][2] = -a[0];
	  	B[2][0] = -a[1];
	  	B[2][1] = a[0];
  	}
  	
  	return B;
  }

  /**
   * Calculate the sum of the given array.
   * @param a an array
   * @return the sum of the array
   */
  public static double sum(double[] a)
  {
  	double sum = 0;
  	
  	for(int i=0; i<a.length; i++){
  		sum += a[i];
  	}
  	
  	return sum;
  }
  
	/**
   * Calculate the mean of an array of values.
   * @param a the array
   * @return the mean of the given array
   */
  public static double mean(double[] a)
  {
    double sum = 0;
    
    for(int i=0; i<a.length; i++){
    	sum += a[i];
    }
    
    return sum / a.length;
  }
  
  /**
   * Calculate the mean of the given vector of points.
   *  @param X the vector of points
   *  @return the mean
   */
  public static double[] mean(Vector<double[]> X)
  {
    double[] M = null;
    
    if(X != null && X.size() > 0){
      int d = X.get(0).length;
      M = vector(d, 0);
      
      for(int i=0; i<X.size(); i++){
        for(int j=0; j<d; j++){
          M[j] += X.get(i)[j];
        }
      }
      
      for(int i=0; i<d; i++){
        M[i] /= (double)X.size();
      }
    }
    
    return M;
  }
  
  /**
   * Calculate the standard deviation of the given array of values.
   * @param a an array
   * @param mean the mean of the array
   * @return the standard deviation
   */
  public static double std(double[] a, double mean)
  {
  	double std = 0;
  	double tmpd;
  	
  	for(int i=0; i<a.length; i++){
  		tmpd = a[i] - mean;
  		std += tmpd * tmpd;
  	}
  	
  	if(a.length > 0) std /= a.length;
  	std = Math.sqrt(std);
  	
  	return std;
  }
  
  /**
   * Return the an array with the absolute values of each element.
   * @param a an array
   * @return an array containing the absolute values
   */
  public static double[] abs(double[] a)
  {
  	double[] b = new double[a.length];
  	
  	for(int i=0; i<a.length; i++){
  		b[i] = Math.abs(a[i]);
  	}
  	
  	return b;
  }
  
  /**
   * Calculate the maximum values between two arrays.
   * @param a the first array
   * @param b the second array
   * @return the element-wise max of the two arrays
   */
  public static double[] max(double[] a, double[] b)
  {
  	int n = a.length;
  	double[] c = new double[n];
  	
  	for(int i=0; i<n; i++){
  		if(a[i] > b[i]){
  			c[i] = a[i];
  		}else{
  			c[i] = b[i];
  		}
  	}
  	
  	return c;
  }
  
  /**
   * Get the maximum value within the given array.
   * @param a an array
   * @return the maximum value within the array
   */
  public static double max(double[] a)
  {
  	double max = -Double.MAX_VALUE;
  	
		for(int i=0; i<a.length; i++){
			if(a[i] > max){
				max = a[i];
			}
		}
  	
  	return max;
  }
  
  /**
   * Get the index of the maximum value within the given array.
   * @param a an array
   * @return the index of the maximum value within the array
   */
  public static int max_index(double[] a)
  {
  	double max = -Double.MAX_VALUE;
  	int maxi = -1;
  	
		for(int i=0; i<a.length; i++){
			if(a[i] > max){
				max = a[i];
				maxi = i;
			}
		}
  	
  	return maxi;
  }
  
  /**
   * Get the maximum value within the given matrix.
   * @param A a matrix
   * @return the maximum value within the matrix
   */
  public static double max(double[][] A)
  {
  	int m = A.length;
  	int n = A[0].length;
  	double max = -Double.MAX_VALUE;
  	
  	for(int j=0; j<m; j++){
  		for(int i=0; i<n; i++){
  			if(A[j][i] > max){
  				max = A[j][i];
  			}
  		}
  	}
  	
  	return max;
  }
  
  /**
   * Logical OR operator.
   * @param a an array
   * @param b an array
   * @return an array containing a 1 at index i if a[i] or b[i] is greater than 0
   */
  public static double[] or(double[] a, double[] b)
  {
  	int n = a.length;
  	double[] c = new double[n];
  	
  	for(int i=0; i<n; i++){
  		if(a[i]>0 || b[i]>0) c[i] = 1;
  	}
  	
  	return c;
  }
  
  /**
   * Compute the extremes of each component in the given vector of points.
   *  @param X the vector of points
   *  @return a vector containing a pair of points (the minimal point and the maximal point)
   */
  public static Vector<double[]> minmax(Vector<double[]> X)
  {
    Vector<double[]> tmpv = new Vector<double[]>();
    
    if(X != null && X.size() > 1){
      int d = X.get(0).length;
      double[] min = vector(d, Double.MAX_VALUE);
      double[] max = vector(d, Double.MIN_VALUE);
      
      for(int i=0; i<X.size(); i++){
        for(int j=0; j<d; j++){
          if(X.get(i)[j] < min[j]) min[j] = X.get(i)[j];
          if(X.get(i)[j] > max[j]) max[j] = X.get(i)[j];
        }
      }
      
      tmpv.add(min);
      tmpv.add(max);
    }
    
    return tmpv;
  }
  
  /**
   * The sum of squared difference between two points.
   *  @param X0 the first point
   *  @param X1 the second point
   *  @return the SSD between the two points
   */
  public static double ssd(double[] X0, double[] X1)
  {
    double sum = 0;
    double tmp;
    
    for(int i=0; i<X0.length; i++){
      tmp = X0[i];     
      tmp -= X1[i];
      tmp *= tmp;
      sum += tmp;
    }
    
    return sum;
  }
  
  /**
   * The Euclidean distance between two points.
   *  @param X0 the first point
   *  @param X1 the second point
   *  @return the distance between the two points
   */
  public static double distance(double[] X0, double[]X1)
  {
    return Math.sqrt(ssd(X0, X1));
  }
  
  /**
   * Dynamic Time Warping
   * @param X0 a collection of feature vectors
   * @param X1 another collection of feature vectors
   * @param distance the distance to use when comparing elements
   * @return a distance between the two
   */
  public static double dtw(Vector<double[]> X0, Vector<double[]> X1, Distance distance)
  {
  	double[][] d = new double[X0.size()+1][X1.size()+1];
		double d0, d1, d2, tmpd;
		
		//Initialize
		d[0][0] = 0;
		
		for(int i=1; i<=X0.size(); i++){
			d[i][0] = Double.MAX_VALUE;
		}
		
		for(int j=1; j<=X1.size(); j++){
			d[0][j] = Double.MAX_VALUE;
		}
		
		//Fill in table
		for(int i=1; i<=X0.size(); i++){
			for(int j=1; j<=X1.size(); j++){
				tmpd = distance.evaluate(X0.get(i-1), X1.get(j-1));
				d0 = d[i-1][j];			//Deletion
				d1 = d[i][j-1];			//Insertion
				d2 = d[i-1][j-1];		//Match/Substitution
				
				d[i][j] = tmpd + Math.min(Math.min(d0, d1), d2);
			}
		}
		
		return d[X0.size()][X1.size()];
  }
  
  /**
   * Dynamic Time Warping
   * @param X0 a collection of feature vectors
   * @param X1 another collection of feature vectors
   * @return a distance between the two
   */
  public static double dtw(Vector<double[]> X0, Vector<double[]> X1)
  {
  	return dtw(X0, X1, new Distance(){
  		public double evaluate(double[] x0, double[] x1){
  			return distance(x0, x1);
  		}
  	});
  }
  
  /**
	 * Compute the cosine similarity between two points.
	 * @param X0 a point
	 * @param X1 a point
	 * @return the cosine similarity
	 */
	public static double cos(double[] X0, double[] X1)
	{
		return dot(X0, X1) / (norm(X0) * norm(X1));
	}

	/**
   * The discrete cosine transform.
   * @param f a function
   * @return an array of coefficients
   */
  public static double[] dct(double[] f)
  {
  	int n = f.length;
  	double[] c = new double[n];
  	
  	for(int i=0; i<n; i++){
  		for(int x=0; x<n; x++){
  			c[i] += f[x]*Math.cos((Math.PI*(2*x+1)*i)/(2*n));
  		}
  		
  		if(i == 0){
  			c[i] *= Math.sqrt(1.0/n);
  		}else{
  			c[i] *= Math.sqrt(2.0/n);
  		}
  	}
  	
  	return c;
  }
  
  /**
   * The inverse discrete cosine transform.
   * @param c an array of coefficients
   * @return the function
   */
  public static double[] idct(double[] c)
  {
  	int n = c.length;
  	double[] f = new double[n];
  	
  	for(int x=0; x<n; x++){
  		for(int i=0; i<n; i++){
  			if(i == 0){
  				f[x] += Math.sqrt(1.0/n)*c[i]*Math.cos((Math.PI*(2*x+1)*i)/(2*n));
  			}else{
  				f[x] += Math.sqrt(2.0/n)*c[i]*Math.cos((Math.PI*(2*x+1)*i)/(2*n));
  			}
  		}
  	}
  	
  	return f;
  }
  
	/**
   * The discrete Fourier transform.
   * @param x a function
   * @return an array of coefficients
   */
  public static double[] dft(double[] x)
  {
  	int n = x.length;
  	double[] f = new double[n];
  	double c, s;
  	
  	for(int i=0; i<n; i++){
  		c = 0;
  		s = 0;
  		
  		for(int j=0; j<n; j++){
  			c += x[j]*Math.cos((2*Math.PI*i*j)/n);
  			s += x[j]*Math.sin((2*Math.PI*i*j)/n);
  		}
  		
  		f[i] = Math.sqrt(c*c + s*s);
  	}
  	
  	return f;
  }
  
  /**
   * Calculate the norm of the given point.
   *  @param X the point
   *  @return the norm/magnitude of the point
   */
  public static double norm(double[] X)
  {
    double sum = 0;
    
    for(int i=0; i<X.length; i++){
      sum += X[i]*X[i];
    }
    
    return Math.sqrt(sum);
  }
  
  /**
	 * Create a 1D gaussian filter.
	 * @param n the length of the filter
	 * @param sigma the standard deviation of gaussian
	 * @return the filter
	 */
	public static double[] fgaussian(int n, double sigma)
	{
		double[] f = new double[n];
		double sigma_sqrd = sigma * sigma;
		double C = 1.0 / (2.0 * Math.PI * sigma_sqrd);
		double center = n/2.0 - 0.5;
		double itmp;
	
		for(int i=0; i<n; i++){
			itmp = i - center;
			f[i] = C * Math.exp(-0.5*itmp*itmp/sigma_sqrd);
		}
	
		f = divide(f, norm(f));
		
		return f;
	}

	/**
	 * Convolve the given vector with the given filter.
	 * @param a a vector
	 * @param f a filter
	 * @return the convolved vector (of the same size!)
	 */
	public static double[] conv(double[] a, double[] f)
	{  	
		double[] b = new double[a.length];
		int fradius = (f.length%2 == 1) ? (f.length-1)/2 : f.length/2;
		double sum = 0;
	
		for(int i=fradius; i<a.length-fradius; i++){
			sum = 0;
	
			for(int j=0; j<f.length; j++){
				sum += f[j] * a[j+i-fradius];
			}
	
			b[i] = sum;
		}
	
		return b;
	}

	/**
   * Perform an element-wise addition to the given array.
   *  @param a an array of doubles
   *  @param b the value to add on
   *  @return the sum
   */
  public static double[] plus(double[] a, double b)
  {
    double[] c = new double[a.length];
    
    for(int i=0; i<a.length; i++){
      c[i] = a[i] + b;
    }
    
    return c;
  }
  
  /**
   * Perform an element-wise addition of the given matrices.
   *  @param A a matrix
   *  @param B a matrix
   *  @return the sum of the two matrices
   */
  public static double[][] plus(double[][] A, double[][] B)
  {
  	int m = A.length;
  	int n = A[0].length;
    double[][] C = new double[m][n];
    
    for(int j=0; j<m; j++){
	    for(int i=0; i<n; i++){
	      C[j][i] = A[j][i] + B[j][i];
	    }
    }
    
    return C;
  }
  
  /**
   * Negate an array's contents.
   * @param a the array to negate
   * @return the negated array
   */
  public static double[] uminus(double[] a)
  {
  	int n = a.length;
  	double[] b = new double[n];
  	
  	for(int i=0; i<n; i++){
  		b[i] = -a[i];
  	}
  	
  	return b;
  }
  
  /**
   * Negate a matrices contents.
   * @param A the matrix to negate
   * @return the negated matrix
   */
  public static double[][] uminus(double[][] A)
  {
  	int m = A.length;
  	int n = A[0].length;
  	double[][] B = new double[m][n];
  	
  	for(int j=0; j<m; j++){
  		for(int i=0; i<n; i++){
  			B[j][i] = -A[j][i];
  		}
  	}
  	
  	return B;
  }  
  
  /**
   * Perform an element-wise subtraction of the two arrays.
   *  @param a an array of doubles
   *  @param b another array of doubles to subtract off the first array
   *  @return the difference
   */
  public static double[] minus(double[] a, double[] b)
  {
    double[] c = new double[a.length];
    
    for(int i=0; i<a.length; i++){
      c[i] = a[i] - b[i];
    }
    
    return c;
  }
  
  /**
   * Perform an element-wise subtraction of the two matrices.
   *  @param A a matrix
   *  @param B another matrix to subtract off the first matrix
   *  @return the difference
   */
  public static double[][] minus(double[][] A, double[][] B)
  {
  	int m = A.length;
  	int n = A[0].length;
    double[][] C = new double[m][n];
    
    for(int j=0; j<m; j++){
	    for(int i=0; i<n; i++){
	      C[j][i] = A[j][i] - B[j][i];
	    }
    }
    
    return C;
  }
  
  /**
   * Multiply an array element wise by the given value.
   * @param A the array to multiply
   * @param d the factor
   * @return the resulting array
   */
  public static double[] times(double[] A, double d)
  {
  	int n = A.length;
  	double[] B = new double[n];
  	
  	for(int i=0; i<n; i++){
  		B[i] = A[i] * d;
  	}
  	
  	return B;
  }

  /**
   * Multiply an array element wise by the given array of factors.
   * @param a the array to multiply
   * @param b the array of factors
   * @return the resulting array
   */
  public static double[] times(double[] a, double[] b)
  {
  	int n = a.length;
  	double[] c = new double[n];
  	
  	for(int i=0; i<n; i++){
  		c[i] = a[i] * b[i];
  	}
  	
  	return c;
  }
  
  /**
   * Multiply a matrix element wise by the given value.
   * @param A the matrix to multiply
   * @param d the factor
   * @return the resulting matrix
   */
  public static double[][] times(double[][] A, double d)
  {
  	int m = A.length;
  	int n = A[0].length;
  	double[][] B = new double[m][n];
  	
  	for(int j=0; j<m; j++){
	  	for(int i=0; i<n; i++){
	  		B[j][i] = A[j][i] * d;
	  	}
  	}
  	
  	return B;
  }
  
  /**
   * Perform a matrix multiplication.
   *  @param A a matrix multiplicand
   *  @param B a matrix multiplicand
   *  @return the result of the multiplication
   */
  public static double[][] mtimes(double[][] A, double[][] B)
  {
  	if(A == null || B == null) return null;
  	
  	int m = A.length;
  	int xa = A[0].length;
  	int xb = B.length;
  	int n = B[0].length;
  	double tmpd;
  	
  	if(xa != xb) return null;
  	
  	double[][] C = new double[m][n];
  	
  	for(int j=0; j<m; j++){
  	  for(int i=0; i<n; i++){
  	  	tmpd = 0;
  	  	
  	  	for(int k=0; k<xa; k++){
  	  		tmpd += A[j][k] * B[k][i];
  	  	}
  	  	
  	  	C[j][i] = tmpd;
  	  }
  	}
  	
  	return C;
  }
  
  /**
   * Perform a matrix multiplication.
   *  @param A a matrix multiplicand
   *  @param b an array multiplicand (assumed to be nx1)
   *  @return the result of the multiplication
   */
  public static double[] mtimes(double[][] A, double[] b)
  {
  	if(A == null || b == null) return null;
  	
  	int m = A.length;
  	int xa = A[0].length;
  	int xb = b.length;
  	int n = 1;
  	double tmpd;
  	
  	if(xa != xb) return null;
  	
  	double[] c = new double[m];
  	
  	for(int j=0; j<m; j++){
  	  for(int i=0; i<n; i++){
  	  	tmpd = 0;
  	  	
  	  	for(int k=0; k<xa; k++){
  	  		tmpd += A[j][k] * b[k];
  	  	}
  	  	
  	  	c[j] = tmpd;
  	  }
  	}
  	
  	return c;
  }

	/**
   * Multiply an array element wise by the given value.
   * @param A the array to multiply
   * @param d the factor
   */
  public static void timesEquals(double[] A, double d)
  {
  	for(int i=0; i<A.length; i++){
  		A[i] *= d;
  	}
  }
  
  /**
   * Multiply a matrix element wise by the given value.
   * @param A the matrix to multiply
   * @param d the factor
   */
  public static void timesEquals(double[][] A, double d)
  {
  	for(int j=0; j<A.length; j++){
  		for(int i=0; i<A[0].length; i++){
  			A[j][i] *= d;
  		}
  	}
  }
  
	/**
   * Perform an element-wise division of the array by the given value.
   *  @param a the array
   *  @param b the divisor
   *  @return the resulting array
   */
  public static double[] divide(double[] a, double b)
  {
    double[] c = new double[a.length];
    
    for(int i=0; i<a.length; i++){
      c[i] = a[i] / b;
    }
    
    return c;
  }
  
	/**
   * Perform an element-wise division of the given arrays.
   *  @param a the array
   *  @param b the array of divisors
   *  @return the resulting array
   */
  public static double[] divide(double[] a, double[] b)
  {
    double[] c = new double[a.length];
    
    for(int i=0; i<a.length; i++){
      c[i] = a[i] / b[i];
    }
    
    return c;
  }
  
	/**
   * Perform an element-wise division of the given matrix.
   *  @param A the matrix
   *  @param b the divisor
   *  @return the resulting array
   */
  public static double[][] divide(double[][] A, double b)
  {
  	int m = A.length;
  	int n = A[0].length;
    double[][] C = new double[m][n];
    
    for(int j=0; j<m; j++){
	    for(int i=0; i<n; i++){
	      C[j][i] = A[j][i] / b;
	    }
    }
    
    return C;
  }
  
  /**
   * Return the dot product of two vectors.
   *  @param a the first vector
   *  @param b the second vector
   *  @return the dot product
   */
  public static double dot(double[] a, double[] b)
  {
    double sum = 0;
    
    for(int i=0; i<a.length; i++){
      sum += a[i]*b[i];
    }
    
    return sum;
  }
  
  /**
   * Calculate the cross product between the given vectors.
   *  @param va
   *  @param vb
   *  @return the cross product
   */
  public static double[] cross(double[] va, double[] vb)
  {
    double[] vc = new double[3];
    
    vc[0] = (va[1] * vb[2]) - (vb[1] * va[2]);
    vc[1] = -(va[0] * vb[2]) + (vb[0] * va[2]);
    vc[2] = (va[0] * vb[1]) - (vb[0] * va[1]);
    
    return vc;
  }

  /**
   * Calculate the trace of a matrix.
   * @param A a square matrix
   * @return the trace of the matrix
   */
  public static double trace(double[][] A)
  {
  	int m = A.length;
  	double sum = 0;
  	
		for(int i=0; i<m; i++){
			sum += A[i][i];
		}
  	
  	return sum;
  }
  
  /**
   * Calculate consecutive element differences.
   * @param a an array
   * @return the approximate derivative
   */
  public static double[] diff(double[] a)
  {
  	double[] b = new double[a.length-1];
  	
  	for(int i=0; i<a.length-1; i++){
  		b[i] = a[i+1] - a[i];
  	}
  	
  	return b;
  }
  
  /**
   * Greater than operator.
   * @param a the array to operate on
   * @param b the value that each element must be greater than
   * @return the resulting array
   */
  public static double[] gt(double[] a, double b)
  {
  	int n = a.length;
  	double[] c = new double[n];
  	
  	for(int i=0; i<n; i++){
  	  if(a[i] > b) c[i] = 1;	
  	}
  	
  	return c;
  }
  
  /**
   * Greater than operator.
   * @param A the matrix to operate on
   * @param b the value that each element must be greater than
   * @return the resulting matrix
   */
  public static double[][] gt(double[][] A, double b)
  {
  	int m = A.length;
  	int n = A[0].length;
  	double[][] C = new double[m][n];
  	
  	for(int j=0; j<m; j++){
	  	for(int i=0; i<n; i++){
	  	  if(A[j][i] > b) C[j][i] = 1;	
	  	}
  	}
  	
  	return C;
  }
  
  /**
   * Less than operator.
   * @param a the array to operate on
   * @param b the value that each element must be less than
   * @return the resulting array
   */
  public static double[] lt(double[] a, double b)
  {
  	int n = a.length;
  	double[] c = new double[n];
  	
  	for(int i=0; i<n; i++){
  	  if(a[i] < b) c[i] = 1;	
  	}
  	
  	return c;
  }
  
  /**
   * Less than operator.
   * @param A the matrix to operate on
   * @param b the value that each element must be less than
   * @return the resulting matrix
   */
  public static double[][] lt(double[][] A, double b)
  {
  	int m = A.length;
  	int n = A[0].length;
  	double[][] C = new double[m][n];
  	
  	for(int j=0; j<m; j++){
	  	for(int i=0; i<n; i++){
	  	  if(A[j][i] < b) C[j][i] = 1;	
	  	}
  	}
  	
  	return C;
  }
  
	/**
   * Save a matrix to an ASCII file.
   * @param filename the filename to save to
   * @param A the matrix to save
   */
  public static void save(String filename, double[][] A)
  {
  	int m = A.length;
  	int n = A[0].length;
    
    try{    
      BufferedWriter outs = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(filename)));
      
      for(int j=0; j<m; j++){
        for(int i=0; i<n; i++){
          outs.write(A[j][i] + " ");
        }
        
        outs.newLine();
      }
      
      outs.close();
    }catch(Exception e){
      e.printStackTrace();
    }
  }
  
  /**
   * Save an array to an ASCII file.
   * @param filename the filename to save to
   * @param a the array to save
   */
  public static void save(String filename, double[] a)
  {    
    try{    
      BufferedWriter outs = new BufferedWriter(new FileWriter(filename));
      
      for(int i=0; i<a.length; i++){
        outs.write(a[i] + " ");
      }
      
      outs.close();
    }catch(Exception e){e.printStackTrace();}
  }
  
  /**
   * Normalize a point.
   * @param x a point
   */
  public static void normalize(double[] x)
  {
  	double tmpd = norm(x);
  	
  	for(int i=0; i<x.length; i++){
  		x[i] /= tmpd;
  	}
  }
  
  /**
   * Normalize a set of points.
   * @param P a set of points
   * @param mean the mean of the group of points
   * @param std the standard deviation of the points
   */
	public static void normalize(double[][] P, double[] mean, double[] std)
	{
		for(int i=0; i<P.length; i++){
			for(int j=0; j<P[i].length; j++){
				P[i][j] = (P[i][j] - mean[j]) / std[j];
			}
		}
	}
	
	/**
	 * Calculate the centroid of a group of points.
	 * @param X a group of points
	 * @return the centroid
	 */
	public static double[] mean(double[][] X)
	{
		double[] M = null;
		
    if(X != null && X.length > 0){
      int d = X[0].length;
      
      M = new double[d];
      
      for(int j=0; j<d; j++){
      	for(int i=0; i<X.length; i++){
          M[j] += X[i][j];
        }
      }
      
      for(int j=0; j<d; j++){
        M[j] /= (double)X.length;
      }
    }
    
    return M;
	}
	
	/**
	 * Calculate the standard deviation of each coordinate in a group of points.
	 * @param X a group of points
	 * @param mean the centroid of the points
	 * @return the standard deviation of each coordinate
	 */
  public static double[] std(double[][] X, double[] mean)
  {
		double[] S = null;
		
    if(X != null && X.length > 0){
      int d = X[0].length;
      double tmpd;
      
      S = new double[d];
              
      for(int j=0; j<d; j++){
      	for(int i=0; i<X.length; i++){
      		tmpd = X[i][j] - mean[j];
          S[j] += tmpd * tmpd;
        }
      }
      
      for(int j=0; j<d; j++){
        S[j] /= (double)X.length;
        S[j] = Math.sqrt(S[j]);
      }
    }
    
    return S;
  }
  
  /**
   * Calculate a specified moment of a group of 2D points.
   * @param X a group of 2D points
   * @param p the moment of the x-coordinate
   * @param q the moment of the y-coordinate
   * @return the moment
   */
	public static double moment(double[][] X, double p, double q)
	{
		double m = 0;
		
		for(int i=0; i<X.length; i++){
			m += Math.pow(X[i][0],p) + Math.pow(X[i][1],q);
    }
		
		m /= X.length;
		
		return m;
	}
	
	/**
   * A simple main for debug purposes only.
   *  @param args arguments to the program
   */
  public static void main(String args[])
  {
  	if(false){
	  	double[][] A = new double[][]{{ 1,  2,  3,  4},
	  		                            { 5,  6,  7,  8},
	  		                            { 9, 10, 11, 12},
	  		                            {13, 14, 15, 16}};
	  	
	  	double[][] B = new double[][]{{ 1,  2,  3,  4},
	        													{ 5,  6,  7,  8}};
	  	
	  	println(A);
	  	println(transpose(A));
	  	println(mtimes(A, A));
	  	println(mtimes(A, transpose(B)));
  	}
  	
  	if(true){
  		double[] f = new double[]{1, 1, 10, 10, 1, 1, 10, 10, 1, 1, 10, 10, 1, 1};
  		double[] c = dct(f);
  		double[] g = idct(c);
  		
  		println(f);
  		println(c);
  		println(g);
  	}
  }
}