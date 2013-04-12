package kgm.image;
import static com.googlecode.javacv.cpp.opencv_highgui.*;
import static com.googlecode.javacv.cpp.opencv_imgproc.*;
import java.io.*;
import java.util.*;
import com.googlecode.javacv.*;
import com.googlecode.javacv.cpp.opencv_core.*;
import mpicbg.imagefeatures.*;
import kgm.image.ImageUtility.*;
import kgm.matrix.*;
import kgm.utility.Pair;

/**
 * A demonstration tool
 * @author Kenton McHenry
 */
public class Demo{
	/**
	 * The main.
	 * @param args arguments to the program
	 */
	public static void main(String args[])
	{
		String filename = "C:/Users/kmchenry/Files/Data/Images/scar1.jpg";
		
		//Convolution
		if(false){
	  	ImageViewer viewer = new ImageViewer();
	  	int[][] tmp = ImageUtility.load(filename);
	  	int h = tmp.length;
	  	int w = tmp[0].length;
	  	int[] image_rgb = ImageUtility.to1D(tmp);
	  	double[] image_gray = ImageUtility.argb2g(image_rgb, w, h);
	  	double[][] filter = ImageUtility.getFilter(Option.EDGE, 2, 1, 45);	filter = MatrixUtility.transpose(filter);
	  	
	  	viewer.add(image_rgb, w, h, true);
	  	viewer.add(image_gray, w, h, true);
	  	viewer.add(MatrixUtility.to1D(filter), filter[0].length, filter.length, true);
	  	viewer.add(ImageUtility.convolve(image_gray, w, h, filter), w, h, true);
		}
		
		//Canny
		if(false){
			int[][] image_rgb = ImageUtility.load(filename);
	  	int h = image_rgb.length;
	  	int w = image_rgb[0].length;
			double[] image_gray = JavaCVUtility.getEdges(image_rgb, 800, 800, 5);

			ImageViewer.show(image_gray, w, h, 1000, "Edges");
		}
		
		//Hough
		if(false){
			int[][] tmp = ImageUtility.load(filename);
			int[] image_rgb = ImageUtility.to1D(tmp);
			int h = tmp.length;
			int w = tmp[0].length;
			
			double[] image_edges = JavaCVUtility.getEdges(tmp, 200, 200, 3);
			Vector<Pair<Double,Double>> lines = JavaCVUtility.getLines(image_edges, w, h, 1, 40);
			Vector<Pixel[]> line_segments = JavaCVUtility.getLineSegments(lines, w, h);
			
			System.out.println("Lines: " + lines.size());
			ImageUtility.drawLines(image_rgb, w, h, line_segments, 0x000000ff);		
			ImageViewer viewer = ImageViewer.show(image_edges, w, h, "Lines");
			viewer.add(image_rgb, w, h, true);
		}
		
		//Harris
		if(false){
			int[][] tmp = ImageUtility.load(filename);
			int h = tmp.length;
			int w = tmp[0].length;
			int[] image_rgb = ImageUtility.to1D(tmp);
			double[] image_gray = ImageUtility.argb2g(tmp);
			
			Vector<Pixel> corners = JavaCVUtility.getCorners(image_gray,  w, h, 5, 3, 0.04, 0.001);
			
			System.out.println("Corners: " + corners.size());
			ImageUtility.drawBoxes(image_rgb, w, h, corners, 4, 0x00ffff00);
			ImageViewer.show(image_rgb, w, h, 1000, "Corners");
		}
		
		//Oversegmentation
		if(true){
	  	ImageViewer viewer = new ImageViewer(1000);
	  	int[][] image_rgb = ImageUtility.load(filename);
	  	int h = image_rgb.length;
	  	int w = image_rgb[0].length;
	  	int[][] image_segments;
	  	
	  	viewer.add(image_rgb, w, h, false);
	  	image_segments = ImageUtility.getSuperPixels(image_rgb, 0.2, true);
	  	viewer.add(ImageUtility.n2argb(image_segments), w, h, true);
		}
		
		//Active Contours
		if(false){
	  	int[][] tmp = ImageUtility.load(filename);
	  	int h = tmp.length;
	  	int w = tmp[0].length;
	  	double[][] image_gray = MatrixUtility.to2D(h, w, ImageUtility.argb2g(tmp));
	  	
  		int h_small = (int)Math.round(0.5*h);
  		int w_small = (int)Math.round(0.5*w);
  		image_gray = ImageUtility.resizeBicubic(image_gray, w_small, h_small);
  		image_gray = ImageUtility.smooth(image_gray, 4);
  		
  		ImageViewer.show(image_gray, w_small, h_small);
  		ImageUtility.getGeodesicActiveContour(image_gray, 0.1, 100, 0.5, 3000);
		}
		
		//SIFT
		if(false){
			int[][] image_rgb = ImageUtility.load(filename);
			Vector<Feature> features = ImageJUtility.getSIFTFeatures(image_rgb);
			
			ImageJUtility.drawFeatures(image_rgb, features, 0x00ffff00);
			ImageViewer.show(image_rgb, -1, -1, 1000, "");
		}
		
		//Thinning
		if(false){
			int[][] image_rgb = ImageUtility.load(filename);
			int h = image_rgb.length;
			int w = image_rgb[0].length;
			double[] image_gray = ImageUtility.argb2g(image_rgb);
			double[] image_bw = ImageUtility.g2bw(image_gray, w, h, 0.5);
			
			image_bw = MatrixUtility.abs(MatrixUtility.plus(image_bw, -1));
			ImageUtility.hilditch(image_bw, w, h);
			ImageViewer.show(image_bw, w, h, "Thinning");
		}
		
		//SURF
		if(false){
			int[][] image_rgb = ImageUtility.load(filename);
			double[] image_gray = ImageUtility.argb2g(image_rgb);
			int h = image_rgb.length;
			int w = image_rgb[0].length;
			
			Pair<Vector<Ellipse>,Vector<double[]>> features = JavaCVUtility.getSURFFeatures(image_gray, w, h);
			
			System.out.println("Features: " + features.first.size());
			ImageUtility.drawEllipses(image_rgb, features.first, 0x00ffff00);
			ImageViewer.show(image_rgb, w, h, 1000, "SURF");
		}
	}
}