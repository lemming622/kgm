package kgm.image;
import kgm.utility.*;
import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.*;

/**
 * Select regions within an image.
 * @author Kenton McHenry
 */
public class SelectRegions extends JFrame implements MouseListener, MouseMotionListener, MouseWheelListener, KeyListener
{
	private Vector<String> images = new Vector<String>();
	private String image_path = null;
	private String segmentation_path = null;
	private String output_path = null;
	private int index = 0;
	
	private int[][] Irgb;
	private int[][] Is = null;
	private int width, height;
  private TreeSet<Integer> selected_regions = new TreeSet<Integer>();
  private boolean DIRTY = false;
  
	private int mouse_button = -1;
  private int offx = 3;
  private int offy = 26;
  
	/**
	 * Class constructor.
	 * @param image_path the path or name of the image file(s)
	 * @param segmentation_path the path to the segmented regions
	 * @param output_path the path to store selected region output
	 */
	public SelectRegions(String image_path, String segmentation_path, String output_path)
	{		
		this.image_path = image_path;
		this.segmentation_path = segmentation_path;
		this.output_path = output_path;
				
		if(!image_path.endsWith("/")){
			images.add(Utility.getFilename(image_path));
			this.image_path = Utility.getFilenamePath(image_path);
		}else{
			File[] files = new File(image_path).listFiles();
			
			for(int i=0; i<files.length; i++){
				images.add(Utility.getFilename(files[i].getName()));
			}
		}
		
    addMouseListener(this);
    addMouseMotionListener(this);
    addMouseWheelListener(this);
    addKeyListener(this);
    
    updateImage(0);
	}
	
	public void updateImage(int increment)
	{
		//Save selected regions
		if(DIRTY){			
		  try{
		    BufferedWriter outs = new BufferedWriter(new FileWriter(output_path + Utility.getFilenameName(images.get(index)) + ".txt", true));
		    
		    for(Iterator<Integer> itr = selected_regions.iterator(); itr.hasNext();){
		    	outs.write(itr.next() + "\n");
		    }

		    outs.close();
		  }catch(Exception ex) {ex.printStackTrace();}
		  
		  DIRTY = false;
		}
		
		//Move the image index
		if(index + increment < 0){
			index = images.size() - 1;
		}else if(index + increment >= images.size()){
			index = 0;
		}else{		
			index += increment;
		}
		
		//Load new image
		Irgb = ImageUtility.load(image_path + images.get(index));
  	height = Irgb.length;
  	width = Irgb[0].length;		
  	    
  	//Load segmentations
		try{
			FileInputStream fis = new FileInputStream(segmentation_path + Utility.getFilenameName(images.get(index)) + ".ser");
			ObjectInputStream ois = new ObjectInputStream(fis);
			Is = (int[][])ois.readObject();
			ois.close();
			fis.close();
		}catch(Exception e) {e.printStackTrace();}
		
		//Re-scale image to match segmentation if need be
		if(Is.length != height || Is[0].length != width){
			height = Is.length;
			width = Is[0].length;
			Irgb = ImageUtility.resize(Irgb, width, height);
		}
  	
		//Load selected regions
		selected_regions.clear();
    
		try{
			Scanner scanner = new Scanner(new File(output_path + Utility.getFilenameName(images.get(index)) + ".txt"));
			
			while(scanner.hasNextInt()){
				selected_regions.add(scanner.nextInt());
			}
		}catch(Exception e) {}    
    
  	setTitle((index+1) + ": " + images.get(index));
  	setSize(width, height+30);
    repaint();
	}
	
	/**
	 * Select the region at the given coordinate.
	 * @param x the x-coordinate
	 * @param y the y-coordinate
	 */
	public void select(int x, int y)
	{
		if(x >= 0 && x < width && y >= 0 && y < height){
			selected_regions.add(Is[y][x]);
			DIRTY = true;
		}
		
		repaint();
	}
	
	/**
	 * Un-select the region at the given coordinate.
	 * @param x the x-coordinate
	 * @param y the y-coordinate
	 */
	public void unselect(int x, int y)
	{
		if(x >= 0 && x < width && y >= 0 && y < height){
			selected_regions.remove(Is[y][x]);
			DIRTY = true;
		}
		
		repaint();
	}
	
	/**
	 * Paint the frame.
	 * @param g the graphics context to paint to
	 */
  public void paint(Graphics g)
  {
  	int[][] Itmp = ImageUtility.n2argb(Is, Irgb);
  	
  	//Color selected regions
  	for(int x=0; x<width; x++){
  		for(int y=0; y<height; y++){
  			if(selected_regions.contains(Is[y][x])){
  				Itmp[y][x] = 0x00ff0000;
  			}
  		}
  	}
  	
		g.drawImage(ImageUtility.argb2image(Itmp), offx, offy, null);
  }
  
  /**
   * Listener for mouse pressed events.
   * @param e the mouse event
   */
  public void mousePressed(MouseEvent e)
  {     
  	mouse_button = e.getButton();
  	
    if(mouse_button == MouseEvent.BUTTON1){
    	select(e.getX()-offx, e.getY()-offy);
    }else if(mouse_button == MouseEvent.BUTTON3){
    	unselect(e.getX()-offx, e.getY()-offy);
    }
  }
  
  /**
   * Listener for mouse drag events.
   * @param e the mouse event
   */
	public void mouseDragged(MouseEvent e) 
	{
    if(mouse_button == MouseEvent.BUTTON1){
    	select(e.getX()-offx, e.getY()-offy);
    }else if(mouse_button == MouseEvent.BUTTON3){
    	unselect(e.getX()-offx, e.getY()-offy);
    }
	}
	
  /**
   * Mouse wheel listener that handles changing the currently displayed image.
   * @param e the mouse wheel event
   */
  public void mouseWheelMoved(MouseWheelEvent e)
  {  	
  	updateImage(e.getWheelRotation() < 0 ? -1 : 1);    
  }
	
  /**
   * Key listener.
   * @param e a key event
   */
	public void keyPressed(KeyEvent e)
	{
		if(e.getKeyChar() == 'q'){
			updateImage(0);
			dispose();
		}
	}
	
	public void mouseMoved(MouseEvent e) {}
	public void mouseClicked(MouseEvent e) {}
	public void mouseReleased(MouseEvent e) {}
	public void mouseEntered(MouseEvent e) {}
	public void mouseExited(MouseEvent e) {}	
	public void keyTyped(KeyEvent e) {}
	public void keyReleased(KeyEvent e) {}	
	
	/**
	 * Main for a program to select regions within an image.
	 * @param args command line arguments
	 */
	public static void main(String args[])
	{
		String image_path = null;
		String segmentation_path = null;
		String output_path = null;
		SelectRegions lr;
		
		if(args.length == 0){
			//args = new String[]{"-s", "tmp/", "-o", "tmp/", "C:/Users/kmchenry/Files/Data/Images/scar1.jpg"};
			args = new String[]{"-s", "//isda.ncsa.uiuc.edu/BP/iReq215/NCSA_Keywords/Ocean-Segmented/", "-o", "tmp/", "//isda.ncsa.uiuc.edu/BP/iReq215/NCSA_Keywords/Ocean/"};
		}
		
		//Process command line arguments
		for(int i=0; i<args.length; i++){
			if(args[i].equals("-?")){
				System.out.println("Usage: SelectRegions [options] [image/images path]");
				System.out.println();
				System.out.println("Options: ");
				System.out.println("  -?: display this help");
				System.out.println("  -s path: set the segmenation path");
				System.out.println("  -o path: set the output path");
				System.out.println();
				System.exit(0);
			}else if(args[i].equals("-s")){
				segmentation_path = args[++i];
			}else if(args[i].equals("-o")){
				output_path = args[++i];
			}else{
				image_path = args[i];
			}
		}
		
		if(segmentation_path == null) segmentation_path = Utility.getFilenamePath(image_path);
		if(output_path == null) output_path = Utility.getFilenamePath(image_path);

		lr = new SelectRegions(image_path, segmentation_path, output_path);
		lr.setVisible(true);
	}
}