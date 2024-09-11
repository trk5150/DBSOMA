package mapScanning;

import java.awt.Color;
import java.awt.Font;
import java.awt.GradientPaint;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Polygon;
import java.awt.AWTException;
import java.awt.Color;
import java.awt.Font;
import java.awt.GradientPaint;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.Robot;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;
import javax.imageio.ImageIO;
import javax.swing.JPanel;

public class DrawSet extends JPanel
{
	public int minDataPoints , maxDataPoints;
	public Color colorA=Color.red;
	public Color colorB=Color.blue;
	public int height = 200;
	public int width = 200;
	public DrawSet d;
	public int xNodes, yNodes, minPts, fontSize;
	public ArrayList<Polygon> nodeHexs;
	public ArrayList<Color> nodeColors;
	public double[][] polyCenterCoords;
	public int[] nodeWeights;
	public boolean scanner;
	public DrawSet(int min, int xN, int yN, int[] w)
	{
		scanner = false;
		minPts = min;
		xNodes =xN;
		yNodes = yN;
		nodeWeights = w;
		nodeColors = new  ArrayList<Color>();
		intitializeColors(xNodes*yNodes);

		polyCenterCoords = new double[xNodes*yNodes][2];
		nodeHexs = new ArrayList<Polygon>();
	}
	
	public void paintComponent(Graphics g)
	{
	    super.paintComponent(g);
	    colorBar(g);
	    //setBackground(Color.BLUE);

	    for(int i = 0; i < nodeHexs.size(); i++)
	    {
		    g.setColor(Color.WHITE);
		    g.fillPolygon(nodeHexs.get(i));
	  		g.setColor(nodeColors.get(i));
	  		g.fillPolygon(nodeHexs.get(i));
	    }
	}
	public void projectSet(int[][] values, boolean scanned)
	{
		int[] projection = new int[values.length];
		int index = 0;
		scanner = scanned;
		if(scanned)
		{
//			colorA = Color.blue;
			index = 1;
		}
		int min = 0;
		int max = 0;
		for(int i = 0; i < values.length; i++)
		{
			int value = values[i][index];
			if(scanned)
			{
				value = 0;
				if(values[i][1]>minPts)
					value = nodeWeights[i];
			}
			
			projection[i] = value;
			if(value> max)
				max = value;
			if(value < min)
				min = value;
		}
		heatMapping(min, max, projection);
		repaint();
		
	}
	public void nodeBuild(int winW, int winH)
	{  
		int xMin = (int)(.1 * winW);
	    int yMin = (int)(.1 * winH);
	    int xMax = (int)(.9 * winW);
	    int yMax = (int)(.9 * winH);
	    fontSize = winH/85;
	    
		
	    double winWidth = xMax-xMin;
	    double winHeight = yMax-yMin;
	    
	    double height = winHeight/yNodes;
	    double width = winWidth/xNodes;
	    height = height/3;
	    width = width/2;
	    int x = (int) (width);
	    int y = (int) (height);
	    int w = (int) (winWidth/width);
	    int h = (int)(winHeight/height);
	    int centerCount = 0;
	    for(int i = 0; i<h; i++)
		{
	    	for(int m = 0; m<w; m++)
			{
				if((i%6==1 && m%2==0)||(i%6==4 && m%2 == 1))
				{
					double xLoc = m*x+xMin;
					double yLoc = i*y+yMin;
					polyCenterCoords[centerCount][0] = xLoc;
					polyCenterCoords[centerCount][1] = yLoc;
					centerCount++;
				}
			}
		}
	    for(int c = 0; c<polyCenterCoords.length; c++)
	    {
		   
	    	int[] xPoints = new int[6];
			int[] yPoints = new int[6];
			  
			xPoints[0] = (int) polyCenterCoords[c][0];
			yPoints[0] = (int) polyCenterCoords[c][1]-(2*y);
			  
			xPoints[1] = (int) polyCenterCoords[c][0]+x;
			yPoints[1] = (int) polyCenterCoords[c][1]-y;
			  
			xPoints[2] = (int) polyCenterCoords[c][0]+x;
			yPoints[2] = (int) polyCenterCoords[c][1]+y;
				
			xPoints[3] = (int) polyCenterCoords[c][0];
			yPoints[3] = (int) polyCenterCoords[c][1]+(2*y);
	
		  	xPoints[4] = (int) polyCenterCoords[c][0]-x;
		  	yPoints[4] = (int) polyCenterCoords[c][1]+y;
				
		  	xPoints[5] = (int) polyCenterCoords[c][0]-x;
		  	yPoints[5] = (int) polyCenterCoords[c][1]-y;
		  
		  	nodeHexs.add(polygonMaker(xPoints, yPoints, 6));
	    }
	}

	public Polygon polygonMaker(int[] x, int[] y, int i)
	{
		Polygon p = new Polygon(x,y,i);
		return p;
	}
	public void colorBar(Graphics g)
	{
		Graphics2D g2d = (Graphics2D)g;
		int x = (int)(getWidth()*.2);
		int y = (int)(getHeight()*.96);
		int width = (int)(getWidth()*.55);
		int height = (int)(getHeight()*.02);
		
		GradientPaint colorbar = new GradientPaint(x, y, Color.white, x+width, y, colorA, false);
		if(scanner)
			colorbar = new GradientPaint(x, y, Color.white, x+width, y, colorB, false);
		g2d.setPaint(colorbar);
		g2d.fillRect(x, y, width, height);
		
		g2d.setPaint(Color.black);
		g2d.setColor(Color.BLACK);
		g2d.drawRect(x, y, width, height);
		g2d.drawString(""+minDataPoints, (int)(getWidth()*.17),(int)(getHeight()*.95));
		g2d.drawString(""+maxDataPoints, (int)(getWidth()*.2+(int)(getWidth()*.55)),(int)(getHeight()*.95));
		
	}
	public void intitializeColors(int num)
	{
		for(int i = 0; i < num; i++)
		{
			Color color = new Color(colorA.getRed(), colorA.getGreen(),colorA.getBlue(), 255);
			nodeColors.add(color);
		}
	}
	public void heatMapping(int min, int max, int[] projections)
	{
		minDataPoints = min;
		maxDataPoints = max;
		for(int i = 0; i<nodeColors.size(); i++)
		{	
			double d = projections[i];
			double doop = max;
			Color color = new Color(colorA.getRed(), colorA.getGreen(),colorA.getBlue(),(int)(255*d/doop));
			if(scanner)
			{
				color = new Color(colorB.getRed(), colorB.getGreen(),colorB.getBlue(),(int)(255*d/doop));
			}
			nodeColors.set(i, color);
		}
	}
	public void saveImg(String name, File outDir)
	{
		BufferedImage image = null;
		try {
			image = new Robot().createScreenCapture(new Rectangle(this.getLocationOnScreen().x+62, this.getLocationOnScreen().y+30, 511, 457));
		} catch (AWTException e1) {
			e1.printStackTrace();
		}
	    try 
	    {
	        File outputfile = new File(outDir.getAbsolutePath()+File.separator+name+".png");
	        ImageIO.write(image, "png", outputfile);
	    } catch (IOException e) {}
	}
	public void saveImg(String name, File outDir, int w, int h)
	{
		this.setVisible(true);
		this.setSize(w, h);
	    BufferedImage bi = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);
	    Graphics2D g = bi.createGraphics();
	    print(g);
	    try 
	    {
	        File outputfile = new File(outDir.getAbsolutePath()+File.separator+name+".png");
	        ImageIO.write(bi, "png", outputfile);
	    } catch (IOException e) {}
	}
}
