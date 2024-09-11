package mapDrawing;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Scanner;

import som_analysis.SingleReader;

public class PaintSelectFiles 
{
	public int radius;
	public int minPts;
	public DrawHex d;
	public String fold, mapper;
	public String[] geneStates;
	public File outDir, trimmedOutDir, imageOutDir;
	public double gini, pVal;
	public int nodes,xs,ys,xNodes,yNodes, maxBins, minBins, colorNum, winW, winH, pValnVal, binSize;
	public ArrayList<DataPoint> bins;
	public ArrayList<MiniNode> nodeList;
	public ArrayList<File> searchers;
	public MiniSystem nodeSystem;
	public boolean addToOldFile, pics, equalWeight, trimmedPic;
	public boolean printGiniCurves=false;
	public double[][] g;
	public double[] sep;
	public double giniWLArea =0;
	public Color heatmapColorA=Color.red, heatmapColorB=Color.blue;
	
	public PaintSelectFiles(String map, ArrayList<File> fillle, File outie, String[] states, int rad, int min)
	{
		radius = rad;
		minPts = min;
		geneStates = states;
		trimmedPic = true;
		mapper = map;
		outDir = outie;
		trimmedOutDir = new File(outDir + "\\Trimmed");
		imageOutDir = new File(outDir + "\\Images");
		String mapp = map;
		yNodes=0;
		xNodes=0;
		mapReader(mapp);
		binSize = bins.get(0).maxLocus - bins.get(0).minLocus;
		searchers = new ArrayList<File>();
		for (int i = 0; i < fillle.size(); i++) 
		{
			File file = fillle.get(i);
		    if (file.isFile() && file.getName().endsWith(".grp") || file.getName().endsWith(".txt")) 
		    {
		    	searchers.add(file);
		    }
		}
		nodes = yNodes*xNodes;
		
		heatMapping();
	}
	public void setColorA(Color c) {heatmapColorA=c;}
	public void setColorB(Color c) {heatmapColorB=c;}
	
	
	public void searchSystem()
	{
		outDir.mkdir();
		imageOutDir.mkdir();
		d = new DrawHex(mapper);
		d.setColorA(heatmapColorA);
		d.setColorB(heatmapColorB);
		d.setEqualWeight(equalWeight);
		d.nodeBuild(2000,2000);
    	//d.colors();
    	d.heatMapping();
    	d.countingDPS(null);
    	d.saveImg("whole_map", outDir,  2000, 2000);
		
		SingleReader r = new SingleReader(false);
		r.readSOM(mapper);
		r.mapMaker();
		
		
		ArrayList<ArrayList<String>> stateLists = new ArrayList<ArrayList<String>>();
		for(String fifi : geneStates)
		{
			File fifo = new File(fifi);
			stateSearch(fifo, r);			
		}
		for (File s:searchers)
		{
			search(s, r, stateLists);
		}
		System.out.println("done");
	}
	
	public void search(File file, SingleReader r, ArrayList<ArrayList<String>> stateLists)
	{
		int shell = radius;
		int requiredCalls = minPts;
		
		for(int i = 0; i < nodeSystem.size(); i++)
		{
			MiniNode mini = nodeSystem.get(i);
			mini.counting.clear();
			mini.weight = 0;
		}
		
		BufferedReader br = null;
		ArrayList<Integer> inds = new ArrayList<Integer>();
		String line = "";
		try 
		{
			br = new BufferedReader(new FileReader(file));
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		ArrayList<String> foundGenes = new ArrayList<String>();
	    try {
            line = br.readLine();

            while (line != null) 
	        {
	        	String[] s = new String[1];
	        	s[0] = line;
	        	if(line.contains(","))
	        	{
	        		s = line.split(",");
	        	}
//	        	if(s.length>0)
//	        		System.out.println(s[0] + " " + s.length);
	        	
	        	for(int k = 0; k<s.length; k++)
	        	{
	        		
	        		if(s[k].length()>1)
	        		{
		        		r.geneCheck(s[k], shell);//for trim analysis
			        	for(int i = 0; i< d.dataPoints.size(); i++)
			        	{
			        		if(d.dataPoints.get(i).name.equalsIgnoreCase(s[k].trim()))
			        		{
			        			d.dataPoints.get(i).myMini.counting.add(d.dataPoints.get(i));
			        			inds.add(i);
			        			foundGenes.add(s[k]);
			        		}
			        		
			        	}
		        	}
	        	}
	        	line = br.readLine();
        	}
//            System.out.println("searched for " + tryCount + " genes, found "+foundCount);   
		    } catch (IOException e)
		    {
				e.printStackTrace();
			}
	 
	  		ArrayList<String> trimmedList = r.countMultipleRepeatedGenes(requiredCalls, false); //for trim analysis
			
			{
				d.search(file);
				d.repaint();
				d.saveImg(file.getName(), imageOutDir, 2000, 2000);
				if(trimmedPic)
				{
					trimmedOutDir.mkdir();
					for(int i = 0; i < nodeSystem.size(); i++)
					{
						MiniNode mini = nodeSystem.get(i);
						mini.counting.clear();
						mini.weight = 0;
					}
					
					//flip colors for the trimmed maps for easier viewing
					d.setColorA(heatmapColorB);
					d.setColorB(heatmapColorA);
					
					//search the Trimmed List of genes
					d.search(trimmedList);
					d.repaint();
					d.saveImg(file.getName() + "_trimmed", trimmedOutDir, 2000, 2000);
					
					//reset colors
					d.setColorA(heatmapColorA);
					d.setColorB(heatmapColorB);
				}
			}
			
			
			r.genio.clear();
	}
	
	public void stateSearch(File file, SingleReader r)
	{
		int shell = radius;
		int requiredCalls = minPts;
		
		for(int i = 0; i < nodeSystem.size(); i++)
		{
			MiniNode mini = nodeSystem.get(i);
			mini.counting.clear();
			mini.weight = 0;
		}
		
		BufferedReader br = null;
		ArrayList<Integer> inds = new ArrayList<Integer>();
		String line = "";
		try 
		{
			br = new BufferedReader(new FileReader(file));
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	    try {
            line = br.readLine();

            while (line != null) 
	        {
	        	String[] s = new String[1];
	        	s[0] = line;
	        	if(line.contains(","))
	        	{
	        		s = line.split(",");
	        	}
//	        	if(s.length>0)
//	        		System.out.println(s[0] + " " + s.length);
	        	for(int k = 0; k<s.length; k++)
	        	{
	        		if(s[k].length()>1)
	        		{
		        		r.geneCheck(s[k], shell);//for trim analysis
			        	for(int i = 0; i< d.dataPoints.size(); i++)
			        	{
			        		if(d.dataPoints.get(i).name.equalsIgnoreCase(s[k].trim()))
			        		{
			        			d.dataPoints.get(i).myMini.counting.add(d.dataPoints.get(i));
			        			inds.add(i);
			        		}
			        		
			        	}
		        	}
	        	}
	        	line = br.readLine();
        	}
//            System.out.println("searched for " + tryCount + " genes, found "+foundCount);   
		    } catch (IOException e)
		    {
				e.printStackTrace();
			}
	    	
	 
	  		ArrayList<String> trimmedList = r.countMultipleRepeatedGenes(requiredCalls, false); //for trim analysis
			
		
			if(trimmedPic)
			{
				outDir.mkdir();
				for(int i = 0; i < nodeSystem.size(); i++)
				{
					MiniNode mini = nodeSystem.get(i);
					mini.counting.clear();
					mini.weight = 0;
				}
				
				//flip colors for the trimmed maps for easier viewing
				d.setColorA(heatmapColorB);
				d.setColorB(heatmapColorA);
				
				//search the Trimmed List of genes
				d.search(trimmedList);
				d.repaint();
				d.saveImg(file.getName() + "_STATE_", outDir, 2000, 2000);
				
				//reset colors
				d.setColorA(heatmapColorA);
				d.setColorB(heatmapColorB);
			}
			
			r.genio.clear();
	}
	
	public void degreesOfSep()
	{
		double[][] dist = new double[nodes][nodes];
		for (int i = 0; i<nodes; i++)
		{
			int xCoordA = i%xNodes;
			int yCoordA = i/yNodes;
			for (int j = 0; j<nodes; j++)
			{
				int xCoordB = j%xNodes;
				int yCoordB = j/yNodes;
				
				double ddd = hexDist(xCoordA, yCoordA, xCoordB, yCoordB);
				dist[i][j] = ddd;
			}
		}
		
		for(int i = 0; i<g[0].length; i++)
		{
			double count = 0;
			for(int j = 0; j<nodes; j++)
			{
				for(int k = 0; k<nodes; k++)
				{
					sep[i] += dist[j][k]*g[j][i]*g[k][i];
					count+= g[j][i]*g[k][i];
				}
			}
			sep[i] /= count;
		}
	}
	public double degreesP()
	{
		double slit = 0;
		double avg = 0;
		for(int i = 1; i < sep.length; i++)
		{
			avg += sep[i];
		}
		avg /= pValnVal;
		double std = 0;
		for(int i = 1; i < pValnVal+1; i++)
		{
			std += (sep[i]-avg)*(sep[i]-avg);
		}
		if(std==0)
			std = 0.00001;
		std /= pValnVal;
		std = Math.sqrt(std);
//		NormalDistribution norm = new NormalDistribution(avg,std);
//		slit = norm.cumulativeProbability(sep[0]);
		return slit;
	}
	public double hexDist(int x1, int y1, int x2, int y2)
	{
		if(x1==x2 && y2 == y1)
			return 0;
		else
		{
			int hm = 0;
			int vm = 0;
			int cm = 0;
			boolean right = x1-x2<0;
			hm = Math.abs(x1-x2);
			vm = Math.abs(y1-y2);
			if(xNodes-x1<xNodes-x2)
				hm = Math.min(xNodes-x1 + x2,hm);
			else
				hm = Math.min(xNodes-x2 +x1, hm);
			
			if(yNodes-y1<yNodes-y2)
				vm = Math.min(yNodes-y1 + y2,vm);
			else
				vm = Math.min(yNodes-y2 +y1, vm);
			
			if(vm %2 !=0)
			{
				if(!right && y1%2 == 0)
					cm = Math.min(hm, (int)(((double)vm)/2 +.5));
				else if(right && y1%2 == 1)
					cm = Math.min(hm, (int)(((double)vm)/2 +.5));
			}
			else
				cm = Math.min(hm, (int)(((double)vm)/2));
			return hm+vm-cm; 
		}
	}
	
	
	public void mapReader(String f)
	{
		ArrayList<String> StringMat = new ArrayList<String>();
		try 
		{
			@SuppressWarnings("resource")
			Scanner in = new Scanner(new FileReader(f));
			String sizer = in.next();
			int xo = Integer.parseInt(sizer.substring(0,sizer.indexOf("x")));
			int yo = Integer.parseInt(sizer.substring(sizer.indexOf("x")+1,sizer.length()));
			in.next();
			while(in.hasNext())
			{
				StringMat.add(in.next());
			}
			createNoder(StringMat, xo, yo);
		} catch (FileNotFoundException e) {e.printStackTrace();}
	}
	
	public void createNoder(ArrayList<String> strings, int x, int y)
	{
		nodeList = new ArrayList<MiniNode>(x*y);
		bins = new ArrayList<DataPoint>();
		MiniNode current = null;
		xNodes = x;
		yNodes = y;
		for(String whole: strings)
		{
			if(whole.contains("["))
			{
				current = new MiniNode(whole);
				nodeList.add(current);
			}
			else
				current.addDP(whole);
		}
		for(int i = 0; i<nodeList.size(); i++)
		{
			for(int j = 0; j < nodeList.get(i).bins.size(); j++)
			{
				bins.add(nodeList.get(i).bins.get(j));
			}
		}
		nodeSystem = new MiniSystem(nodeList,xNodes,yNodes);
	}
	
	public void heatMapping()
	{
		minBins = (int) (nodeList.get(0).bins.size());
		maxBins = (int) (nodeList.get(0).bins.size());
		for(int i = 0; i<nodeList.size(); i++)
		{
			if(nodeList.get(i).bins.size()>maxBins)
				maxBins = (int) (nodeList.get(i).bins.size());
			if(nodeList.get(i).bins.size()<minBins)
				minBins = (int) (nodeList.get(i).bins.size());
		}
		if(maxBins == minBins) maxBins++;
	}
	public void writeFile(ArrayList<String> writer)
	{
		BufferedWriter b;
		try 
		{
			FileWriter ff;
			if(pics)
				ff = new FileWriter(outDir+File.separator+"lorenz_analysis.txt",false);  
			else
				ff = new FileWriter("lorenz_analysis.txt",false);
			b = new BufferedWriter(ff);
			PrintWriter printer = new PrintWriter(b);
			for(int i = 0; i< writer.size(); i++)
			{
				printer.println(writer.get(i));
			}
			printer.close();
		}
		catch (IOException e) 
		{
			e.printStackTrace();
		}
	}
}
