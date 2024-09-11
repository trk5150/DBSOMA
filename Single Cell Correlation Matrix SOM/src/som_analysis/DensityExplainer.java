package som_analysis;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import mapScanning.DataPoint;
import mapScanning.DrawHex;
import mapScanning.MiniNode;
import mapScanning.SOMViewer;

public class DensityExplainer 
{
	public int minIndex;
	public ArrayList<int[]> combos;
	public int[][] dists;
	public String[] namesBest;
	public int[] indexBest;
	public boolean differentCenterFile;
	public int centers;
	public String SOMfile, densityFile, centerFile;
	public DrawHex som;
	public ArrayList<DataPoint> densityGenes, centerGenes;
	public DensityExplainer(int i, String file1, String file2, String file3)
	{
		differentCenterFile = true;
		centers = i;
		SOMfile = file1;
		densityFile = file2;
		centerFile = file3;
		som = new DrawHex(file1);
	}
	public DensityExplainer(int i, String file1, String file2)
	{
		differentCenterFile = false;
		centers = i;
		SOMfile = file1;
		densityFile = file2;
		centerFile = file2;
		som = new DrawHex(file1);
	}
	public void readFiles()
	{
		combos = new ArrayList<int[]>();
		readMap();
		densityGenes = new ArrayList<DataPoint>();
		readDensities();
		centerGenes = new ArrayList<DataPoint>();
		readCenters();
		//System.out.println(som.dataPoints.size() + " " + densityGenes.length);
	}
	public void readMap()
	{
		som.reader(SOMfile);
	}
	public void readDensities()
	{
		BufferedReader br = null;
		
		String line = "";
		try 
		{
			br = new BufferedReader(new FileReader(densityFile));
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	    try {
            line = br.readLine();

        	int tryCount=0;
        	int foundCount=0;
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
	        			tryCount++;
		        		boolean found = false;
			        	for(int i = 0; i< som.dataPoints.size(); i++)
			        	{
			        		if(som.dataPoints.get(i).name.equalsIgnoreCase(s[k].trim()))
			        		{
			        			found = true;
			        			foundCount++;
			        			densityGenes.add(som.dataPoints.get(i));
			        			//System.out.println(s[k].trim());
			        		}
			        		
			        	}
//			        	if(found == false)
//		        		{
//		        			System.out.println(s[k] + " not found");
//		        		}
		        	}
	        	}
	        	line = br.readLine();
        	}
            System.out.println("searched for " + tryCount + " genes, found "+foundCount);   
	    } catch (IOException e)
	    {
			e.printStackTrace();
		}
	    		
	}
	public void readCenters()
	{
		if(differentCenterFile == false)
			centerGenes = densityGenes;
		else
		{
			ArrayList<String> g = new ArrayList<String>();
			BufferedReader br = null;
			
			String line = "";
			try 
			{
				br = new BufferedReader(new FileReader(centerFile));
				
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		    try {
	            line = br.readLine();

	        	int tryCount=0;
	        	int foundCount=0;
	            while (line != null) 
		        {
		        	String[] s = new String[1];
		        	s[0] = line;
		        	if(line.contains(","))
		        	{
		        		s = line.split(",");
		        	}
//		        	if(s.length>0)
//		        		System.out.println(s[0] + " " + s.length);
		        	for(int k = 0; k<s.length; k++)
		        	{
		        		if(s[k].length()>1)
		        		{
		        			tryCount++;
			        		boolean found = false;
				        	for(int i = 0; i< som.dataPoints.size(); i++)
				        	{
				        		if(som.dataPoints.get(i).name.equalsIgnoreCase(s[k].trim()))
				        		{
				        			found = true;
				        			foundCount++;
				        			centerGenes.add(som.dataPoints.get(i));
				        			//System.out.println(s[k].trim());
				        		}
				        		
				        	}
//				        	if(found == false)
//			        		{
//			        			System.out.println(s[k] + " not found");
//			        		}
			        	}
		        	}
		        	line = br.readLine();
	        	}
	            System.out.println("searched for " + tryCount + " genes, found "+foundCount);   
		    } catch (IOException e)
		    {
				e.printStackTrace();
			}
		}
	}
	public void distanceToCenters()
	{
		dists = new int[densityGenes.size()][centerGenes.size()];
		for(int i = 0; i<densityGenes.size(); i++)
		{
			for(int j = 0; j<centerGenes.size(); j++)
			{
				int x1 = densityGenes.get(i).myMini.xCoord;
				int y1 = densityGenes.get(i).myMini.yCoord;
				int x2 = centerGenes.get(j).myMini.xCoord;
				int y2 = centerGenes.get(j).myMini.yCoord;
				int ddd = som.hexDist(x1, y1, x2, y2);
//				System.out.println(ddd);
				dists[i][j] = ddd;
			}
		}

		//System.out.println(som.hexDist(0,0,49,49));
	}
	public void densityAnalysisPrep()
	{
		for(int counter = 0; counter < 1000; counter++)
		{
			int[] kk = new int[centers];
			int i = 0;
			while(i<centers)
			{
				kk[i] = (int)(Math.random()*centerGenes.size());
				//System.out.println(kk[i]);
				i++;
			}
			combos.add(kk);
		}
		explainDensity();
	}
	public void explainDensity()
	{
		int[] totalMinDist = new int[combos.size()];
		//for each combination:
		for(int i = 0; i< combos.size(); i++)
		{
			//cycle through each gene
			for(int j = 0; j < densityGenes.size(); j++)
			{
				int tempMin = 1000000000;//dists[j][0]+25;  //just needs to be a "max" type value
				//cycle through each TF in the combination for it's distance to that gene
				for(int k = 0; k < combos.get(i).length; k++)
				{
					if(dists[j][combos.get(i)[k]]<tempMin)
					{
						//System.out.println("tempMin = " +tempMin +" and new value = " + dists[j][combos.get(i)[k]]);
						tempMin = dists[j][combos.get(i)[k]];
					}
				}
				totalMinDist[i] += tempMin;
			}
		}
		int min = totalMinDist[0];
		minIndex = 0;
		for(int i = 0; i<totalMinDist.length; i++)
		{
			if(totalMinDist[i] < min)
			{
				minIndex = i;
				min = totalMinDist[i];
			}
			//System.out.println(totalMinDist[i] + " at "+ i);
		}
		//System.out.print(min + " is lowest, at index " + minIndex +  " includes: ");
		namesBest = new String[combos.get(minIndex).length];
		indexBest = new int[combos.get(minIndex).length];
		
		for(int kk = 0; kk < combos.get(minIndex).length; kk++)
		{
			String n = centerGenes.get(combos.get(minIndex)[kk]).name;
			namesBest[kk] = n;
			indexBest[kk] = combos.get(minIndex)[kk];
			//System.out.println(n);
			System.out.print(n + "\t");
			//System.out.println(namesBest[kk] + " "+ indexBest[kk] + " \t");
		}
		//System.out.println(min);
		System.out.println();
	}
	public void amongTheTop(double portion)
	{
		
	}
	public void printGenesAssigned()
	{
		System.out.println();
		for(int k = 0; k< indexBest.length; k++) //each center
		{
			System.out.println(namesBest[k] + " is center for: ");
			for(int i = 0; i < dists[k].length; i++) //check each gene
			{
				
				int dd = dists[i][combos.get(minIndex)[k]];
				int minDist = dd;
				for(int j = 0; j< indexBest.length; j++)
				{
					int distToCenter = dists[i][combos.get(minIndex)[j]];
					if(distToCenter<minDist)
						minDist= distToCenter;
				}
				if(dd<=minDist)
				{
					//System.out.println(dists[k][i]);
					System.out.print(densityGenes.get(i).name + ", ");
				}
			}
			System.out.println("\n");
		}
		
	}
	//takes: number of density centers (int), file location for a SOM, file location for a list of genes (defines a density), optional 3rd file location for specific density centers (otherwise uses genes from density file)
	public static void main(String[] args)
	{
		boolean many = false; 
		
		DensityExplainer d;
		if(args.length == 4)
			d = new DensityExplainer(Integer.parseInt(args[0]), args[1], args[2], args[3]);
		else
			d= new DensityExplainer(Integer.parseInt(args[0]), args[1], args[2]);
		d.readFiles();
		d.distanceToCenters();
		int kk = d.centers;
		if(many)
		{
		//for(int j = 0; j<kk-2; j++)
			{
				//System.out.println(d.centers);
				for(int i = 0; i < 200; i++)
				{
					d.densityAnalysisPrep();
					d.combos.clear();
				}
			//	d.centers -= 1;
			}
		}
		else
		{
			d.densityAnalysisPrep();
			d.printGenesAssigned();
		}
	}
}