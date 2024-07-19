package som_analysis;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import som.Map;

public class SingleReader 
{ 
	String[][] som; //[Node][ArrayList coord]  Will be uneven? Needs to be ArrayList?
	ArrayList<String[]> nodes; //each node's genes are stored here as a String[]
	int somSize;
	int xNodes, yNodes;
	int[] xCoords;
	int[] yCoords;
	public int[][] geneIndicesByShell;
	Map map;
	public boolean notFound, madeShell, printShells, printAll, printNotFounds, doNotPrint;
	public ArrayList<String> genio;
	
	public SingleReader(boolean all)
	{

		printNotFounds = false;
		somSize = 0;
		nodes = new ArrayList<String[]>();
		notFound = false;
		madeShell = false;
		
		genio = new ArrayList<String>();
		printShells = true;
		printAll = all;
		doNotPrint = true;
	}
	
	public void readSOM(String file)
	{
		BufferedReader br = null;
		String line = "";
		
		try 
		{
			br = new BufferedReader(new FileReader(file));
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	    try 
	    {
	    	int i = 0;
	    	int lineNum = 0;
	    	
	    	String line1 = br.readLine();
	    	yNodes = Integer.parseInt(line1.substring(line1.indexOf("x")+1));
	    	xNodes = Integer.parseInt(line1.substring(0,line1.indexOf("x")));
	    	somSize = xNodes * yNodes;
//	    	System.out.println(somSize);
//	    	System.out.println(line1);
	    	br.readLine(); //reads off the sigma information line, not used.
	    	xCoords = new int[somSize];
    		yCoords = new int[somSize];
	    	while (line != null)
	        {
	    		line = br.readLine();
	    		String[] gener = new String[0];
	    		boolean realLine = true;
	    		if(line !=null && line.contains(","))
	    		{
	    			/**gets the coordinate*/
	    			xCoords[i] = Integer.parseInt(line.substring(line.indexOf("[")+1, line.indexOf(",")));
	    			//System.out.println(xCoords[i]);
	    			yCoords[i] = Integer.parseInt(line.substring(line.indexOf(",")+1, line.indexOf("]")));
	    			//System.out.println(yCoords[i]);
	    			
	    			/** gets a list of the genes*/
	    			String geneList = line.substring(line.indexOf("(")+1,line.indexOf(")"));
	    			//System.out.println(geneList);
		    		gener = geneList.split(", ");
		    		//System.out.println(gener[0]);
		    		realLine = true;
	    		}
	    		else
	    			realLine = false;
	    		
	    		if(realLine)
	    		{
	    			nodes.add(gener);
	    			//System.out.println(genes.get(i)[0]);
	    			i++;
	    		}
	        }
	    	//System.out.println(genes.size()+"     "+ i);
	    } catch (IOException e)
	    {
			e.printStackTrace();
		} finally 
		{if(br!=null)
			{
	        	try {
	        		br.close();
	        	} catch (IOException e) 
				{
				e.printStackTrace();
				}
			}
		}
	}
	public void mapMaker()
	{
		map = new Map(xNodes, yNodes, false, 0);
		map.initializeMap();
//		int dd = map.hexDist(xCoords[0], yCoords[0], xCoords[3], yCoords[3]);
//		System.out.println(dd);
	}
	public void geneCheck(String genie, int shellLim)
	{
//		System.out.println(genie);
//		System.out.println(genes.get(1304)[2]);
		int geneShellIndex = -1;
		for(int i = 0; i<nodes.size(); i++)
		{
			String[] checking = nodes.get(i);
			for(int j = 0; j<checking.length;j++)
			{
				if(checking[j].equalsIgnoreCase(genie))
				{
					geneShellIndex = i;
					break;
				}
			}
		}
		if(geneShellIndex == -1)
		{
			notFound = true;
			if(printNotFounds)
				System.out.println("not found: " + genie);
		}
		else
		{
			geneIndicesByShell = new int[shellLim][shellLim*6];
			if(printShells)
				System.out.println(genie + " is at ("+ xCoords[geneShellIndex] + ", " +yCoords[geneShellIndex] + "); Index " + geneShellIndex);
			int[] currentShell = new int[1];
			currentShell[0] = geneShellIndex;
			geneIndicesByShell[0] = currentShell;
			
			for (int k = 1; k<shellLim; k++) 
			{
				currentShell = new int[k*6];
				int currentShellCount = 0;
				for(int i = 0; i<nodes.size(); i++)
				{
					if(map.neighborly[i][geneShellIndex] == k)
					{
						currentShell[currentShellCount] = i;
						//System.out.println(i);
						currentShellCount++;
					}
				}
				geneIndicesByShell[k] = currentShell;
			}
			
			
			//System.out.println("shell analysis for " + genie);
			consoleWriter();
			madeShell = true;
		}
	}
	public void consoleWriter()
	{
		for (int i = 0; i < geneIndicesByShell.length; i++)
		{
			if(printShells)
				System.out.println("shell " + i);
			for(int j = 0; j < geneIndicesByShell[i].length; j++)
			{
				//String[] generator = genes.get(shells[i][j]);
				//System.out.print(shells[i][j]);
				
				for(int k = 0; k < nodes.get(geneIndicesByShell[i][j]).length; k++)
				{
					genio.add(nodes.get(geneIndicesByShell[i][j])[k]);
					if(printAll && !doNotPrint)
						System.out.print(nodes.get(geneIndicesByShell[i][j])[k] + "\n");
					//System.out.print(generator[k] = ", ");
				}
				if(printShells)
					System.out.println();
			}
			//System.out.println();
		}
	}
	public String[][] indexGenesByDistance() /** returns a 2d string matrix with column one being gene names, column two being the string of the integer distance of that gene*/
	{
		ArrayList<String> geneAdder = new ArrayList<String>();
		ArrayList<Integer> distances = new ArrayList<Integer>();
		if(madeShell == false)
			return null;
		else
		{	
			for (int i = 0; i < geneIndicesByShell.length; i++)
			{
				int distance = i;
				for(int j = 0; j < geneIndicesByShell[i].length; j++)
				{
					for(int k = 0; k < nodes.get(geneIndicesByShell[i][j]).length; k++)
					{
						String gene = nodes.get(geneIndicesByShell[i][j])[k];
						geneAdder.add(gene);
						distances.add(distance);
						//System.out.print(gene + ", ");
						//System.out.print(generator[k] = ", ");
					}
					//System.out.println();
				}
				//System.out.println();
			}
			System.out.println(geneAdder.size() + " " +distances.size());
			String[][] geneAndDist = new String[geneAdder.size()][2];
			for(int i = 0; i< geneAdder.size(); i++)
			{
				geneAndDist[i][0] = geneAdder.get(i);
				geneAndDist[i][1] = Integer.toString(distances.get(i));
			}
			return geneAndDist;
		}	
	}
	public ArrayList<String> countMultipleRepeatedGenes(int n, boolean lined)
	{
		HashMap<String, Integer> entryCountMap = new HashMap<>();
		ArrayList<String> includedGenes = new ArrayList<String>();
		for (String entry : genio) 
		{
            // If the entry is already in the map, increment its count
            if (entryCountMap.containsKey(entry)) {
                int count = entryCountMap.get(entry);
                entryCountMap.put(entry, count + 1);
            } else 
            {
                // If the entry is not in the map, add it with a count of 1
                entryCountMap.put(entry, 1);
            }
		}
		for (String entry : entryCountMap.keySet()) 
		{
            int count = entryCountMap.get(entry);
            if(count>n)
            {
            	if(doNotPrint)
            	{
            		for(int i = 0; i< count; i++)
            			includedGenes.add(entry);
            	}
            	else
            	{
	            	if(lined)
	            		System.out.println(entry + "\t" + count);
	            	else
	            		System.out.print(entry + ", ");
            	}
            }
        }
		return includedGenes;
	}
	public void setDoNotPrint(boolean lean)
	{
		doNotPrint = lean;
	}
	public static void main(String[] args)
	{
		int requiredCalls = 1;
		int start = 3;
		try {
            // Attempt to convert the string to an integer
			requiredCalls = Integer.parseInt(args[2]);

            // If the conversion is successful, you can use the 'number' variable here
            System.out.println();
        } catch (NumberFormatException e) {
            System.err.println("No required calls threshold detected, set to 1 by default");
            start = 2;
        }
		
		boolean lined = false;
		SingleReader r = new SingleReader(requiredCalls <= 1);
		r.readSOM(args[0]);
		r.mapMaker();
		r.setDoNotPrint(false);
		int shell = 0;
		if(args.length>0)
			shell = Integer.parseInt(args[1]);
		//remaining args should be a list of genes
		for(int i = start; i< args.length; i++)
		{
			r.geneCheck(args[i], shell);
			if(lined)
				System.out.println("\n");
		}
		r.countMultipleRepeatedGenes(requiredCalls, lined);
	}
}
