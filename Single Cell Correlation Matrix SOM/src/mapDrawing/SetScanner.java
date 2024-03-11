package mapDrawing;

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


import java.util.HashMap;

/**
 * All searching is done as follows:
 * 1) For each file an empty array of integers is created
 * 		Each row: Genes called in a list, gene counts from DBscan
 * 		There is a separate weight matrix calculated as sum of genes in each node
 * 2) gene lists is run through the hash map, 1 gene at a time
 * 3) The genes called values are populated by the node index in the hashmap
 * 4) the DBscan method uses the nodesWithinRadius method to increment the number of calls that node has received
 * 5) The DBscan column of the array can then be used to compare to state files, which have received same treatment as above
 * */

public class SetScanner 
{
	public boolean saveImages = true;
	public boolean estimateFalsePositiveRate = true;
	public int falsePositiveTrials = 100000;
	//Should change to asking user for approximate avg size of input files, estimating the false positive rate, then basing cutoffs off that.
	public double structureCutOff = 0.067; //~2x max structure from randoms (0.4)
	public double qualityCutOff = 0.05; //~2x max quality from randoms (0.17)
	public double stateOverlapCutOff = 0.2; //Somewhere to start?
	public boolean writeFalsePositiveTests = false;
	public boolean useNodeWeightAtCuttoff = false; //used in the boolean comparison method (+1 for node overlap vs +nodeweight)
	public int roundToInt = 100000;
	public int radius;
	public int minPts;
	public double threshold, maxx, maxq, avgFileNum;
	public String[] geneStates;
	public int nodes,xNodes,yNodes, filesOverThresholds;
	public File outDir, trimmedOutDir;
	public ArrayList<File> searchers;
	public String fold, mapper;
	public ArrayList<DataPoint> bins;
	public MiniSystem nodeSystem;
	public int[][] nodesWithinRadius;
	public int[] nodeWeights;
	public HashMap<String, Integer> geneToNode;
	public ArrayList<String> mappedGenes;
	
	public boolean debug = true;

	public ArrayList<MiniNode> nodeList;
	
	public SetScanner(String map, String fillle, String outdir, String[] states)
	{
		radius = 1;
		minPts = 5;
		
		geneStates = states;
		mapper = map;
		fold = fillle;
		File folder = new File(fillle);
		outDir = new File(outdir);
		trimmedOutDir = new File(outdir + "\\Trimmed");
		String mapp = map;
		yNodes=0;
		xNodes=0;
		
		mapReader(mapp);

		nodes = yNodes*xNodes;
		geneToNode = new HashMap<String, Integer>();
		hashMapper();
		populateDistanceMatrix();
		findSearchFiles(folder);
	}
	public void findSearchFiles(File folder)
	{
		searchers = new ArrayList<File>();
		try {
			System.out.println("Analyzing files in "+folder.getCanonicalPath());
		} catch (IOException e) {
			e.printStackTrace();
		}
		File[] listOfFiles = folder.listFiles();

		int count = 0;
		for (File file : listOfFiles) 
		{
			count++;
			if(debug && count%100000 == 0)
			{
				System.out.println("scanned files count " + count);
			}
		    if (file.isFile() && file.getName().endsWith(".grp") || file.getName().endsWith(".txt")) 
		    {
		    	searchers.add(file);
		    }
		}
	}
	public void hashMapper()
	{
		mappedGenes = new ArrayList<String>();
		nodeWeights = new int[nodes];
		for(int i = 0; i <nodes; i++)
		{
//			System.out.println(nodeList.get(i).bins.size());
			for(int j = 0; j < nodeList.get(i).bins.size(); j++)
			{
				String s = nodeList.get(i).bins.get(j).name;
				mappedGenes.add(s);
				nodeWeights[i]++;
				geneToNode.put(s, i);
//				System.out.println(nodeList.get(i).bins.get(j).name);
				
			}
		}
		if(debug)
		{
			System.out.println("Gene-to-node hashmap created");
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
	}
	public void populateDistanceMatrix()
	{
		int nodesWithinNum = 0;
		int increment = 6;
		for(int i = 0; i<radius; i++)
		{
			nodesWithinNum += increment; 
			increment += 6;
		}
		if(debug)
			System.out.println("Nodes within radius " + radius + " = " + nodesWithinNum);
		
		double[][] allDists = degreesOfSep();
		nodesWithinRadius = new int[nodes][nodesWithinNum];
		for(int i = 0; i< nodes-1; i++)
		{
			int count = 0;
			for(int j = 0; j < nodes; j++)
			{
				if(i!=j && allDists[i][j] <=radius)
				{
					nodesWithinRadius[i][count] = j;
					count++;
				}
			}

			if(count>nodesWithinNum)
				System.err.println("too many nodes within radius, something has gone wrong");	
		}
		
		if(debug)
		{
			int locus = 220;

			System.out.print("Sample within radius " + radius + " for node " + locus + ": ");
			for(int i = 0; i< nodesWithinRadius[locus].length; i++)
			{
				System.out.print(nodesWithinRadius[locus][i]+", ");
			}
			System.out.println();
		}
		
	}
	public double[][] degreesOfSep()
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
		return dist;		
	}
	public int hexDist(int x1, int y1, int x2, int y2)
	{
		//System.out.println(xNodes + "  " + yNodes);
		
		if(x1 == x2 && y1 == y2)
			return 0;
		else
		{
			int hm = 0;
			int vm = 0;
			int cm = 0;
			boolean right = x1-x2<0;
			
			hm = Math.abs(x1-x2);
			vm = Math.abs(y1-y2);
			if(xNodes-x1<xNodes-x2 && xNodes-x1 + x2 < hm)
			{
				right = true;
				hm =  xNodes-x1 + x2;	
			}
			else if (xNodes-x1>xNodes-x2 && xNodes-x2 +x1 < hm)
			{
				right = false;
				hm = xNodes-x2 +x1;
			}
			
			if(yNodes-y1<yNodes-y2 && yNodes-y1 + y2 < vm)
			{
				vm = yNodes-y1 + y2;
			}
			else if(yNodes-y1>yNodes-y2 && yNodes-y2 + y1 < vm)
			{
				vm = yNodes-y2 + y1;
			}
			
			if(vm %2 !=0)
			{
				if(!right && y1%2 == 0)
				{
					cm = (int) Math.min(hm, (((double)vm)/2 +.5));
				}
				else if(!right && y1%2 == 1)
				{
					cm = (int) Math.min(hm, (((double)vm)/2));
				}
				else if(right && y1%2 == 0)
				{
					cm = (int) Math.min(hm, (((double)vm)/2));
				}
				else //if(right && y1%2 == 1)
				{
					cm = (int) Math.min(hm, (((double)vm)/2 +.5));
				}
			}
			else
				cm = (int) Math.min(hm, ((((double)vm)/2)));
			return hm+vm-cm; 
		}
	}
	public void searchSystem()
	{
		ArrayList<int[][]> forPainting = new ArrayList<int[][]>();
		ArrayList<File> overlappingFiles = new ArrayList<File>();
//		ArrayList<File> structuredFiles = new ArrayList<File>();
		ArrayList<String> writer = new ArrayList<String>();
		ArrayList<String> overlapWriter = new ArrayList<String>();
		ArrayList<String> structureWriter = new ArrayList<String>();
		ArrayList<File> stateFiles = new ArrayList<File>();
		String stately = "";
		for(int i = 0; i < geneStates.length; i++)
		{
			stately = stately.concat("\tpercent overlap with " + geneStates[i]);
		}
		String headers = "fileName" + "\t" + "found genes" + "\t" + "dbScan nodes" +"\t" + "Scan:found ratio" + "\t" + "quality" +stately;
		writer.add(headers);
		overlapWriter.add(headers);
		
		//System.out.println(stately);
		
		
		//Produces ArrayList of state search/scan arrays
		ArrayList<int[][]> stateLists = new ArrayList<int[][]>();
		ArrayList<boolean[]> boolStates = new ArrayList<boolean[]>();
		for(String fifi : geneStates)
		{
			File fifo = new File(fifi);
			stateFiles.add(fifo);
//			System.out.println(fifo.getName());
			int[][] scanned = searchFile(fifo);
			stateLists.add(scanned);
			boolStates.add(booleanize(scanned));
		}
		
		//Looks through all search files and creates search/scan arrays for each
		int fileCount = 0;
		int totalFoundCount = 0;
		for (File s:searchers)
		{
			if(debug && fileCount%10000==0)
				System.out.println("file count = " + fileCount);
			fileCount++;
			
			
			int[][] scanned = searchFile(s);
			
			int sum = sumArrays(scanned);
			boolean[] boolList = booleanize(scanned);
			int boolTot = boolSum(boolList);
			double[] overlaps = compareToStates(boolList, boolStates);
			String name = s.getName();
			totalFoundCount+= sum;
			
			double boolVSum = roundDouble((double)boolTot/(double)sum);
//			System.out.println(boolVSum + ", " + boolTot + ", " + sum + ", " + overlaps[0]);
			
			

			double quality = callsInScannedNodesRatio(boolList, scanned);
			String written = name + "\t" + sum +"\t"+ boolTot + "\t" +boolVSum + "\t" + quality;
			for(double d : overlaps)
			{
				written = written.concat("\t" +roundDouble(d));				
			}
			
			boolean overlapped = false;
			boolean structured = false;
			
			if(boolVSum > structureCutOff && quality>qualityCutOff)
			{
				structured = true;
				if(overlaps.length>0 &&overlaps[0]>stateOverlapCutOff)
					overlapped = true;
			}
			
			if(structured)
			{
//				structuredFiles.add(s);
				structureWriter.add(written);
				
				if(overlapped)
				{
					forPainting.add(scanned);
					overlappingFiles.add(s);
					overlapWriter.add(written);
				}
			}	
			writer.add(written);
		}
		avgFileNum = totalFoundCount/fileCount;
		if(debug)
			System.out.println(avgFileNum);
		
		if(estimateFalsePositiveRate)
		{
			double maxStructure = 0;
			double maxQuality = 0;
//			double[] maxOverlaps = new double[geneStates.length];
			System.out.println("estimating false positives");
			
			//produces a given number of random sets of genes then calculates the structure and overlap of each
			for(int i =0; i < falsePositiveTrials; i++)
			{
				ArrayList<String> randomGenes = new ArrayList<String>();
				for(int j = 0; j < avgFileNum; j++)
				{
					randomGenes.add(mappedGenes.get((int) (mappedGenes.size()* Math.random())));
				}
				int[][] scanned = searchStrings(randomGenes);
				
				
				int sum = sumArrays(scanned);
				boolean[] boolList = booleanize(scanned);
				int boolTot = boolSum(boolList);
				double[] overlaps = compareToStates(boolList, boolStates);
				String name = "random set number " + i;
				double boolVSum = roundDouble((double)boolTot/(double)sum);
				double quality = callsInScannedNodesRatio(boolList, scanned);
				
				if (boolVSum >maxStructure)				
					maxStructure = boolVSum;
					
				if(quality>maxQuality)
					maxQuality = quality;
				if(writeFalsePositiveTests&&boolTot>0)
				{
					boolVSum = roundDouble(boolVSum);
					String written = name + "\t" + sum +"\t"+ boolTot + "\t" + boolVSum + "\t" + quality;
					for(double d : overlaps)
					{
						written = written.concat("\t" +roundDouble(d));
					}
//					System.out.println(written);
					writer.add(written);
				}
			}
		

			maxx = maxStructure;
			maxq = maxQuality;
			System.out.println("Max Structure = " + maxStructure + "; Max quality = " + maxQuality);
		}
		
		writeSelectFiles(overlapWriter, "hits");
		writeSelectFiles(structureWriter, "structured");
		writeFile(writer);
		writeInfoFile();
		System.out.println("Files meeting thresholds: " + overlappingFiles.size());
		if(saveImages)
		{
			goodFilePaint(stateLists, stateFiles, true);
			goodFilePaint(forPainting, overlappingFiles, false);
		}
//			goodFiles(selectFiles);
//		System.out.println(fileCount);
//		writeFile(writer);
		System.out.println("done");
	}
	public double callsInScannedNodesRatio(boolean[] bools, int[][] calls)
	{
		if(bools.length != calls.length)
			return -1;
		
		double assignedToAValidNode = 0;
		double assignedToAnyNode = 0;
		
		for(int i = 0; i < bools.length; i++)
		{
			if(calls[i][0]>0)
			{
				assignedToAnyNode++;
				if(bools[i])
					assignedToAValidNode++;
			}
			
		}
		return assignedToAValidNode/assignedToAnyNode;
	}
			
	public double roundDouble(double d)
	{
		double rounded = d*roundToInt;
		int ri = (int)rounded;
		rounded = ((double)ri)/roundToInt;
		return rounded;
	}
	
	public void goodFilePaint(ArrayList<int[][]> s, ArrayList<File> names, boolean stateFiles)
	{
		PaintSetScan painter = new PaintSetScan(s, minPts, xNodes, yNodes, nodeWeights, outDir, names);
		painter.searchAll(s, stateFiles);
	}

	public double[] compareToStates(boolean[] comp, ArrayList<boolean[]> states)
	{
		double[] overlaps = new double[states.size()];
		
		for(int i = 0; i< states.size(); i++)
		{
			overlaps[i] = compareToState(comp, states.get(i));
		}
		
		return overlaps;
	}
	public double compareToState(boolean[] comp, boolean[] state)
	{
		double overlap = 0;
		double compCount = 0;
		if(comp.length != state.length)
			System.err.println("comparison lengths are not equal, matrix handling is likely broken");
		for(int i = 0; i < comp.length; i++)
		{
			if(useNodeWeightAtCuttoff)
			{
				if(comp[i])
					compCount+=nodeWeights[i];
				if(comp[i] && state[i])
					overlap += nodeWeights[i];
			}
			else
			{
				if(comp[i])
					compCount++;
				if(comp[i] && state[i])
				{
					overlap ++;
				}
			}
			
		}
		overlap = overlap/compCount;
//		System.out.println(overlap);
		return overlap;
	}
	public boolean[] booleanize(int[][] scanned)
	{
		boolean[] b = new boolean[scanned.length];
		for(int i = 0; i < b.length; i++)
		{
			b[i] = scanned[i][1] > minPts;
		}
		return b;
	}
	public int boolSum(boolean[] b)
	{
		int sum = 0;
		for(boolean bb: b)
		{
			if(bb)
				sum++;
		}
		return sum;
	}
	public int sumArrays(int[][] s)
	{
		int sums = 0;
		for(int i = 0; i < s.length; i++)
		{
			sums += s[i][0];
		}
		return sums;
	}
	public int[][] searchFile(File f)
	{
		int [][] counts = new int[nodes][2]; //found genes in the first value, dbscan values in the second value;
		for(String s : parseGenesFromFile(f)) //parseGenesFromFile method already filters out non-found genes
		{
			counts[geneToNode.get(s)][0]++;
		}
		int [][] scanned = dbScan(counts);
		return scanned;
	}
	public int[][] searchStrings(ArrayList<String> names)
	{
		int [][] counts = new int[nodes][2];
		for(String s : names)
		{
			counts[geneToNode.get(s)][0]++;
		}
		int [][] scanned = dbScan(counts);
		return scanned;
	}
	public int[][] dbScan(int[][] input)
	{
		//each node checks whether it has value, then, if it does, adds its value to the DBscan column of each node within the given radius
		for(int i = 0; i < nodes; i++)
		{
			int value = input[i][0];
			if(value >0)
			{
				input[i][1] += value;
				for(int j = 0; j < nodesWithinRadius[i].length; j++)
				{
					int index = nodesWithinRadius[i][j];
					input[index][1]+=value;
				}
			}
		}
		return input;
	}
	public ArrayList<String> parseGenesFromFile(File f)
	{
		BufferedReader br = null;
		String line = "";
		try 
		{
			br = new BufferedReader(new FileReader(f));
			
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
	        			if(geneToNode.containsKey(s[k].trim()))
	        			{
	        				foundGenes.add(s[k]);
	        			}
		        	}
	        	}
	        	
	        	
	        	line = br.readLine();
        	}
//            if(debug)
//            	System.out.println("searched for " + tryCount + " genes, found "+foundCount);   
		    } catch (IOException e)
		    {
				e.printStackTrace();
			}
	    return foundGenes;
	}
	public void writeFile(ArrayList<String> writer)
	{
		BufferedWriter b;
		try 
		{
			FileWriter ff;
			ff = new FileWriter(outDir+File.separator+"SOM_Analysis.txt",false);
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
	public void writeInfoFile()
	{
		BufferedWriter b;
		try 
		{
			FileWriter ff;
			ff = new FileWriter(outDir+File.separator+"Analysis_metadata.txt",false);
			b = new BufferedWriter(ff);
			PrintWriter printer = new PrintWriter(b);
			printer.println("Radius = " + radius);
			printer.println("minPts = " + minPts);
			printer.println("Map File: " + mapper);
			printer.println("Projection directory: " + fold);
			printer.println("Output Directory: " + outDir);
			printer.println("False positive number: " + falsePositiveTrials);
			printer.println("False positive structure = " + maxx);
			printer.println("False positive quality = " + maxq);
			printer.println("random set size (avg found genes) =" + avgFileNum);
			printer.println("using node weight for cutoff = " + useNodeWeightAtCuttoff);
			printer.println("structureCutOff =" + structureCutOff);
			printer.println("stateOverlapCutOff =" + stateOverlapCutOff);
			printer.println("Files meeting both thresholds = " + filesOverThresholds);
			
			
			
			printer.println();
			printer.close();
		}
		catch (IOException e) 
		{
			e.printStackTrace();
		}
	}
	public void writeSelectFiles(ArrayList<String> writer, String s)
	{
		BufferedWriter b;
		try 
		{
			FileWriter ff;
			ff = new FileWriter(outDir+File.separator+s+"_Analysis.txt",false);
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
