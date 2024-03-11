package correlationMatrixMaker;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

//reads in the common format for sparse matrices, which is a file with gene names, then a second file with coordinates+count (often comes along with a file of cell identifiers, which we ignore)
//matrix file is [gene][cell][count]

/**
 * To do:
 * 	ensure the genes file saves in agreement with the matrix
 * 		Do things known to correlate appear that way
 * 
 * 	Complete "trained from: " statements
 * */

public class SparseParse 
{
	public int roundToInt = 10000;
	public double threshold;
	public boolean[] inclusions;
	public int genes, cells, totalCount;
	public String matFile;
	public String genFile;
	public Correlator correlator;
	public SparseParse(String matrixFile, String genesFile)
	{
		threshold = 0.01;
		matFile = matrixFile;
		genFile = genesFile;
		correlator = new Correlator();
		inclusions = new boolean[0];
	}
	public ArrayList<int[]> readMat()
	{

    	ArrayList<int[]> readout = new ArrayList<int[]>();
		BufferedReader br = null;
		String line = "";
		
		try 
		{
			br = new BufferedReader(new FileReader(matFile));
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
	    try 
	    {
//	    	if(line.contains("%"))
	    		line=br.readLine();

	    	//reading in sizes
	    	line = br.readLine();
	    	String[] lines = line.split(" ");
	    	genes = Integer.parseInt(lines[0]);
	    	cells = Integer.parseInt(lines[1]);
	    	totalCount = Integer.parseInt(lines[2]);
	    	System.out.println("cells: "+cells);
	    	System.out.println("Genes in file header: "+genes);
	    	System.out.println("Total Counts in file: "+totalCount);
	    	
    		line = br.readLine();
	    	while(line != null)
	    	{
	    		lines = line.split(" ");
	    		int g = Integer.parseInt(lines[0]);
		    	int c = Integer.parseInt(lines[1]);
		    	int count = Integer.parseInt(lines[2]);
	    		int[] data = new int[3];
	    		data[0] = g;
	    		data[1] = c;
	    		data[2] = count;
	    		readout.add(data);
	    		line = br.readLine();
	    	}
	    	System.out.println("lines: " + readout.size());
	    	
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
	    return readout;
	}
	public ArrayList<String> readGenes()
	{
		ArrayList<String> names = new ArrayList<String>();
		BufferedReader br = null;
		String line = "";
		
		try 
		{
			br = new BufferedReader(new FileReader(genFile));
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		int count = 0;
	    try 
	    {
     		line = br.readLine();
 	    	while(line != null)
 	    	{
 	    		count++;
 	    		names.add(line.trim());
 	    		line = br.readLine();
 	    	}
 	    	
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
	    System.out.println("Genes file lines: " + names.size()+ " "  + count);
 	    return names;
	}
	public double roundDouble(double d)
	{
		double rounded = d*roundToInt;
		int ri = (int)rounded;
		rounded = ((double)ri)/roundToInt;
		return rounded;
	}
	public double[][] correlate(ArrayList<int[]> coords, ArrayList<int[]> counts, int[] totals)
	{
		double[][] cors = new double[coords.size()][coords.size()];
		double[] bars = new double[totals.length];
		
		//calc avg for each gene
		for(int i = 0; i< bars.length; i++)
			bars[i] = (double)totals[i]/(double)cells;
		
		
		
		double[] stdevs = new double[totals.length];
		for(int i = 0; i< stdevs.length; i++)
		{
			stdevs[i] = stdev(coords.get(i), counts.get(i), bars[i]);
//			System.out.println(stdevs[i]);
		}
		System.out.println("stdevs calculated");
		
		for(int i = 0; i<coords.size(); i++)
		{
			for (int j = i; j<coords.size(); j++)
			{
				if(i == j)
					cors[i][j] = 1;
				else
				{
					int[] aCoords = coords.get(i); 
					int[] bCoords = coords.get(j);
					int[] aMags = counts.get(i);
					int[] bMags = counts.get(j);
					double pears = pearsonFast(aCoords, bCoords, aMags, bMags, bars[i], bars[j], stdevs[i], stdevs[j]);
					double r = roundDouble(pears);
					cors[i][j] = r;
					cors[j][i] = r;
				}
			}
			if(i%250 == 0)
				System.out.println("count: " + i);
		}
		
		//Old, slower code
//		for(int i = 0; i<coords.size(); i++)
//		{
//			for (int j = i; j<coords.size(); j++)
//			{
//				if(i == j)
//					cors[i][j] = 1;
//				else
//				{
//					int aTotal = totals[i];
//					int bTotal = totals[j];
//					int[] aCoords = coords.get(i); 
//					int[] bCoords = coords.get(j);
//					int[] aMags = counts.get(i);
//					int[] bMags = counts.get(j);
//					double pears = pearson(aCoords, bCoords, aMags, bMags, aTotal, bTotal);
//					cors[i][j] = pears;
//					cors[j][i] = pears;
//				}
//			}
//			System.out.println("count: " + i);
//		}
		return cors;
	}
	
	public double stdev(int[] xCoords, int[] xCounts, double xBar)
	{
		double n = cells;
		
		
		double sX = 0;
		
		sX = (n-xCoords.length)*(xBar*xBar); //equal to the sum of (0-values minus average)^2
		for(int i = 0; i < xCoords.length; i++)
		{
			sX+= ((xCounts[i]-xBar)*(xCounts[i]-xBar));
		}
//		System.out.println(sX);
		sX /= (n-1);
		sX = Math.sqrt(sX);
		
		return sX;
	}
	public double pearson(int[] xCoords, int[] yCoords, int[] xCounts, int[] yCounts, int xTotal, int yTotal)
	{
		double r = 0;
		if(xTotal == 0 || yTotal == 0)
			System.out.println("A zero total survived pruning: " + xTotal + ", " +yTotal);
		
		int n = cells;
		double xBar = (double)xTotal/(double)n;
		double yBar = (double)yTotal/(double)n;
		
		double sX = 0;
		sX = (n-xCoords.length)*(xBar*xBar); //equal to the sum of (0-values minus average)^2
		for(int i = 0; i < xCoords.length; i++)
		{
			sX+= ((xCounts[i]-xBar)*(xCounts[i]-xBar));
		}
		sX /= (n-1);
		sX = Math.sqrt(sX);
		
		double sY = 0;
		sY = (n-yCoords.length)*yBar*yBar; 
		for(int i = 0; i < yCoords.length; i++)
		{
			sY+= ((yCounts[i]-yBar)*(yCounts[i]-yBar));
		}
		sY /= (n-1);
		sY = Math.sqrt(sY);
		
		double num = 0;
		for(int i = 0 ; i < xCoords.length; i++)
		{
			int compare = xCoords[i];
			for(int j = 0; j < yCoords.length; j++)
			{
				if(yCoords[j] == compare)
				{
					num += xCounts[i] * yCounts[j];
				}
			}
		}
		num -= n*xBar*yBar;
		
		double denom = 0;
		denom = (n-1)*sX*sY;
		
		r = num/denom;
//		
//		System.out.println(r);
		return r;
	}
	public double pearsonFast(int[] xCoords, int[] yCoords, int[] xCounts, int[] yCounts, double xBars, double yBars, double sX, double sY)
	{
		double r = 0;
		
		double n = cells;
		double num = 0;

		int lastJ = 0;
    	for(int i = 0 ; i < xCoords.length; i++)
		{
			int compare = xCoords[i];
			for(int j = lastJ; j < yCoords.length; j++)
			{
				if(yCoords[j] == compare)
				{
					lastJ = j; //coordinates are listing small to large, so next search can start from the previous find
					num += xCounts[i] * yCounts[j];
					break; //the above boolean operation can only ever occur once per inner loop, break if found
				}
				if(yCoords[j]> compare)
					break; //if yCoords[j] has already overshot xCoords[i], no reason to increment upwards as values will only increase
			}
		}
		num -= n*xBars*yBars;
		double denom = (n-1)*sX*sY;
		r = num/denom;
//		System.out.println(r);
		return r;
	}
	public void writeFiles(double[][] corMatrix, ArrayList<String> genio)
	{		
		BufferedWriter b;
		File f = new File(matFile);
		String s = f.getParent();
		s = s +"\\";
		String corFilename = s+"_Correlations.txt";
		String geneFilename = s+"_Genes.txt";
		//Correlations file
		try 
		{
			System.out.println("Printing correlations to "+corFilename);
			FileWriter ff = new FileWriter(corFilename,true);
			b = new BufferedWriter(ff);
			PrintWriter printer = new PrintWriter(b);
			printer.print("trained from: ");

			for(int i = 0; i<corMatrix.length; i++)
			{
				printer.print("\n");
				for(int j = 0; j<corMatrix[i].length; j++)
				{
					if(Double.isNaN(corMatrix[i][j]))
					{
						System.err.println("NaN detected at line i=" + i + ";  position j =" + j + "\t set to zero by default to allow clustering to attempt to continue. Recheck input matrix to ensure NaNs are not an issue");
	    				corMatrix[i][j]= 0;
					}
					printer.print(corMatrix[i][j]+"\t");
				}
			}
			printer.print("\n");
			printer.close();
		}catch (IOException e){
			e.printStackTrace();
			System.err.print("Error printing SOM file "+corFilename);
		}
		
		try 
		{
			System.out.println("Printing correlations to "+geneFilename);
			FileWriter ff = new FileWriter(geneFilename,true);
			b = new BufferedWriter(ff);
			PrintWriter printer = new PrintWriter(b);
			for(int i = 0; i<genio.size(); i++)
			{
				printer.print(genio.get(i));
				printer.print("\n");
			}
			printer.close();
		}catch (IOException e){
			e.printStackTrace();
			System.err.print("Error printing gene file "+geneFilename);
		}
		
	}
	public int[] getCoordsCalc(ArrayList<int[]> r)
	{
		int max = 0;
		int[] counts = new int[genes];
		for(int j = 0; j< r.size(); j++)
		{
			int[] cur = r.get(j);
			counts[cur[0]-1]++;
		}
		for(int i = 0; i < counts.length; i++)
		{
			if(counts[i]> max)
				max = counts[i];
		}
		System.out.println("Max coords" + max);
		return counts;
	}
	public double[][] populateVectors(ArrayList<int[]> r)
	{
		int[] cellCoordCounts  = getCoordsCalc(r);
		ArrayList<int[]> allCoords = new ArrayList<int[]>();
		ArrayList<int[]> allCounts = new ArrayList<int[]>();
		int[] totals = new int[genes];
		
		//create 2x jagged 2d array for gene x cell to store coordinates and their counts
		for(int i =0; i < cellCoordCounts.length; i++)
		{
			allCoords.add(new int[cellCoordCounts[i]]);
			allCounts.add(new int[cellCoordCounts[i]]);
		}
		
		int[] currentIndex = new int[genes];
		for(int i = 0; i < r.size(); i++)
		{
			int[] cur = r.get(i);
			int gene = cur[0]-1;
			int cell = cur[1]-1;
			int count = cur[2];
			
			int index = currentIndex[gene];
			currentIndex[gene]++;
			
			totals[gene] += count;
			allCoords.get(gene)[index] = cell;
			allCounts.get(gene)[index] = count;
		}
		int count = 0;
		double sum = 0;
		for(int i = 0; i < totals.length; i++)
		{
			sum+=totals[i];
		}
		double avg = sum/totals.length;
		double newThreshold = threshold*avg;
		
		boolean[] includer = new boolean[genes];
		for(int i = 0; i < totals.length; i++)
		{
			boolean bb = totals[i]>newThreshold;
			if(bb)
			{
				count++;
				includer[i] = bb;
			}
		}
		int[] prunedTotals = new int[count];
		System.out.println(count + " above threshold, " + 0.05*(sum/totals.length) +" avg");
		int index = 0;
		for(int i = 0; i < totals.length; i++)
		{
			if(includer[i])
			{
				if(totals[i] == 0)
					System.out.println("its broken");
				prunedTotals[index] = totals[i];
				index++;
			}
		}
		
		ArrayList<int[]> prunedCoords = pruneMatrix(allCoords, includer);
		ArrayList<int[]> prunedCounts = pruneMatrix(allCounts, includer);
		setInclusions(includer);
		
		
		System.out.println(prunedCoords.size()+" " + prunedCounts.size() + " " + prunedTotals.length);
		return correlate(prunedCoords, prunedCounts, prunedTotals);
		
	}
	public void setInclusions(boolean[] a)
	{
		inclusions = a;
	}
	public ArrayList<int[]> pruneMatrix(ArrayList<int[]> fixIt, boolean[] includer)
	{
		if(fixIt.size()!= includer.length)
			System.out.println("Inequal size");
		ArrayList<int[]> a = new ArrayList<int[]>();
		for(int i = 0 ; i < fixIt.size(); i++)
		{
			if(includer[i])
				a.add(fixIt.get(i));
		}
		return a;
	}
	public ArrayList<String> prunedMatrix(ArrayList<String> g)
	{
		if(g.size()!= inclusions.length)
			System.out.println("Inequal size");
		ArrayList<String> a = new ArrayList<String>();
		for(int i = 0 ; i < g.size(); i++)
		{
			if(inclusions[i])
				a.add(g.get(i));
		}
		System.out.println("gene" + a.size());
		return a;
	}
	public static void main(String args[])
	{
		SparseParse sp = new SparseParse(args[0], args[1]);
		ArrayList<int[]> read = sp.readMat();
		ArrayList<String> gigi = sp.readGenes();
		double[][] cors = sp.populateVectors(read);
		sp.writeFiles(cors, sp.prunedMatrix(gigi));
	}
}
