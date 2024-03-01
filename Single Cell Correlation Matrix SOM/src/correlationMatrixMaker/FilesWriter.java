package correlationMatrixMaker;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collections;

public class FilesWriter 
{
	/** should move the top hits calculation machinery to it's own class */
	public void writeTopHits(int[] queries, double[][] corMatrix, String[] genes, String outPrefix)
	{
		//System.out.println(queries[0]);
		BufferedWriter b;
		for(int i = 0; i< queries.length; i++)
		{
			if(queries[i] == -1) //gene of interest not found in gene list
				i++;
			if (i >= queries.length)
				break;
			
			double[] corsToRank = corMatrix[queries[i]];
			
			int[] topHits = topHitRanker(corsToRank);
			
			String filePrefix = outPrefix + "\\";
			String topHitsFileName = filePrefix+ genes[queries[i]]+".txt";
			try 
			{
				System.out.println("Printing top hits to "+topHitsFileName);
							
				FileWriter ff = new FileWriter(topHitsFileName,true);
				b = new BufferedWriter(ff);
				PrintWriter printer = new PrintWriter(b);
				printer.print(genes[queries[i]]+"\n");
				for(int j =0; j <topHits.length; j++)
				{
					//System.out.println(topHits[j]);
					printer.print(genes[topHits[j]]+"\t" + corMatrix[queries[i]][topHits[j]] + "\n");
				}
				
				printer.close();
			}catch (IOException e){
				e.printStackTrace();
				System.err.print("Error printing file "+topHitsFileName);
			}
		}
	}
	public int[] topHitRanker(double[] cors)
	{
		int[] sortedIndicies = new int[cors.length];
		double[] sorted = new double[cors.length];
		for(int i = 0; i<cors.length; i++)
		{
			sorted[i] = cors[i];
		}
		Arrays.sort(sorted);
		double[] reversedSort = new double[sorted.length];
		int n = reversedSort.length;
		for(int i = 0; i < n; i++)
		{
			reversedSort[n-i-1] = sorted[i];
		}
		double[] checker = new double[cors.length];
		for(int i = 0; i<checker.length; i++)
		{
			checker[i] = cors[i];
		}
		for(int i = 0; i< reversedSort.length; i++)
		{
			//System.out.println(reversedSort[i]);
			for (int j = 0; j< reversedSort.length; j++)
			{
				if(reversedSort[i] == checker[j])
				{
					checker[j] = -1000; 
					//System.out.println(reversedSort[i]);
					//System.out.println("sorted[i] = " + sorted[i] +" and cors[j] = "+ cors[j] + " so " +i + " set to: " + j);
					sortedIndicies[i] = j;
					break;
				}
			}
		}
		
		return sortedIndicies;
	}
	public void writeFiles(double[][] corMatrix, String[] genes, String sourceFile, int cells, int geneNum, String outPrefix, boolean trimmed, boolean subCell)
	{		
		BufferedWriter b;
		String filePrefix;
		if(trimmed && subCell)
			filePrefix = outPrefix + "\\Training_set_TRIMMED_SUBCELL";
		filePrefix = outPrefix + "\\Training_set";
		if(trimmed)
			filePrefix+="_TRIMMED";
		if(subCell)
			filePrefix+="_SubCell";
		String corFilename = filePrefix+"_Correlations.txt";
		String genesFilename = filePrefix+"_Genes.txt";
		
		//Correlations file
		try 
		{
			System.out.println("Printing correlations to "+corFilename);
						
			FileWriter ff = new FileWriter(corFilename,true);
			b = new BufferedWriter(ff);
			PrintWriter printer = new PrintWriter(b);
			
			printer.print("trained from: "+sourceFile+ ", with "+geneNum + " genes and "+ cells +" cells \n");
			for(int i = 0; i<corMatrix.length; i++)
			{
				for(int j = 0; j<corMatrix[i].length; j++)
				{
					if(Double.isNaN(corMatrix[i][j]))
					{
						System.err.println("NaN detected at line i=" + i + ";  position j =" + j + "\t set to zero by default to allow clustering to attempt to continue. Recheck input matrix to ensure NaNs are not an issue");
	    				corMatrix[i][j]= 0;
					}
					printer.print(corMatrix[i][j]+"\t");
				}
				printer.print("\n");
			}
			printer.print("\n");
			printer.close();
		}catch (IOException e){
			e.printStackTrace();
			System.err.print("Error printing SOM file "+corFilename);
		}
		
		//Genes file
		try 
		{
			System.out.println("Printing genes to "+genesFilename);
						
			FileWriter ff = new FileWriter(genesFilename,true);
			b = new BufferedWriter(ff);
			PrintWriter printer = new PrintWriter(b);
			for(int i = 0; i<genes.length; i++)
				printer.print(genes[i]+"\n");
			printer.close();
		
		}catch (IOException e){
			e.printStackTrace();
			System.err.print("Error printing info file "+genesFilename);
		}
		
	}
}
