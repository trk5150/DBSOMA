package CorMatAnalyzer;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Arrays;

import correlationMatrixMaker.FilesWriter;
import correlationMatrixMaker.SingleCellParserTSV;
import file_interaction.CorrelationMatrixReader;

public class GeneLookUp 
{
	public String[] geneList;
	
	public GeneLookUp(String[] genes)
	{
		geneList = genes;
	}
	public int[] findGeneIndicies(String[] queries)
	{
		ArrayList<Integer> indexer = new ArrayList<Integer>();
		
		for(int i = 0; i< queries.length; i++)
		{
			int k = -1;
			for(int j = 0; j <geneList.length; j++)
			{
				if(queries[i].equalsIgnoreCase(geneList[j]))
				{
					k=j;
					System.out.println(geneList[j] + " " + j);
				}
				
			}
			if(k>-1)
				indexer.add(k);
		}
		int[] indy = new int[indexer.size()];
		for(int i = 0; i < indexer.size(); i++)
		{
			indy[i] = indexer.get(i);
		}
		return indy;
	}
	
	public static void main(String[] args)
	{
		CorrelationMatrixReader r = new CorrelationMatrixReader(args[0], args[1]);
		String[] genie = r.getGenes();
		double[][] d= r.getMatrix();
		String[] ofInterest = r.geneNames(args[2]);
		GeneLookUp g = new GeneLookUp(genie);
		int[] interestIndexer = g.findGeneIndicies(ofInterest);
		FilesWriter w = new FilesWriter();
		
		w.writeTopHits(interestIndexer, d, genie, args[3]);
	}
}
