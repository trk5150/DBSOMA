package correlationMatrixMaker;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

/** class contains the methods required to build a gene x gene expression correlation matrix from single cell RNA-seq data
 * 
 *  Takes 2 arguments: arg 1 = the RNA-seq file location; arg 2 = output folder
 *  
 *  the RNA-seq file must be:
 *  	-tab delimited
 *  	-array of integer read-counts or doubles of normalized counts with the genes in the columns
 *  	-gene names listed in the first row
 *  	-each row after the first the cells must start with the cell/barcode identifier followed by the counts for each gene
 *  
 *  Have other classes used for RNA-seq file wrangling to prep them for this input format
 *  
 *  **/


public class SingleCellParserTSV 
{
	public boolean trimLowCount;
	public boolean subCellPop;
	public double cellCutOffPercentage;
	public double countCutOffPercentAvg;
	private String fileLoc;
	public int cells;
	public double avgCount;
	public int geneNum;
	public double cutoffPortion;
	public String[] geneList;
	public double[] readCountTotals;
	public double[][] parsedCounts;
	public double[][]parsedCountTranspose;
	
	SingleCellParserTSV (String s)
	{
		trimLowCount = false;
		subCellPop = false;
		cellCutOffPercentage = 20; //must be above this % of average count for...
		cutoffPortion = 0.75;      //for this portion of the identity genes
		
		if(trimLowCount)
			countCutOffPercentAvg = 0.5;
		
		fileLoc = s;
		System.out.println(fileLoc);
		//int sizzle = matsSizer(fileLoc);
		//size = sizzle;
		cells = cellCounter(fileLoc);
			System.out.println("cells counted: " + cells);
		geneList =  geneNamer(fileLoc);
		geneNum = geneList.length;
			System.out.println("genes counted: " + geneNum);
		readCountTotals = new double[geneNum];
			System.out.println("array sizes calculated");
		parsedCounts = new double[cells][geneNum];
		parsedCounts = fillMat(fileLoc, parsedCounts);
			System.out.println("matrix filled");
		parsedCountTranspose = transposeMatrix(parsedCounts);
			System.out.println("matrix transposed");
		
		readCountTotals = totaler();
		
		if(subCellPop)
		{
			System.out.println("refining cells");
			cutCells();
			if(trimLowCount)
				readCountTotals = totaler(); //need new totals if we're going to be trimming after
		}
			
		if(trimLowCount)
		{
			trim();
			//readCountTotals = totaler();
		}
		
		
				
	}
	public void cutCells()
	{
		
		/**manual entry of identity genes should be moved to args!*/
		String[] IDMarks = new String[5]; IDMarks[4] = "GCK"; IDMarks[0] = "ins"; IDMarks[1] = "SLC2A1"; IDMarks[2] = "iapp"; IDMarks[3] = "PCSK2";		
		
		
		
		int[] idindexes = new int[IDMarks.length];
		for(int i = 0; i < idindexes.length; i++)
		{
			idindexes[i] = -1;
		}
		
		for(int i =0; i< IDMarks.length; i++)
		{
			for(int j = 0; j < geneList.length; j++)
			{
				if(IDMarks[i].equalsIgnoreCase(geneList[j]))
				{
					idindexes[i] = j;
				}
			}
			if(idindexes[i] ==-1)
			{
				System.out.println("one or more of your input cell identity markers, including at least: \"" + IDMarks[i]+ "\" is not observed in the data set. Refine you input set");
				System.out.println("no cells will be removed from the dataset");
				return;
			}
			
		}
		
		//testing code:
//		for(int i = 0; i < idindexes.length; i++)
//		{
//			double a = readCountTotals[idindexes[i]];
//			String b = geneList[idindexes[i]];
//			System.out.println(a + " " + b + " " + IDMarks[i] + " " + idindexes[i]);
//		}

		//store whether each cell has count above threshold for each given cell identifier gene
		int[][] boo = new int[cells][idindexes.length];
		double check = cellCutOffPercentage/100;
		//System.out.println(check);
		for(int i = 0; i < idindexes.length; i++)
		{
			int cGI = idindexes[i]; //currentGeneIndex
			
			for(int j = 0; j < parsedCountTranspose[cGI].length; j++) // parsedCountTranspose matrix is [genes][cells]
			{
				//System.out.println(parsedCountTranspose[cGI][j] + " < " + check*readCountTotals[cGI]/cells);
				if(parsedCountTranspose[cGI][j] > check*readCountTotals[cGI]/cells) // is the read count in this cell higher that x% of the average for this gene
					boo[j][i] = 1;
				else
					boo[j][i] = 0;
				
			}
		}
		//decide whether to include or exclude each cell based on number of genes above threshold it has
		boolean[] include = new boolean[cells];
		int inclusionThreshold = (int)(idindexes.length*cutoffPortion);
		//System.out.println(inclusionThreshold);
		int includedNo = 0;
		for(int i = 0 ; i < boo.length; i++)
		{
			int toto = 0;
			for(int j = 0; j < boo[i].length; j++)
			{
				toto += boo[i][j];
			}
			if(toto >= inclusionThreshold)
			{
				include[i] = true;
				includedNo ++;
			}
			else
				include[i] = false;
		}
		System.out.println("including " + includedNo + " cells");
		double[][] refined = new double[includedNo][parsedCounts[0].length]; //[cells][genes]
		int refinedI = 0;
		for(int i = 0; i < parsedCounts.length; i++)  //for each cell
		{
			if(include[i])
			{
				refined[refinedI] = parsedCounts[i];
				refinedI++;
			}
		}
		cells = includedNo;
		parsedCounts = refined;
		parsedCountTranspose = transposeMatrix(refined);
//		System.out.println("reached exit command");
//		System.exit(0);
	}
	
	public double[] totaler()
	{
		double[] tots = new double[geneNum];
		int sum;
		double runningSum = 0;
		for(int i= 0; i<geneNum; i++) //for each gene
		{
			sum = 0;
			for(int j = 0; j < cells; j++) //go through all cells for each gene
			{
				sum += parsedCountTranspose[i][j]; //matrix is [genes][cells]
			}
			tots[i] = sum;
			runningSum+=sum;
			//System.out.println(sum);
		}
		//System.out.println(parsedCountTranspose.length + " "+parsedCountTranspose[0].length);
		avgCount = runningSum/geneNum;
		System.out.println("avg count is "+ avgCount);
		return tots;
	}
	public void trim()
	{
		double trimCut = avgCount*(countCutOffPercentAvg/100);
		System.out.println("cutoff is: "+trimCut);
		
		//get gene x cell matrix index for genes below cutoff
		ArrayList<Integer> lowIndexes = new ArrayList<Integer>();
		ArrayList<Integer> highIndexes = new ArrayList<Integer>();
		int lowGeneCount = 0;
		int highGeneCount = 0;
		for(int i = 0; i< readCountTotals.length; i++)
		{
			if(readCountTotals[i] < trimCut)
			{
				lowIndexes.add(i);
				lowGeneCount++;
			}
			else
			{
				highIndexes.add(i);
				highGeneCount++;
			}
				
		}
		System.out.println("genes below cutoff: "+ lowGeneCount + "; genes above = " + highGeneCount);
		
		//remove genes below cutoff from gene x cell matrix
		String[] geneSubList = new String[highGeneCount];
		double[][] subMatrix = new double[highGeneCount][cells];
		
		int k = 0;
		for(int i = 0; i< highIndexes.size(); i++)
		{
			geneSubList[k] = geneList[highIndexes.get(i)];
			subMatrix[k] = parsedCountTranspose[highIndexes.get(i)];
			k++;
		}
		geneList = geneSubList;
		parsedCountTranspose = subMatrix;
		geneNum = highGeneCount;
		
	}
	/**unfinished method aimed at normalizing read counts before calculating correlations*/
//	public double[][] normalizer(double[][] raw, double[] tots)
//	{
//		double[][] norms = new double[raw.length][raw[0].length];
//		//for(int i = 0; i)
//		
//		
//		return norms;
//	}
	public int cellCounter(String loc)
	{
		SingleCellReader r = new SingleCellReader(loc);
		return r.matSizer(loc);
	}
	public double[][] transposeMatrix(double[][] matrix)
	{
		System.out.println("transposing");
	    int m = matrix.length;
	    int n = matrix[0].length;

	    double[][] transposedMatrix = new double[n][m];

	    //System.out.println("m = " + m + "; n = " + n);
	    //System.out.println(transposedMatrix.length + "\t" + transposedMatrix[0].length);
	    
	    for (int i = 0; i < matrix.length; i++)
	    {
	    	for(int j = 0; j <matrix[i].length; j++)
	    	{
	    		double k = matrix[i][j];
	    		transposedMatrix[j][i] = k;
//	    		transposedMatrix[j][i] = matrix[i][j];
	    	}
	    }
	    return transposedMatrix;
	}

	public String[] geneNamer(String loc)
	{
		SingleCellReader r = new SingleCellReader(loc);
		String[] sizzle = r.geneNames();
		int n = sizzle.length;
		String[] sizer = new String[n-1];
		for(int i = 1; i<n; i++)
		{
			sizer[i-1] = sizzle[i];
			//System.out.print(sizer[i-1]+", "); //test the gene list is correct
//			if(sizer[i-1].equalsIgnoreCase("c1orf127")||sizer[i-1].equalsIgnoreCase("ins")||sizer[i-1].equalsIgnoreCase("iapp"))
//			{
//				System.out.println(sizer[i-1]+" "+i);
//			}
		}
		return sizer;
	}
	public double[][] fillMat(String loc, double[][] d)
	{
		SingleCellReader r = new SingleCellReader(loc);
		double[][] read = r.parseTSV(d);
//		for(int i = 0; i<read.length; i++)
//		{
//			if(read[i].length != read[0].length)
//			{
//				System.out.println(i);
//			}
//		}
		return read;
	}
	public double[] findGeneIndicies(String[] queries)
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
		double[] indy = new double[indexer.size()];
		for(int i = 0; i < indexer.size(); i++)
		{
			indy[i] = indexer.get(i);
		}
		return indy;
	}
	
	public static void main(String [] s)
	{
		File f = new File(s[0]);
		String fileLoc = f.getParent();
		System.out.println(fileLoc);
		SingleCellParserTSV parse = new SingleCellParserTSV(s[0]);
		//System.out.println(parse.pairwiseCorrelation.length + " cells, " + parse.pairwiseCorrelation[0].length+ " genes");
		
		System.out.println("Correlation start");
		Correlator cors = new Correlator();
		double[][] correlationMatrix = cors.buildCorrelationMatrix(parse.parsedCountTranspose);
		System.out.println("finished Correlation");
		System.out.println("printing correltation matrix and gene list files");	
		FilesWriter w = new FilesWriter();
		w.writeFiles(correlationMatrix,parse.geneList,parse.fileLoc,parse.cells,parse.geneNum,fileLoc, parse.trimLowCount, parse.subCellPop);
		
//		String[] ofInterest = new String[7]; 
//		ofInterest[0] = "c1orf127";
//		ofInterest[1] = "ins";
//		ofInterest[2] = "iapp";
//		ofInterest[3] = "gcg";
//		ofInterest[4] = "ppy";
//		ofInterest[5] = "sst";
//		ofInterest[6] = "ubc";
		
		//System.out.println("printing top hits");
		//w.writeTopHits(parse.findGeneIndicies(ofInterest), correlationMatrix, parse.geneList, s[1]);
		
		System.out.println("\nargot");
	}
}
