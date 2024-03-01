package som;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

public class InputMatrixReader 
{
	/** 
	 * reads a matrix file, from given pathway, and generates a 2D array of doubles. Matrix file must be symmetric after a first line which is QC information.
	 * */
	
	public String fileLoc;
	public String matrixQC;
	public boolean debug;
	
	public InputMatrixReader(String a, boolean deb) 
	{
		debug = deb;
		matrixQC = "file not yet read";
		fileLoc = a;
	}
	public double[][] readMat()
	{

		double[][] mat = new double[0][0];
		BufferedReader br = null;
		String line = "";
		try 
		{
			br = new BufferedReader(new FileReader(fileLoc));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	    try 
	    {
	    	//first line is QC
	    	matrixQC = br.readLine();
	    	
	    	//Second line is the initial reading of data, needs to be handled differently
	    	line = br.readLine();
    		String[] lines = line.split("\\t");
    		
    		if(debug) {System.out.println("genes num = "+lines.length);}
	    
    		int geneNum = lines.length;
    		mat = new double[geneNum][geneNum];
    		
	    	//fill the rest of the array
	    	for(int i =0; i<geneNum; i++)
	    	{
	    		for(int j = 0; j<geneNum; j++)
	    		{
	    			double ddd = (double)Double.parseDouble(lines[j]);
	    			if(Double.isNaN(ddd))
	    			{
	    				System.err.println("NaN detected at line i=" + i + ";  position j =" + j + "\t set to zero by default to allow clustering to attempt to continue. Recheck input matrix to ensure NaNs are not an issue");
	    				ddd= 0;
	    			}
	    			mat[i][j] = ddd;
	    		}
	    		if(i<geneNum-1)
	    		{
	    			line = br.readLine();
	    			lines = line.split("\\t");
	    		}
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
	    
	    //reads out the matrix, for testing purposes
	    if(debug)
	    {
		    for(int p = 0; p<mat.length; p++)
		    {
		    	System.out.println(mat[p][p]);
//		    	for(int j = 0; j< mat.length; j++)
//		    	{
//					System.out.print(mat[p][j]+"\t");
//		    	}
		    	System.out.print("\n");
		    }
		    System.out.println(matrixQC);
		    System.out.println(mat.length);
	    }
	    
		return mat;
	}
}
