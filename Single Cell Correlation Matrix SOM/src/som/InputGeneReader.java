package som;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

public class InputGeneReader 
{
	public String fileLoc;
	public String geneQC;
	public boolean debug;
	
	public InputGeneReader(String s, boolean deb)
	{
		debug = deb;
		geneQC = "file not yet read";
		fileLoc = s;
	}
	
	public String[] readGenes(int geneNum)
	{
		//file contains a line split list of genes with no header information, so no need to skip any initial lines
		String[] s = new String[geneNum];
		
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
	    	geneQC = "file read";
	    
	    	for(int i =0; i<geneNum; i++)
	    	{
	    		s[i] = br.readLine();
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
		    for(int p = 0; p<s.length; p++)
		    {
		    	System.out.println(s[p]);
		    }
		    System.out.println(geneQC);
	    }
		return s;
	}
	
}
