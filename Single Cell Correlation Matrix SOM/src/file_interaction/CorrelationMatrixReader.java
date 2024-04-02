package file_interaction;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

public class CorrelationMatrixReader 
{
	public boolean skipLineOne = true;
	public int count;
	public int geneNum;
	public String matFile;
	public String geneNamesFile;
	public String lookup;
	public String[] genes;
	public double[][] matrix; 
	
	public CorrelationMatrixReader(String mat, String geneFile)
	{
		matFile = mat;
		geneNamesFile = geneFile;
		count = 0;
		genes = geneNames(geneFile);
		geneNum = genes.length;
		matrix = parseCorMatFile(mat,geneNum);
		
	}
	public double[][] parseCorMatFile(String file, int sizer)
	{
		BufferedReader br = null;
		int i = 0;
		double[][] d = new double[sizer][sizer];
		String line = "";
		
		try 
		{
			br = new BufferedReader(new FileReader(file));
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	    try {
	    	if(skipLineOne)
	    		br.readLine(); //skip line 1
	    	line = br.readLine();//start on line 2
	    	while (line != null)// && i<=0) 
	        {
	    		if(line.isBlank())
	        		break;
	        	d[i] = parseLine(line.split("\\t"));
	        	i++;
	        	if(i%1000 == 0)
	        		System.out.println(i);
	        	if(line.isBlank())
	        		break;
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
		return d;
	}
	
	
	public double[] parseLine(String[] s)
	{
		double[] parsed = new double[s.length];
		//if(s.length>0 && !s[0].isEmpty())
		{
			for(int i = 0; i<s.length; i++)
			{
				parsed[i] = Double.parseDouble(s[i]);
				//System.out.print(parsed[i-1]+ ", ");
			}
		}
		return parsed;
	}
	
	
	public String[] geneNames(String file)
	{
		BufferedReader br = null;
		String names = "";
		String line = "";
		try 
		{
			br = new BufferedReader(new FileReader(file));
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	    try {
	    	line = br.readLine();
	    	while(line != null)
    		{
	    		names = names + line;
	    		line = br.readLine();
	    		names = names+"\t";
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
		System.out.println(names);
	    String [] namer = names.split("\\t");
		return namer;
	}
	public double[][] getMatrix()
	{
		return matrix;
	}
	public String[] getGenes()
	{
		return genes;
	}
	
}
