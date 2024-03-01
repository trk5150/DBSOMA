package correlationMatrixMaker;

import java.io.*;
import java.util.ArrayList;

public class SingleCellReader 
{
	public int count;
	public String file;
	public String lookup;
	
	public SingleCellReader(String fileName)
	{
		file = fileName;
		count = 0;
	}
	public double[][] parseTSV(double[][] d)  //[cells][genes]
	{
		BufferedReader br = null;
		int i = 0;
		String line = "";
		try 
		{
			br = new BufferedReader(new FileReader(file));
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	    try {
        br.readLine();//skip first line
        line = br.readLine();//start on second line
	        while (line != null)// && i<=0) 
	        {
	        	d[i] = parseLine(line.split("\\t"));
	        	i++;
	        	line = br.readLine();
//	        	if(i%100 == 0)
//	        		System.out.println(i);
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
		double[] parsed = new double[s.length-1];
		for(int i = 1; i<s.length; i++)
		{
			parsed[i-1] = Double.parseDouble(s[i]);
			//System.out.print(parsed[i-1]+ ", ");
		}
		return parsed;
	}
	
	
	public String[] geneNames()
	{
		BufferedReader br = null;
		String line = "";
			try 
			{
				br = new BufferedReader(new FileReader(file));
				
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		    try {
	        
		    line = br.readLine();
		    
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
	    String [] names = line.split("\\t");
		return names;
	}
	public int matSizer(String fileLoc)
	{
		BufferedReader br = null;
		int i = 0;
		String line = "";
		try 
		{
			br = new BufferedReader(new FileReader(file));
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	    try {
           // line = br.readLine();
	        while (line != null) //&& i<=0) 
	        {
	        	line = br.readLine();
//	        	if(i == 0)
//	        		System.out.println(line);
	        	i++;	
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
		return i;
		
	}
}
