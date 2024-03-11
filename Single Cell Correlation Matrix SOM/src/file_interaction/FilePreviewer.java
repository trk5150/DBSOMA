package file_interaction;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

public class FilePreviewer 
{
	public static void main(String[] args)
	{
		BufferedReader br = null;
		String line = "";
		
		try 
		{
			br = new BufferedReader(new FileReader(args[0]));
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	    try 
	    {
    	  /**checks for genes in a data set, gets read count*/
//	    	int i = 0;
//	    	int lineNum = 0;
//	    	
//	    	String[] genes = br.readLine().split("\t");
//	    	String check = "c1orf127";
//	    	int position = 0;
//	    	for(int cc = 0; cc < genes.length; cc++)
//	    	{
//	    		if(genes[cc].equalsIgnoreCase(check))
//    			{
//	    			position = cc;	
//	    			System.out.println(check + " at " + position);
//    			}
//	    	}
//	    	double countSum = 0;
//	    	while (line != null)// && i<=0) 
//	        {
//	    		line = br.readLine();
//	    		if(line !=null)
//	    		{
//		    		String[] lines = line.split("\\t");					/** change check to equal the lookup*/
//		    		if(lines.length>position)
//		    			countSum += Double.parseDouble(lines[position]);
//	    		}
//	        }
//	    	System.out.println(countSum);
//	    	System.out.println("total counts = " + countSum);
	    	
//	    	while (line != null)// && i<=0) 
//	        {
//	    		line = br.readLine();
//	    		if(line !=null)
//	    		{
//		    		String[] lines = line.split("\\t");
//		    		String check = "c1orf127";						/** change check to equal the lookup*/
//		    		if(lines[0].equalsIgnoreCase(check))
//		    		{
//		    			System.out.println(check + " at " + i);
//		    			System.out.println(line);
//		    		}
//	
//		    		i++;
//		    		lineNum++;
//	    		}
//	        }
	    	/** checks that the matrix is square*/
//	    	int finalLine = 0;
//	    	while (line != null)// && i<=0) 
//	        {
//	    		line = br.readLine();
//	    		if(line !=null)
//	    		{
//		    		if(i==0)
//		    		{
//		    			String[] lines = line.split("\\t");
//		    			System.out.println(lines[0]);
//		    		}
//		    		if(i==1)
//		    		{
//		    			String[] lines = line.split("\\t");
//		    			System.out.println(lines.length);
//		    			finalLine = lines.length;
//		    		}
//		    		if(lineNum>=finalLine)
//		    		{
//		    			//System.out.println(line);
//		    		}
//		    		i++;
//		    		lineNum++;
//	    		}
//	        }
//	    	System.out.println(lineNum-2); //lineNum -2 because first line is matrix information, last line is blank
/***	    	reads first x lines*/
	    	for(int i1 =0; i1<=1; i1++)
	    	{
	    		line = br.readLine();
	    		System.out.println(line);
	    		String[] lines = line.split("\\t");
	    		//System.out.println(lines.length);
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
	}
	
}
