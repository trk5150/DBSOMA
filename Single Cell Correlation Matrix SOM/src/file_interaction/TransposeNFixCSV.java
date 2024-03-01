package file_interaction;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class TransposeNFixCSV 
{
	public int[][] transposeMatrix(int[][] matrix)
	{
		System.out.println("transposing");
	    int m = matrix.length;
	    int n = matrix[0].length;

	    int[][] transposedMatrix = new int[n][m];

	    System.out.println("m = " + m + "; n = " + n);
	    System.out.println(transposedMatrix.length + "\t" + transposedMatrix[0].length);
	    
	    for (int i = 0; i < matrix.length; i++)
	    {
	    	for(int j = 0; j <matrix[i].length; j++)
	    	{
	    		int k = matrix[i][j];
	    		transposedMatrix[j][i] = k;
//	    		transposedMatrix[j][i] = matrix[i][j];
	    	}
	    }
	    return transposedMatrix;
	}
	
	public int matSizer(String fileLoc)
	{
		BufferedReader br = null;
		int i = 0;
		String line = "";
		try 
		{
			br = new BufferedReader(new FileReader(fileLoc));
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	    try {
           // line = br.readLine();
	        while (line != null) //&& i<=0) 
	        {
	        	if(line.length()>10)
	        		i++;	
	        	line = br.readLine();
//	        	if(i == 0)
//	        		System.out.println(line);
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
	public String[] parseCSV(int[][] d, String[] g, String file)  //returns a comma split first line string array
	{
		BufferedReader br = null;
		int i = 0;
		String line = "";
		String[] cells = new String[0];
		try 
		{
			br = new BufferedReader(new FileReader(file));
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	    try {
        String l1 = br.readLine();//special first line
        cells = l1.split(",");
        line = br.readLine();//start on second line
	        while (line != null)// && i<=0) 
	        {
	        	String[] s = line.split(",");
	        	if(s[0].contains("_"))
	        		g[i] = s[0].substring(1, s[0].indexOf("_"));
	        	else
	        		g[i] = s[0];
	        	d[i] = parseLine(s);
	        	i++;
	        	line = br.readLine();
//	        	if(i%1000 == 0)
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
	    return cells;
	}
	
	
	public int[] parseLine(String[] s)
	{
		int[] parsed = new int[s.length-1];
		for(int i = 1; i<s.length; i++)
		{
			parsed[i-1] = (int)Double.parseDouble(s[i]);
			//System.out.print(parsed[i-1]+ ", ");
		}
		return parsed;
	}
	public static void main(String[] args)
	{
		File f = new File(args[0]);
		String fileRoot = f.getName().substring(0,f.getName().indexOf("."));
		//System.out.println(fileRoot + " from \n" +f.getName() + " from \n" + f.getPath() + " from \n" + args[0]);
		TransposeNFixCSV transpo = new TransposeNFixCSV();
		int lines = transpo.matSizer(args[0]);
		String[] genes = new String[lines-1];
		int[][] data = new int[lines-1][0];
		String[] cells = transpo.parseCSV(data,genes,args[0]);
		
		
		int [][] flippedData = transpo.transposeMatrix(data);
		
		transpo.writeFiles(flippedData, genes, cells, args[1], fileRoot);
		
	}

	public void writeFiles(int[][] data, String[] genes, String[] cells, String outPrefix, String originalFile)
	{		
		BufferedWriter b;
		String path = outPrefix+"\\transposed_from_"+originalFile+".txt";
		try 
		{
			System.out.println("Printing correlations to "+path);
						
			FileWriter ff = new FileWriter(path,true);
			b = new BufferedWriter(ff);
			PrintWriter printer = new PrintWriter(b);
			for(int k = 0; k <genes.length; k++)
			{
				printer.print("\t" + genes[k]);
			}
			printer.print("\n");
			for(int i = 0; i<data.length; i++)
			{
				printer.print(cells[i+1]+"\t");
				for(int j = 0; j<data[i].length; j++)
				{
					printer.print(data[i][j]+"\t");
				}
				printer.print("\n");
			}
			printer.print("\n");
			printer.close();
		}catch (IOException e){
			e.printStackTrace();
			System.err.print("Error printing SOM file "+path);
		}	
	}
}
