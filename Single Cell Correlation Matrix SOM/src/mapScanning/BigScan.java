package mapScanning;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

/* To do: 
 * 	pass the radius and minPts.
 * 	change the structure cutoff based on quality and structure metrics of false positive first
 * 		The above will allow me to more quickly parameterize without needing to worry about input file sizes
 * 		calculate false positive rate for generic file sizes 200, 250, 500, 750, 1000... etc. Then for any file, whichever size is closest, use a cutoff of 1.5x that false positive rate
 * 
 * 	Change so that the state by which files are saved is not just first alphabetical
 * */

/**
 * This class takes 4 directories as arguments: 0: SOM files, 1: state files, 2: projection files, 3:output directory address
 * 
 * for each SOM file, it runs an instance of SetScan with all state files and projection files, which saves the output in an output sub-directory named after the SOM file
 */

public class BigScan 
{
	
	public String somDir, stateDir, projDir, outDir;
	public BigScan(String soms, String projections, String out, String states)
	{
		outDir = out;
		ArrayList<File> somFiles = prepDirectory(soms); 
		ArrayList<File> stateFiles = prepDirectory(states);
//		System.out.println("state number " + stateFiles.size());
		String[] stateAddresses =new String[stateFiles.size()];
		for(int i = 0; i< stateAddresses.length; i++)
		{
			stateAddresses[i] = stateFiles.get(i).getAbsolutePath();
		}
		
		for(int i = 0; i< somFiles.size(); i++)
		{
			String subName = somFiles.get(i).getName().substring(0,somFiles.get(i).getName().indexOf("."));

			System.out.println("SOM: " + subName + " being analyzed");
			
			File subOutDir = new File(outDir + "//" +subName);
			subOutDir.mkdir();
			
			String[] args = new String[5+stateAddresses.length];
			args[0] = "1";
			args[1] = "5";
			args[2] = somFiles.get(i).getAbsolutePath();
			args[3] = projections;
			args[4] = subOutDir.getAbsolutePath();
			for(int j =0; j< stateAddresses.length; j++)
			{
				args[5+j] = stateAddresses[j];
			}
			SOMAnalysis.main(args);
		}
	}
	public BigScan(String soms, String projections, String out, String states, String rad, String min)
	{
		outDir = out;
		ArrayList<File> somFiles = prepDirectory(soms); 
		ArrayList<File> stateFiles = prepDirectory(states);
//		System.out.println("state number " + stateFiles.size());
		String[] stateAddresses =new String[stateFiles.size()];
		for(int i = 0; i< stateAddresses.length; i++)
		{
			stateAddresses[i] = stateFiles.get(i).getAbsolutePath();
		}
		
		for(int i = 0; i< somFiles.size(); i++)
		{
			String subName = somFiles.get(i).getName().substring(0,somFiles.get(i).getName().indexOf("."));

			System.out.println("SOM: " + subName + " being analyzed");
			
			File subOutDir = new File(outDir + "//" +subName);
			subOutDir.mkdir();
			
			String[] args = new String[5+stateAddresses.length];
			args[0] = rad;
			args[1] = min;
			args[2] = somFiles.get(i).getAbsolutePath();
			args[3] = projections;
			args[4] = subOutDir.getAbsolutePath();
			for(int j =0; j< stateAddresses.length; j++)
			{
				args[5+j] = stateAddresses[j];
			}
			SOMAnalysis.main(args);
		}
	}
	public ArrayList<File> prepDirectory(String fold)
	{
		File folder = new File(fold);
		ArrayList<File> files = new ArrayList<File>();
		try {
			System.out.println("getting files from "+folder.getCanonicalPath());
		} catch (IOException e) {
			e.printStackTrace();
		}
		File[] listOfSoms = folder.listFiles();

		for (File file : listOfSoms) 
		{
		    if (file.isFile() && (file.getName().endsWith(".som") || file.getName().endsWith(".txt")|| file.getName().endsWith(".gmt"))) 
		    {
		    	files.add(file);
		    }
		}
		return files;
	}
	//need a special arg for the lead state?
	public static void main(String args[])
	{
		BigScan b;
		if(args.length>4)
			b = new BigScan(args[0], args[1], args[2], args[3], args[4], args[5]);
		else
			b = new BigScan(args[0], args[1], args[2], args[3]);
		System.out.println(b.toString());
	}
}
