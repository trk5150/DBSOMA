package mapScanning;

import java.awt.Polygon;
import java.io.File;
import java.util.ArrayList;

public class PaintSetScan 
{
	public int height = 200;
	public int width = 200;
	public DrawSet d, dBlue;
	ArrayList<int[][]> nodeScans;
	ArrayList<Polygon> nodeHexs;
	ArrayList<File> files;
	double[][] polyCenterCoords;
	private File outDir;
	private File trimmedOutDir;
	private File imageOutDir;
	public int[] weights;
	
	public PaintSetScan(ArrayList<int[][]> scans, int min, int xN, int yN, int[] nodeWeights, File outie, ArrayList<File> fileNames)
	{
		weights = nodeWeights;
		nodeScans = scans;
		files = fileNames;
		d = new DrawSet(min, xN, yN, nodeWeights);
		dBlue = new DrawSet(min, xN, yN, nodeWeights);
		outDir = outie;
		trimmedOutDir = new File(outDir + "\\Trimmed");
		imageOutDir = new File(outDir + "\\Images");
		
	}
	public void searchAll(ArrayList<int[][]> scans, boolean stateFiles)
	{
		outDir.mkdir();
		
		imageOutDir.mkdir();
		trimmedOutDir.mkdir();
		
		
		d.nodeBuild(2000,2000);
    	//d.colors();
		if(stateFiles)
		{
			int[][] wholeMap = new int[weights.length][2];
			for(int i = 0; i < weights.length; i++)
				wholeMap[i][0] = weights[i];
			
			d.projectSet(wholeMap, false);
        	d.saveImg("whole map", outDir,  2000, 2000);
			
			for(int i = 0; i < scans.size(); i++)
	    	{
	        	d.projectSet(scans.get(i), true);
	        	d.saveImg(files.get(i).getName(), outDir,  2000, 2000);
	    	}
		}
		else
		{
	    	for(int i = 0; i < scans.size(); i++)
	    	{
	    		
	    		d.projectSet(scans.get(i), false);;
	        	d.saveImg(files.get(i).getName(), imageOutDir,  2000, 2000);
	        	
	        	d.projectSet(scans.get(i), true);
	        	d.saveImg(files.get(i).getName(), trimmedOutDir,  2000, 2000);
	    	}
		}
	}
	
}
