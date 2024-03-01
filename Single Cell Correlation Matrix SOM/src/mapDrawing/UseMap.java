package mapDrawing;
import java.awt.Color;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Scanner;

import som_analysis.SingleReader;



public class UseMap
{
	public DrawHex d;
	public String fold, mapper;
	public String[] geneStates;
	public File outDir, trimmedOutDir;
	public double gini, pVal;
	public int nodes,xs,ys,xNodes,yNodes, maxBins, minBins, colorNum, winW, winH, pValnVal, binSize;
	public ArrayList<DataPoint> bins;
	public ArrayList<MiniNode> nodeList;
	public ArrayList<File> searchers;
	public MiniSystem nodeSystem;
	public boolean addToOldFile, pics, equalWeight, trimmedPic;
	public boolean printGiniCurves=false;
	public double[][] g;
	public double[] sep;
	public double giniWLArea =0;
	public Color heatmapColorA=Color.red, heatmapColorB=Color.blue;
	
	public UseMap(String map, String fillle, String outdir, boolean pic, boolean equalWeighting, String[] states)
	{
		geneStates = states;
		trimmedPic = true;
		pics = pic;
		mapper = map;
		fold = fillle;
		File folder = new File(fillle);
		outDir = new File(outdir);
		trimmedOutDir = new File(outdir + "\\Trimmed");
		pValnVal= 1000;
		gini=0;
		String mapp = map;
		equalWeight = equalWeighting;
		yNodes=0;
		xNodes=0;
		searchers = new ArrayList<File>();
		mapReader(mapp);
		pVal = 0;
		binSize = bins.get(0).maxLocus - bins.get(0).minLocus;
		try {
			System.out.println("Analyzing files in "+folder.getCanonicalPath());
		} catch (IOException e) {
			e.printStackTrace();
		}
		File[] listOfFiles = folder.listFiles();

		for (File file : listOfFiles) 
		{
		    if (file.isFile() && file.getName().endsWith(".grp") || file.getName().endsWith(".txt")) 
		    {
		    	searchers.add(file);
		    }
		}
		nodes = yNodes*xNodes;
		
		heatMapping();
	}
	public void setColorA(Color c) {heatmapColorA=c;}
	public void setColorB(Color c) {heatmapColorB=c;}
	
	public void searchSystemDouble()
	{
		if(pics)
		{
			outDir.mkdir();
			d = new DrawHex(mapper);
			d.setColorA(heatmapColorA);
			d.setColorB(heatmapColorB);
			d.setEqualWeight(equalWeight);
			d.nodeBuild(2000,2000);
	    	d.heatMapping();
	    	d.countingDPS(null);
	    	d.saveImg("whole_map", outDir,  2000, 2000);
		}
		ArrayList<String> writer = new ArrayList<String>();
		writer.add("File\tgini\tginiW/WeightedBaseline(n="+pValnVal+")\tzScore");
		for (int i = 0; i < searchers.size()-1; i++)
		{
			for(int j = i+1; j < searchers.size(); j++)
			{
				g = new double[nodes+2][pValnVal+1];
				for(int x=0; x<nodes+2; x++) { for(int y=0; y<=pValnVal; y++) {g[x][y]=0;} } 
				gini = 0;
				sep = new double[pValnVal+1];
				multiSearch(searchers.get(i), searchers.get(j));
				writer.add(searchers.get(i) + "_vs_" + searchers.get(j) + coClustering(searchers.get(i), searchers.get(j)));
			}
		}
		writeFile(writer);
		System.out.println("done");
	}
	public void searchSystem()
	{
		if(pics)
		{
			outDir.mkdir();
			d = new DrawHex(mapper);
			d.setColorA(heatmapColorA);
			d.setColorB(heatmapColorB);
			d.setEqualWeight(equalWeight);
			d.nodeBuild(2000,2000);
	    	//d.colors();
	    	d.heatMapping();
	    	d.countingDPS(null);
	    	d.saveImg("whole_map", outDir,  2000, 2000);
		}
		ArrayList<String> writer = new ArrayList<String>();
		String stately = "";
		for(int i = 0; i < geneStates.length; i++)
		{
			stately = stately.concat("\tpercent overlap with " + geneStates[i]);
		}
		//System.out.println(stately);
		writer.add("File\tgini\tginiW/WeightedBaseline(n="+pValnVal+")\tzScoregini\tSearchCount\tFoundCount" + "\t" + "TrimmedCount" + "\t" + "TrimmedCount/FoundCount" + stately);
		SingleReader r = new SingleReader(false);
		r.readSOM(mapper);
		r.mapMaker();
		
		
		ArrayList<ArrayList<String>> stateLists = new ArrayList<ArrayList<String>>();
		for(String fifi : geneStates)
		{
			//System.out.println(s.getName());
			g = new double[nodes+2][pValnVal+1]; 
			for(int x=0; x<nodes+2; x++) { for(int y=0; y<=pValnVal; y++) {g[x][y]=0;} }
			gini = 0;
			sep = new double[pValnVal+1];
			//above code seems like it shouldn't be necessary, but causes errors when not included
			
			File fifo = new File(fifi);
			System.out.println(fifo.getName());
			stateLists.add(stateSearch(fifo, r));			
		}
		int fileCount = 0;
		for (File s:searchers)
		{
			//System.out.println(s.getName());
			g = new double[nodes+2][pValnVal+1]; 
			for(int x=0; x<nodes+2; x++) { for(int y=0; y<=pValnVal; y++) {g[x][y]=0;} }
			gini = 0;
			sep = new double[pValnVal+1];
			
			
			double[] nf = search(s, r, stateLists);
			
			String report = "";
			for(int i = 0; i < nf.length; i++)
			{
				//System.out.println(nf[i]);
				report = report +("\t" + nf[i]);
			}
			//System.out.println(report);
			writer.add(s.getName() + "\t" + g[nodes+1][0] +"\t" + g[nodes][0] + "\t" + g[nodes+1][1] + report);
			if(fileCount%10000==0)
				System.out.println(fileCount);
			fileCount++;
		}
		writeFile(writer);
		System.out.println("done");
	}
	public void multiSearch(File file1, File file2)
	{
		if(pics)
		{
			System.out.println(file1.getName()+" vs. "+file2.getName());
			d.multiSearch(file1, file2);
			d.repaint();
			d.saveImg(file1.getName()+"_vs_"+file2.getName(), outDir, 2000, 2000);
		}
	}
	public double[] search(File file, SingleReader r, ArrayList<ArrayList<String>> stateLists)
	{
		int shell = 2;
		int requiredCalls = 5;
		
		for(int i = 0; i < nodeSystem.size(); i++)
		{
			MiniNode mini = nodeSystem.get(i);
			mini.counting.clear();
			mini.weight = 0;
		}
		
		BufferedReader br = null;
		ArrayList<Integer> inds = new ArrayList<Integer>();
		String line = "";
		int foundCount=0;
    	int tryCount=0;
		try 
		{
			br = new BufferedReader(new FileReader(file));
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		ArrayList<String> foundGenes = new ArrayList<String>();
	    try {
            line = br.readLine();

            while (line != null) 
	        {
	        	String[] s = new String[1];
	        	s[0] = line;
	        	if(line.contains(","))
	        	{
	        		s = line.split(",");
	        	}
//	        	if(s.length>0)
//	        		System.out.println(s[0] + " " + s.length);
	        	
	        	for(int k = 0; k<s.length; k++)
	        	{
	        		
	        		if(s[k].length()>1)
	        		{
		        		r.geneCheck(s[k], shell);//for trim analysis
	        			tryCount++;
		        		boolean found = false;
			        	for(int i = 0; i< d.dataPoints.size(); i++)
			        	{
			        		if(d.dataPoints.get(i).name.equalsIgnoreCase(s[k].trim()))
			        		{
			        			found = true;
			        			foundCount++;
			        			d.dataPoints.get(i).myMini.counting.add(d.dataPoints.get(i));
			        			inds.add(i);
			        			foundGenes.add(s[k]);
			        		}
			        		
			        	}
		        	}
	        	}
	        	
	        	
	        	line = br.readLine();
        	}
//            System.out.println("searched for " + tryCount + " genes, found "+foundCount);   
		    } catch (IOException e)
		    {
				e.printStackTrace();
			}
	    	int numberOfStates = stateLists.size();
	    	
    		//count the number of genes that were found on the SOM by the current search file that are also included in each state
	    	double[] percentStates = new double[numberOfStates];
	    	for(int i = 0; i< percentStates.length; i++)
	    	{
	    		percentStates[i] = 0; 
	    		ArrayList<String> currentState = stateLists.get(i);
	    		double foundAndInState = 0;
	    		if(currentState.size()>0 && foundGenes.size()>0)
	    		{
	    			for(String geneCheck: foundGenes)
	    			{
	    				if(currentState.contains(geneCheck))
	    				{
	    					foundAndInState++;
	    				}
	    			}

		    		percentStates[i] = foundAndInState / (double)foundGenes.size();
	    		}
	    	}
	 
	  		ArrayList<String> trimmedList = r.countMultipleRepeatedGenes(requiredCalls, false); //for trim analysis
			
			//Random assignment of the same counts and weights to nodes
			for(int p = 1; p<g[0].length; p++){
				for(int z=0; z<foundCount; z++){
					int ree = (int) (Math.random()*bins.size());
					ree = nodeList.indexOf(bins.get(ree).myMini);
					g[ree][p] ++;
				}
			}
		
			//Actual assignment of the weight to all (count) assigned nodes
			for(Integer i: inds)
			{
				bins.get(i).myMini.counting.add(bins.get(i));
				bins.get(i).myMini.weight ++;
				g[nodeList.indexOf(bins.get(i).myMini)][0] ++;

			}
			
			giniWithWeightedBaseline();
			gini();
			pVal();
			double[] genio = new double[4];
			int trimmedCount = trimmedList.size();
			genio[0] = tryCount;
			genio[1] = foundCount;
			genio[2] = trimmedCount;
			genio[3] = trimmedCount / (double)foundCount*Math.log(foundCount);
			double zzz = Math.abs(g[nodes+1][1]);
			//System.out.println(zzz);

			double percentThreshold = -1;
			boolean states = true;
			if(states)
				percentThreshold = 16;
			
			
			//if(pics && foundCount > 15 && foundCount < 3500 && ((zzz >= 3 && genio[3] > 3) || (zzz>=5) || genio[3]>7)) //ad hoc cluster calling metrics (rare enough gini and tuned ratio or exceptional individual value)

			if(pics && percentStates[0]*100 > percentThreshold && foundCount > 25 && foundCount < 3500 && ((zzz >= 3 && genio[3] > 3) || (zzz>=5) || genio[3]>7)) //ad hoc cluster calling metrics (rare enough gini and tuned ratio or exceptional individual value)
			{
				System.out.println(file.getName());
				d.equalWeight = equalWeight;
				d.search(file);
				d.repaint();
				d.saveImg(file.getName(), outDir, 2000, 2000);
				if(trimmedPic)
				{
					trimmedOutDir.mkdir();
					for(int i = 0; i < nodeSystem.size(); i++)
					{
						MiniNode mini = nodeSystem.get(i);
						mini.counting.clear();
						mini.weight = 0;
					}
					
					//flip colors for the trimmed maps for easier viewing
					d.setColorA(heatmapColorB);
					d.setColorB(heatmapColorA);
					
					//search the Trimmed List of genes
					d.search(trimmedList);
					d.repaint();
					d.saveImg(file.getName() + "_trimmed", trimmedOutDir, 2000, 2000);
					
					//reset colors
					d.setColorA(heatmapColorA);
					d.setColorB(heatmapColorB);
				}
			}
			
			
			
			r.genio.clear();
			double[] returner = combineArrays(genio, percentStates);
			return returner;
	}
	public static double[] combineArrays(double[] array1, double[] array2) 
	{
        int length1 = array1.length;
        int length2 = array2.length;

        // Create a new array with enough space for elements from both arrays
        double[] combinedArray = new double[length1 + length2];

        // Copy elements from the first array
        System.arraycopy(array1, 0, combinedArray, 0, length1);

        // Copy elements from the second array
        System.arraycopy(array2, 0, combinedArray, length1, length2);

        return combinedArray;
    }
	public ArrayList<String> stateSearch(File file, SingleReader r)
	{
		int shell = 2;
		int requiredCalls = 5;
		
		for(int i = 0; i < nodeSystem.size(); i++)
		{
			MiniNode mini = nodeSystem.get(i);
			mini.counting.clear();
			mini.weight = 0;
		}
		
		BufferedReader br = null;
		ArrayList<Integer> inds = new ArrayList<Integer>();
		String line = "";
		int foundCount=0;
    	int tryCount=0;
		try 
		{
			br = new BufferedReader(new FileReader(file));
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	    try {
            line = br.readLine();

            while (line != null) 
	        {
	        	String[] s = new String[1];
	        	s[0] = line;
	        	if(line.contains(","))
	        	{
	        		s = line.split(",");
	        	}
//	        	if(s.length>0)
//	        		System.out.println(s[0] + " " + s.length);
	        	for(int k = 0; k<s.length; k++)
	        	{
	        		if(s[k].length()>1)
	        		{
		        		r.geneCheck(s[k], shell);//for trim analysis
	        			tryCount++;
		        		boolean found = false;
			        	for(int i = 0; i< d.dataPoints.size(); i++)
			        	{
			        		if(d.dataPoints.get(i).name.equalsIgnoreCase(s[k].trim()))
			        		{
			        			found = true;
			        			foundCount++;
			        			d.dataPoints.get(i).myMini.counting.add(d.dataPoints.get(i));
			        			inds.add(i);
			        		}
			        		
			        	}
		        	}
	        	}
	        	line = br.readLine();
        	}
//            System.out.println("searched for " + tryCount + " genes, found "+foundCount);   
		    } catch (IOException e)
		    {
				e.printStackTrace();
			}
	    	
	 
	  		ArrayList<String> trimmedList = r.countMultipleRepeatedGenes(requiredCalls, false); //for trim analysis
			
			//Random assignment of the same counts and weights to nodes
			for(int p = 1; p<g[0].length; p++){
				for(int z=0; z<foundCount; z++){
					int ree = (int) (Math.random()*bins.size());
					ree = nodeList.indexOf(bins.get(ree).myMini);
					g[ree][p] ++;
				}
			}
		
			//Actual assignment of the weight to all (count) assigned nodes
			for(Integer i: inds)
			{
				bins.get(i).myMini.counting.add(bins.get(i));
				bins.get(i).myMini.weight ++;
				g[nodeList.indexOf(bins.get(i).myMini)][0] ++;

			}
			
			giniWithWeightedBaseline();
			gini();
			pVal();
			double[] genio = new double[4];
			int trimmedCount = trimmedList.size();
			genio[0] = tryCount;
			genio[1] = foundCount;
			genio[2] = trimmedCount;
			genio[3] = trimmedCount / (double)foundCount*Math.log(foundCount);
			double zzz = Math.abs(g[nodes+1][1]);
			//System.out.println(zzz);

			
			/**Here is where I'll add the code to only capture sets that include a pre-defined genetic state*/
			
			if(pics && foundCount > 15 && foundCount < 3500 && ((zzz >= 3 && genio[3] > 3) || (zzz>=5) || genio[3]>7)) //ad hoc cluster calling metrics (rare enough gini and tuned ratio or exceptional individual value)
			{
//				d.equalWeight = equalWeight;
//				d.search(file);
//				d.repaint();
//				d.saveImg(file.getName(), outDir, 2000, 2000);
				if(trimmedPic)
				{
					outDir.mkdir();
					for(int i = 0; i < nodeSystem.size(); i++)
					{
						MiniNode mini = nodeSystem.get(i);
						mini.counting.clear();
						mini.weight = 0;
					}
					
					//flip colors for the trimmed maps for easier viewing
					d.setColorA(heatmapColorB);
					d.setColorB(heatmapColorA);
					
					//search the Trimmed List of genes
					d.search(trimmedList);
					d.repaint();
					d.saveImg(file.getName() + "_STATE_", outDir, 2000, 2000);
					
					//reset colors
					d.setColorA(heatmapColorA);
					d.setColorB(heatmapColorB);
				}
			}
			
			
			
			r.genio.clear();
			return trimmedList;
	}
	public String coClustering(File file1, File file2)
	{
		return "";
	}
	public void pVal()
	{
		double[] pValFind = new double[g[nodes].length-1];
		double mean = 0;
		for(int i = 0; i< pValFind.length; i++)
		{
			mean+= g[nodes][i+1];
			pValFind[i] = g[nodes][i+1];
		}
		mean/=pValnVal;
		double std = 0;
		for(int i = 0; i< pValFind.length; i++)
		{
			//if(i==0)
			//	System.out.println(g[nodes][i+1]+"\t"+mean+"\t"+giniWLArea+"\t"+g[nodes][0]);
			std += (g[nodes][i+1] - mean)*(g[nodes][i+1] - mean);
		}
		if(std!=0)
		{	//System.out.println(std);
			std/=pValnVal;
			std = Math.sqrt(std);
			g[nodes+1][1] = (giniWLArea-mean)/std;
		}
	}
	public void degreesOfSep()
	{
		double[][] dist = new double[nodes][nodes];
		for (int i = 0; i<nodes; i++)
		{
			int xCoordA = i%xNodes;
			int yCoordA = i/yNodes;
			for (int j = 0; j<nodes; j++)
			{
				int xCoordB = j%xNodes;
				int yCoordB = j/yNodes;
				
				double ddd = hexDist(xCoordA, yCoordA, xCoordB, yCoordB);
				dist[i][j] = ddd;
			}
		}
		
		for(int i = 0; i<g[0].length; i++)
		{
			double count = 0;
			for(int j = 0; j<nodes; j++)
			{
				for(int k = 0; k<nodes; k++)
				{
					sep[i] += dist[j][k]*g[j][i]*g[k][i];
					count+= g[j][i]*g[k][i];
				}
			}
			sep[i] /= count;
		}
	}
	public double degreesP()
	{
		double slit = 0;
		double avg = 0;
		for(int i = 1; i < sep.length; i++)
		{
			avg += sep[i];
		}
		avg /= pValnVal;
		double std = 0;
		for(int i = 1; i < pValnVal+1; i++)
		{
			std += (sep[i]-avg)*(sep[i]-avg);
		}
		if(std==0)
			std = 0.00001;
		std /= pValnVal;
		std = Math.sqrt(std);
//		NormalDistribution norm = new NormalDistribution(avg,std);
//		slit = norm.cumulativeProbability(sep[0]);
		return slit;
	}
	public double hexDist(int x1, int y1, int x2, int y2)
	{
		if(x1==x2 && y2 == y1)
			return 0;
		else
		{
			int hm = 0;
			int vm = 0;
			int cm = 0;
			boolean right = x1-x2<0;
			hm = Math.abs(x1-x2);
			vm = Math.abs(y1-y2);
			if(xNodes-x1<xNodes-x2)
				hm = Math.min(xNodes-x1 + x2,hm);
			else
				hm = Math.min(xNodes-x2 +x1, hm);
			
			if(yNodes-y1<yNodes-y2)
				vm = Math.min(yNodes-y1 + y2,vm);
			else
				vm = Math.min(yNodes-y2 +y1, vm);
			
			if(vm %2 !=0)
			{
				if(!right && y1%2 == 0)
					cm = Math.min(hm, (int)(((double)vm)/2 +.5));
				else if(right && y1%2 == 1)
					cm = Math.min(hm, (int)(((double)vm)/2 +.5));
			}
			else
				cm = Math.min(hm, (int)(((double)vm)/2));
			return hm+vm-cm; 
		}
	}
	public void giniWithWeightedBaseline()
	{
		// ~equality baseline
		double eArea = 0;
		for(int num = 1; num<=pValnVal; num++)
		{
			double maxWeight = 0;
			for(int i = 0; i<nodeList.size(); i++)
			{
				if(g[i][num]>maxWeight)
					maxWeight = (int) g[i][num];
			}
			double[] equal = new double[(int)maxWeight+1];
			for(int i=0; i<nodeList.size(); i++)
			{
				int o = (int) g[i][num];
				equal[o]++;
			}
			
			int eQTotalWealth = 0;
			for(int i=0; i<equal.length; i++)
			{

				//if(num==1 && equal[i]>0)
				//	System.out.println(i+"\t"+equal[i]);
				eQTotalWealth += i*equal[i];
			}
			
			double[] eQCumPop = new double[(int)maxWeight+1];
			double[] eQCumWealth = new double[(int)maxWeight+1];
			
			eQCumPop[0] = equal[0]/(double)nodeList.size();
			eQCumWealth[0] = 0;
			if(num==1  && printGiniCurves)
				System.out.println("Random assignment curve");
			for(int i = 1; i<equal.length; i++)
			{
				eQCumPop[i] = equal[i]/(double)nodeList.size() + eQCumPop[i-1];
				eQCumWealth[i] = (equal[i]*i)/eQTotalWealth + eQCumWealth[i-1];
				if(num==1 && printGiniCurves)
					System.out.println(eQCumPop[i]+"\t"+eQCumWealth[i]);
			}		
			//Equality area
			for(int i=0; i<equal.length-1; i++)
			{
				g[nodes][num]+=(((eQCumPop[i+1]-eQCumPop[i])*eQCumWealth[i])+((eQCumPop[i+1]-eQCumPop[i])*eQCumWealth[i+1]))/2;
				eArea += (((eQCumPop[i+1]-eQCumPop[i])*eQCumWealth[i])+((eQCumPop[i+1]-eQCumPop[i])*eQCumWealth[i+1]))/2;
			}
		}
		eArea /= pValnVal;		
		
		//inequality
		double maxWeight = 0;
		for(int j =0; j<nodeList.size(); j++)
		{
			if(g[j][0] > maxWeight)
				maxWeight=g[j][0];
		}
		
		/**
		//casting weights to ints may put all to zero in some data types?
		**/
		
		int[] unequal = new int[(int)maxWeight+1];
		for(int i1=0; i1<nodeList.size(); i1++)
		{
			int o = (int)g[i1][0];
			unequal[o] ++;
		}
		
		double uETotalWealth = 0;
		
		for(int i1=0; i1<unequal.length; i1++)
		{
			uETotalWealth += i1*unequal[i1];
		}

		double[] uECumPop = new double[unequal.length];
		double[] uECumWealth = new double[unequal.length];
		
		uECumPop[0] = unequal[0]/(double)nodeList.size();
		uECumWealth[0] = 0;
		
		if(printGiniCurves)
			System.out.println("\nObserved assignment curve");
		for(int i1 = 1; i1<unequal.length; i1++)
		{
			uECumPop[i1] = unequal[i1]/(double)nodeList.size() + uECumPop[i1-1];
			
			uECumWealth[i1] = (unequal[i1]*i1)/uETotalWealth + uECumWealth[i1-1];
			
			if(printGiniCurves)
				System.out.println(uECumPop[i1]+"\t"+uECumWealth[i1]);
		}
		
		//Lorenz area
		double lArea = 0;
		for(int i1=0; i1<unequal.length-1; i1++)
		{
			lArea += (((uECumPop[i1+1]-uECumPop[i1])*uECumWealth[i1])+((uECumPop[i1+1]-uECumPop[i1])*uECumWealth[i1+1]))/2;
		}
		//difference in area
		double dArea = eArea - lArea;
		gini = dArea / (lArea + dArea);
		g[nodes][0] = gini;
		giniWLArea = lArea;
		gini=0;
	}
	public void gini()
	{
		// ~equality baseline
		int mdp = nodeList.get(0).bins.size();
		for(int i = 0; i<nodeList.size(); i++)
		{
			if(nodeList.get(i).bins.size()>mdp)
				mdp = (int) (nodeList.get(i).bins.size());
		}
		double[] equal = new double[mdp+1];
		for(int i=0; i<nodeList.size(); i++)
		{
			int o = nodeList.get(i).bins.size();
			equal[o]++;
		}
		
		int eQTotalWealth = 0;
		for(int i=0; i<equal.length; i++)
			eQTotalWealth += i*equal[i];
		
		double[] eQCumPop = new double[mdp+1];
		double[] eQCumWealth = new double[mdp+1];
		eQCumPop[0] = equal[0]/nodeList.size();
		eQCumWealth[0] = 0;
		for(int i = 1; i<equal.length; i++)
		{
			eQCumPop[i] = equal[i]/nodeList.size() + eQCumPop[i-1];
			eQCumWealth[i] = (equal[i]*i)/eQTotalWealth + eQCumWealth[i-1];
		}		
		//Equaliy area
		double eArea = 0;
		for(int i=0; i<equal.length-1; i++)
			eArea += (((eQCumPop[i+1]-eQCumPop[i])*eQCumWealth[i])+((eQCumPop[i+1]-eQCumPop[i])*eQCumWealth[i+1]))/2 ;
		
		//inequality
		double maxWeight = 0;
		for(int j =0; j<nodeList.size(); j++)
			if(g[j][0] > maxWeight)
				maxWeight=g[j][0];
		
		/**
		//casting weights to ints may put all to zero in some data types?
		**/
		
		int[] unequal = new int[(int)maxWeight+1];
		for(int i1=0; i1<nodeList.size(); i1++)
		{
			int o = (int)g[i1][0];
			unequal[o] ++;
		}
		
		double uETotalWealth = 0;
		
		for(int i1=0; i1<unequal.length; i1++)
		{
			//if(unequal[i1]>0)
				//System.out.println("unequal:" + unequal[i1]);
			uETotalWealth += i1*unequal[i1];
		}

		double[] uECumPop = new double[unequal.length];
		double[] uECumWealth = new double[unequal.length];
		
		uECumPop[0] = unequal[0]/(double)nodeList.size();
		uECumWealth[0] = 0;
		
		for(int i1 = 1; i1<unequal.length; i1++)
		{
			uECumPop[i1] = unequal[i1]/(double)nodeList.size() + uECumPop[i1-1];
			
			uECumWealth[i1] = (unequal[i1]*i1)/uETotalWealth + uECumWealth[i1-1];
		}
		
		//Lorenz area
		double lArea = 0;
		for(int i1=0; i1<unequal.length-1; i1++)
		{
			lArea += (((uECumPop[i1+1]-uECumPop[i1])*uECumWealth[i1])+((uECumPop[i1+1]-uECumPop[i1])*uECumWealth[i1+1]))/2;
		}	
		//difference in area
		double dArea = eArea - lArea;
		gini = dArea / (lArea + dArea);
		g[nodes+1][0] = gini;
		gini=0;
	}
	
	
	public void mapReader(String f)
	{
		ArrayList<String> StringMat = new ArrayList<String>();
		try 
		{
			@SuppressWarnings("resource")
			Scanner in = new Scanner(new FileReader(f));
			String sizer = in.next();
			int xo = Integer.parseInt(sizer.substring(0,sizer.indexOf("x")));
			int yo = Integer.parseInt(sizer.substring(sizer.indexOf("x")+1,sizer.length()));
			in.next();
			while(in.hasNext())
			{
				StringMat.add(in.next());
			}
			createNoder(StringMat, xo, yo);
		} catch (FileNotFoundException e) {e.printStackTrace();}
	}
	
	public void createNoder(ArrayList<String> strings, int x, int y)
	{
		nodeList = new ArrayList<MiniNode>(x*y);
		bins = new ArrayList<DataPoint>();
		MiniNode current = null;
		xNodes = x;
		yNodes = y;
		for(String whole: strings)
		{
			if(whole.contains("["))
			{
				current = new MiniNode(whole);
				nodeList.add(current);
			}
			else
				current.addDP(whole);
		}
		for(int i = 0; i<nodeList.size(); i++)
		{
			for(int j = 0; j < nodeList.get(i).bins.size(); j++)
			{
				bins.add(nodeList.get(i).bins.get(j));
			}
		}
		nodeSystem = new MiniSystem(nodeList,xNodes,yNodes);
	}
	
	public void heatMapping()
	{
		minBins = (int) (nodeList.get(0).bins.size());
		maxBins = (int) (nodeList.get(0).bins.size());
		for(int i = 0; i<nodeList.size(); i++)
		{
			if(nodeList.get(i).bins.size()>maxBins)
				maxBins = (int) (nodeList.get(i).bins.size());
			if(nodeList.get(i).bins.size()<minBins)
				minBins = (int) (nodeList.get(i).bins.size());
		}
		if(maxBins == minBins) maxBins++;
	}
	public void writeFile(ArrayList<String> writer)
	{
		BufferedWriter b;
		try 
		{
			FileWriter ff;
			if(pics)
				ff = new FileWriter(outDir+File.separator+"lorenz_analysis.txt",false);  
			else
				ff = new FileWriter("lorenz_analysis.txt",false);
			b = new BufferedWriter(ff);
			PrintWriter printer = new PrintWriter(b);
			for(int i = 0; i< writer.size(); i++)
			{
				printer.println(writer.get(i));
			}
			printer.close();
		}
		catch (IOException e) 
		{
			e.printStackTrace();
		}
	}
}

