package som;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class MultiThreadTraining 
{
	IterateThread[] iThreads;
	AssignThread[] aThreads;
	boolean debug, cos;
	double[][] dpMatrix, nodes;
	String[] genes;
	Map map;
	int threadNum;
	int[][] threadLims;
	boolean[][] assignments;
	int[] dpcount;
	int dpSize, iterations;
	double[] dpMag, nMag, weightsByHood;
	double sgm, sgmStop;
	public double pearsonQuality;
	public double cosineQuality;
	
	public MultiThreadTraining(Map mm, String[] g, double[][] mat, boolean db, int threadNo)
	{
		sgm =3;
		sgmStop = 0.2;
		iterations = 100;
		cos = true;
		debug = db;
		map = mm;
		genes = g;
		dpMatrix = mat;
		dpSize = mat.length;
		threadNum = threadNo;
		iThreads = new IterateThread[threadNo];
		aThreads = new AssignThread[threadNo];
		weightsByHood = map.weightsByHood;
	}
	public class IterateThread extends Thread
	{
		public double[][] partOfNeighbors;
		public int thisT;
		public IterateThread(int i)
		{
			super();
			thisT = i;
			partOfNeighbors = new double[threadLims[i][3]-threadLims[i][2]+1][nodes[0].length];
		}
		public void run()
		{
			for(int i = threadLims[thisT][2]; i <= threadLims[thisT][3]; i++) 										//Updates each Node
			{
				double[] numVec = new double[nodes[i].length]; 						//initialize a vector for sum of all datapoints*weight for weighted average
				//System.out.print(nodes[i].length);
				double sumDenom = 0;												//SumDenom = sum of weight factors used as denominator in weighted average
				
				for(int j = 0; j < map.neighborly[i].length; j++)        				//Cycles through neighborhood  (all nodes) for each node
				{
					//System.out.print(map.neighborly[i][j]);
					if(map.neighborly[i][j]<map.maxDist)
					{
						double nj =dpcount[j];
						double h = weightsByHood[map.neighborly[i][j]]; 							//neighborly = degree of separation between nodes				
						double weight= nj*h;
						//System.out.println(weight + "  " + h +" " + nj +  "  ");
						
						for(int k = 0; k<dpSize; k++)									//Cycles through the data points assigned to each Node in each shell
						{
							if(assignments[k][j])										//if given node has given data point assigned to it, calculate:
							{
								for(int l = 0; l<dpMatrix[k].length; l++)
								{
									numVec[l] += dpMatrix[k][l]*weight;
								}
								sumDenom += weight;
							}
						}
					}
				}
				for(int l = 0; l<numVec.length; l++) 						//Takes weighted average = updated vector for node
				{
					numVec[l] = numVec[l]/sumDenom;
				}
				partOfNeighbors[i-threadLims[thisT][2]] = numVec;								//Node vectors are updated in order, but reassignment happens after all are updated -> batched SOM
			}
		}		
	}
	public class AssignThread extends Thread
	{
		public int[] pp;
		public boolean[][] pAssign;
		public int thisT;
		public AssignThread(int i)
		{
			super();
			thisT = i;
			pp = new int[map.nodeNum];
//			System.out.println("nodeNum = " + map.nodeNum);
			pAssign = new boolean[threadLims[i][1]-threadLims[i][0]+1][map.nodeNum];
		}
		public void run()
		{

			int errCount = 0;
			for(int i = threadLims[thisT][0]; i <= threadLims[thisT][1]; i++)
			{
				double initdist = -1;
				
				if(cos)
					initdist = cosineSim(0, i);
				else
					initdist = pearson(0,i);
				
				int bestNode = 0;
				for(int j = 0; j < map.nodeNum; j++) 					//Each data point goes through each Node looking for most similar
				{
					double dist = 0;
					if(cos)
						dist = cosineSim(j,i);
					else
						dist = pearson(j,i);
					
					if(dist>initdist)  								//Node with highest cosine similarity value to data point = most similar
					{
						initdist = dist;
						bestNode = j;
					}
				}
				if(initdist == 0)
				{
					errCount++;
//					System.err.println("error initdist==0; thisThread =  " + thisT + "\t threadLims " + threadLims[thisT][0] + ", " + threadLims[thisT][1] + "\t i = " + i + "\t bestNode= " + bestNode);
				}
				pAssign[i-threadLims[thisT][0]][bestNode] = true;
				pp[bestNode]++;
				//System.out.println("error intidist == 0 count: "+errCount);
			}
		}
	}
	
	//assigns datapoints to threads
		//Node vectors are each set equal to a randomly selected data point's vector
	public void initialize() 
	{
		int nodeNum = map.nodeNum;
		System.out.println("initial node num = " + map.nodeNum + " dpSize = " + dpSize + " ");
		nodes = new double[nodeNum][dpSize];
		assignments = new boolean[dpSize][nodeNum];    							/** Be Careful, this has opposite i x j convention (dp by node instead of node by dp)*/
		dpcount = new int[nodeNum];
		nMag = new double[nodeNum];
		
		for(int i = 0; i<nodeNum; i++)
			nodes[i] = dpMatrix[(int)(Math.random()*dpSize)];
		threadLims = new int[threadNum][4];
		int dpNum = (dpSize / threadNum)-1;
		int dpExtra = dpSize % threadNum;
		int nNum = (nodeNum / threadNum)-1;
		int nExtra = nodeNum % threadNum;
		int curN = 0;
		int curD=0;
		for(int i = 0; i< threadNum; i++)
		{
			threadLims[i][0]=curD;
			curD += dpNum;
			if(dpExtra > i)
				curD++;
			threadLims[i][1]=curD;
			
			threadLims[i][2]=curN;
			curN += nNum;
			if(nExtra > i)
				curN++;
			threadLims[i][3]=curN;

			curN++;
			curD++;
		}
		dpMag = new double[dpSize];
		for(int i = 0; i< dpSize; i++)
		{
			dpMag[i] = updateMag(dpMatrix[i]);
			if(dpMag[i]==0)
				System.err.println("d  "+dpMag[i]);
		}
		assignNodes();
	}
		
	public double updateMag(double[] g)
	{
		double magD = 0;
		for(int i = 0; i<g.length; i++)
		{
			if(Double.isNaN(g[i]))
				System.out.println("g[i] = NaN at " + i);
			magD += (g[i])*(g[i]);
//			if(Double.isNaN(magD))
//			{
//				System.out.println(magD + " = " + g[i] + "^2 \t from i =" +i);
//			}
			
		}
		magD = Math.sqrt(magD);
		return magD;
	}
	//Cosine similarity
	public double cosineSim(int iNode, int jDP)							//Cosine Similarity = (A  B)/(||A||*||B||)
	{																			//Higher values indicate more similarity (ranges from [-1, 1])
		double dot = 0;
		double magN = nMag[iNode];
		double magD = dpMag[jDP];
		for(int i = 0; i < nodes[iNode].length; i++)
		{
			dot += (nodes[iNode][i]) * (dpMatrix[jDP][i]);
		}
		dot = dot/(magN*magD);
		if(Double.isNaN(magN))
		{
			//System.err.println("NaNa  @"+ iNode+ magN + "   "+jDP + magD);
			dot = 0;
		}
		
		if(dot > 1 && dot< 1.0000001)
	    	dot=1;
	    
		return dot;
	}
	public double pearson(int iNode, int jDP)
	{
		double sx = 0.0;
	    double sy = 0.0;
	    double sxx = 0.0;
	    double syy = 0.0;
	    double sxy = 0.0;

	    int n = nodes[iNode].length;
	    for(int i = 0; i < n; i++) 
	    {
	      double x = nodes[iNode][i];
	      double y = dpMatrix[jDP][i];

	      sx += x;
	      sy += y;
	      sxx += x * x;
	      syy += y * y;
	      sxy += x * y;
	    }

	    // covariation
	    double cov = sxy / n - sx * sy / n / n;
	    // standard error of x
	    double sigmax = Math.sqrt(sxx / n -  sx * sx / n / n);
	    // standard error of y
	    double sigmay = Math.sqrt(syy / n -  sy * sy / n / n);

	    // correlation is just a normalized covariation
	    
	    double rr = cov / sigmax / sigmay;
	    if(rr > 1 && rr< 1.0000001)
	    	rr=1;
	    return rr;
	}
	
	public void iterater()
	{
		for(double i = 0; i < iterations; i++)
		{
			System.out.println(i);
			//This code writes SOMs at different stages in training
//			if(i%(iterations/4) == 0)
//			{
//				finishTraining();
//				writeFile((int)i);
//			}
			nFactor(i);
			iterate(i);
		}
	}
	
	//Finds the weight factor of a given data point based on learning rate and neighborhood function
	public void nFactor(double itercount) 
	{   
		double sig = sgm-sgmStop;
		sig = (sig - (sig*(itercount/iterations)))+sgmStop; 												//sigma = width of Guassian kernel (sgm ->1)
		double r = 0;
		double lrStop = 0.01;
		double lr = 1-lrStop;
		lr = (lr - (lr*(itercount/iterations)))+lrStop;
		for(double o = 0; o < weightsByHood.length; o++)
		{
			r = lr*Math.exp(-1*((o*o)/(2*sig*sig)));
			weightsByHood[(int)o] = r;
//			System.out.println(r);
		}
	}
	
	
	//Where the magic happens
	public void iterate(double itercount)
	{
		for(int i = 0; i< threadNum; i++)
		{
			IterateThread t = new IterateThread(i);
			iThreads[i] = t;
			iThreads[i].start();
		}
		for(Thread t : iThreads){try {t.join();} catch (InterruptedException e) {}}
		
		//stitch the threads 
		synchronized (this)
		{
			for(int i = 0; i < threadNum; i++)
			{
				for(int j = 0; j < iThreads[i].partOfNeighbors.length; j++)
				{
					nodes[threadLims[i][2]+j] = iThreads[i].partOfNeighbors[j];
				}
			}
		}

//		System.out.println(nodes[2][2]);
//		count++;
//		timerIterate += System.currentTimeMillis()-init;
		assignNodes();
	}
	public void assignNodes()
	{
		int nodeNum = map.nodeNum;
		//System.out.println(map.);
		for(int i = 0; i<nodeNum; i++) 						//calculate node magnitudes  
		{
			
			nMag[i] = updateMag(nodes[i]);
			if(debug)
			{
				//System.out.println(nodes[i][2]);
				//System.out.println(nMag[i]);
			}
			if(nMag[i]==0)
				System.out.println("nMag " + i + " = " + nMag[i]);
		}
		for(int i = 0; i< threadNum; i++)
		{
			AssignThread t = new AssignThread(i);
			aThreads[i] = t;
			aThreads[i].start();
		}
		for(Thread t : aThreads){try {t.join();} catch (InterruptedException e) {}}				//join all threads at completion
		assignments = new boolean[dpSize][nodeNum];
		dpcount = new int[nodeNum];
		
		synchronized (this)
		{
			for(int i = 0; i<threadNum; i++)
			{
				for(int j = 0; j < dpcount.length; j++)
				{
					dpcount[j] += aThreads[i].pp[j];
				}
				for(int j = 0; j < aThreads[i].pAssign.length; j++)
				{
					assignments[threadLims[i][0]+j] = aThreads[i].pAssign[j];
				}
			}
		}
		//timerAssign += System.currentTimeMillis()-init;
	}	
	
	public double[] finishTraining()       
	{
		//writeFile();
		double[] rtn = new double[3];
		findPearsonQuality();
		findCosineQuality();
		rtn[0] = pearsonQuality;
		rtn[1] = cosineQuality;
		rtn[2] = percentUnoccupied();
		return rtn;
	}
	//Quality is defined as the average similarity between dataPoints and their assigned nodes 
	public void findPearsonQuality()
	{
		pearsonQuality = 0;
		for(int i = 0; i < assignments.length; i++)					//dp i
		{
			for(int j =0; j< assignments[i].length; j++)			//node j
			{
				if(assignments[i][j])
					pearsonQuality += pearson(j, i);
			}
		}

		pearsonQuality /= dpSize;
	}
	public void findCosineQuality()
	{
		cosineQuality = 0;
		for(int i = 0; i < assignments.length; i++)					//dp i
		{
			for(int j =0; j< assignments[i].length; j++)			//node j
			{
				if(assignments[i][j])
					cosineQuality += cosineSim(j, i);
			}
		}
		
		cosineQuality /= dpSize;
	}
	public double percentUnoccupied()
	{
		double percent = 0;
		for (int i = 0; i< map.nodeNum; i++)
		{
			if (dpcount[i]== 0)
			{
				percent++;
			}
		}
		percent /= map.nodeNum;
		percent = ((double)((int)(percent*10000))/10000);
		return percent;
	}
	public void writeFile()
	{
		String filePrefix = "C:\\Users\\tik105\\Desktop\\mRNA\\single_cell_read_counts\\Stage 5\\SOM\\Stage5"; //"C:\\Users\\tik105\\Desktop\\SubCellTRIMMED_SC-BSom5050";
		BufferedWriter b;
		
		String somFilename = filePrefix+".som";
		String infoFilename = filePrefix+".info";
		String nodevecFilename = filePrefix+".nodevec";
		
		//SOM file
		try 
		{
			System.out.println("Printing SOM to "+somFilename);
						
			FileWriter ff = new FileWriter(somFilename,true);
			b = new BufferedWriter(ff);
			PrintWriter printer = new PrintWriter(b);
			
			printer.print(map.xNodes+"x"+map.yNodes+ "\nSigma:" + sgm + "->" + sgmStop +"\n");
			for(int i = 0; i<map.xNodes; i++)
			{
				for(int j = 0; j<map.yNodes; j++)
				{
					printer.print("["+j +","+i+"]\t"+"(");
					int count = 0;
					int nono = i*map.xNodes + j%map.xNodes;
					for(int k = 0; k< dpSize; k++)
					{	
						if(assignments[k][nono])
						{
							printer.print(genes[k]);
							count++;
							if(count < dpcount[nono])
								printer.print(", ");
						}
					}
					printer.print(")" + "\n");
				}
				printer.print("\n");
			}
			printer.print("\n");
			printer.close();
		}catch (IOException e){
			e.printStackTrace();
			System.err.print("Error printing SOM file "+somFilename);
		}
		
		//Info file
		try 
		{
			System.out.println("Printing info to "+infoFilename);
						
			FileWriter ff = new FileWriter(infoFilename,true);
			b = new BufferedWriter(ff);
			PrintWriter printer = new PrintWriter(b);
			
			String simMetric = cos ? "cosine" : "pearson";
			
			printer.print("Nodes:\t"+map.xNodes+"x"+map.yNodes+"\n" +
						"Sigma:\t" + sgm + "->" + sgmStop +"\n" +
						"Iterations:\t" + iterations +"\n" +
						"SimilarityMetric:\t" + simMetric +"\n" +
						"CosineQuality:\t" + cosineQuality +"\n" +
						"PearsonQuality:\t" + pearsonQuality +"\n" +
						"PercentUnoccupied:\t"+percentUnoccupied() +"\n");
			printer.close();
		
		}catch (IOException e){
			e.printStackTrace();
			System.err.print("Error printing info file "+infoFilename);
		}
		
//		try 
//		{
//			System.out.println("Printing node vectors to "+nodevecFilename);
//			
//			FileWriter ff = new FileWriter(nodevecFilename,true);
//			b = new BufferedWriter(ff);
//			PrintWriter printer = new PrintWriter(b);
//			
//			for(int i = 0; i < nodes.length; i++)
//			{
//				int x = i%map.xNodes;
//				int y = i/map.yNodes;
//				printer.print("("+ x+","+y +")\t");
//				for(int j = 0; j < nodes[i].length-1; j++)
//				{
//					printer.print(nodes[i][j]+"\t");
//				}
//				printer.print(nodes[i][nodes[i].length-1]+"\n");
//			}
//			printer.close();
//		}
//		catch (IOException e) 
//		{
//			e.printStackTrace();
//			System.err.print("Error printing node vector file "+nodevecFilename);
//		}
	}
}
