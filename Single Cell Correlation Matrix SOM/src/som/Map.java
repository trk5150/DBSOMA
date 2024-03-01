package som;

public class Map 
{
	boolean debug;
	int xNodes, yNodes, nodeNum;
	int farthestNeighbor,maxDist;
	int[] dpcount;
	public int[][] neighborly;
	double zero, sgm;
	double[] dpMag, nMag, weightsByHood;
	
	public Map(int x, int y, boolean b, double sig)
	{
		sgm = sig;
		xNodes = x;
		yNodes = y;
		debug = b;
		zero = Math.pow(10, -7);
		nodeNum = xNodes * yNodes;
	}
	public void initializeMap()
	{

		//nodeNum = xNodes * yNodes; moved to constructor 5/3/2023, hope nothing broke <3
		farthestNeighbor = 0;
		dpcount = new int[nodeNum];
		nMag = new double[nodeNum];
		assignNeighbors();
		
		if(debug)
		{
			for(int i = 0; i< xNodes; i++)
			{
				for(int j = 0; j < yNodes; j++)
				{
					System.out.print(neighborly[i][j]+ "\t");
				}
				System.out.println();
			}
		}
	}
	
	//Assigns neighbors based on hex dist
	public void assignNeighbors() 								
	{
		neighborly = new int[nodeNum][nodeNum];
		for(int i = 0; i< nodeNum; i++)
		{
			int xA = i%xNodes;
			int yA = i/yNodes;
			for(int j = i; j < nodeNum; j++)
			{
				int xB = j%xNodes;
				int yB = j/yNodes;
				int dist = hexDist(xA, yA, xB, yB);
				neighborly[i][j] = dist;
				if(dist > farthestNeighbor)
					farthestNeighbor = dist;
				neighborly[j][i] = neighborly[i][j];
			}
		}
		maxDist = farthestNeighbor;
		weightsByHood = new double[farthestNeighbor];
	}
	
//	public int nFactorSize()
//	{
//		int i = 0;
//		double sig = sgm;
//		double lr = 1;
//		double r = 1;
//		while(r > zero)
//		{
//			r = lr*Math.exp(-1*((i*i)/(2*sig*sig)));
//			if(r>zero)
//				i++;
//		}
//		return i;
//	}
	public double degreesOfSep()
	{
		int nodes = nodeNum;
		double sepp = 0;
		
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
			double count = 0;
			for(int j = 0; j<nodes; j++)
			{
				for(int k = 0; k<nodes; k++)
				{
					sepp += dist[j][k]*dpcount[k]*dpcount[j];
					count+= dpcount[k]*dpcount[j];
				}
			}
			sepp /= count;
			return sepp;		
	}
	public int hexDist(int x1, int y1, int x2, int y2)
	{
		//System.out.println(xNodes + "  " + yNodes);
		
		if(x1 == x2 && y1 == y2)
			return 0;
		else
		{
			int hm = 0;
			int vm = 0;
			int cm = 0;
			boolean right = x1-x2<0;
			
			hm = Math.abs(x1-x2);
			vm = Math.abs(y1-y2);
			if(xNodes-x1<xNodes-x2 && xNodes-x1 + x2 < hm)
			{
				right = true;
				hm =  xNodes-x1 + x2;	
			}
			else if (xNodes-x1>xNodes-x2 && xNodes-x2 +x1 < hm)
			{
				right = false;
				hm = xNodes-x2 +x1;
			}
			
			if(yNodes-y1<yNodes-y2 && yNodes-y1 + y2 < vm)
			{
				vm = yNodes-y1 + y2;
			}
			else if(yNodes-y1>yNodes-y2 && yNodes-y2 + y1 < vm)
			{
				vm = yNodes-y2 + y1;
			}
			
			if(vm %2 !=0)
			{
				if(!right && y1%2 == 0)
				{
					cm = (int) Math.min(hm, (((double)vm)/2 +.5));
				}
				else if(!right && y1%2 == 1)
				{
					cm = (int) Math.min(hm, (((double)vm)/2));
				}
				else if(right && y1%2 == 0)
				{
					cm = (int) Math.min(hm, (((double)vm)/2));
				}
				else //if(right && y1%2 == 1)
				{
					cm = (int) Math.min(hm, (((double)vm)/2 +.5));
				}
			}
			else
				cm = (int) Math.min(hm, ((((double)vm)/2)));
			return hm+vm-cm; 
		}
	}
}
