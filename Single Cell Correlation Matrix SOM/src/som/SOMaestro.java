package som;

public class SOMaestro 
{
	public static int matSize = 50;
	
	public static void main(String[] args)
	{
		int threads = 4;
		double sig =3; //need to pass this to training as well as map
		int matrixSizeArg = matSize;
		
		if(args.length >2)
			threads = Integer.parseInt(args[2]);
		if(args.length > 3)
			matrixSizeArg = Integer.parseInt(args[3]);
		
		boolean debugReaders = false;
		InputMatrixReader m = new InputMatrixReader(args[0], debugReaders);
		double[][] matrix = m.readMat();
		System.out.println("finished reading matrix");
		
		InputGeneReader g = new InputGeneReader(args[1],debugReaders);
		String[] genes = g.readGenes(matrix.length);
		System.out.println("finished genes");
		
		boolean debugNeighbors = false;
		Map map = new Map(matrixSizeArg, matrixSizeArg, debugNeighbors, sig);
		map.initializeMap();
		System.out.println("finished building neighborhood");
		
		System.out.println("map stuff: " +map.xNodes + " " + map.yNodes + " " +map.farthestNeighbor + " " + map.maxDist + " " + map.sgm);
		
		
		boolean debugTrainingPt1 = false;
		MultiThreadTraining t = new MultiThreadTraining(map, genes, matrix, debugTrainingPt1, threads, args[0]);
		//System.out.println("pre init trainer stuff: " + t.dpSize + " " + t.dpMatrix[1][1] + " " + t.dpMatrix[0][1] + " " + t.dpMag[0] + " " + t.dpMag[1] + " " + t.nMag[0] + " " + t.nMag[1]);
		t.initialize();
		System.out.println("post init trainer stuff: " + t.dpSize + " " + t.dpMatrix[1][1] + " " + t.dpMatrix[0][1] + " " + t.dpMag[0] + " " + t.dpMag[1] + " " + t.nMag[0] + " " + t.nMag[1]);
		System.out.println("initialized");
		t.iterater();
		TrainedMapWriter tmw = new TrainedMapWriter(t, map);
		tmw.write();
	}
}
