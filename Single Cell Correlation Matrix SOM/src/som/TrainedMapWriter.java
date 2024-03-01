package som;

public class TrainedMapWriter 
{
	MultiThreadTraining som;
	public TrainedMapWriter(MultiThreadTraining t, Map m)
	{
		som = t;
		double[] r = t.finishTraining();
		System.out.println("Percent unoccupied = " + r[2] + "\t Pearson Quality = " + r[0] + "\t Cosine Quality = " + r[1]);
	}
	public void write()
	{
		som.writeFile();
	}
	
}
