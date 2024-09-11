package mapDrawing;
import java.awt.Color;
import java.util.ArrayList;


public class SOMAnalysis
{
	public UseMap u;
	
	public SOMAnalysis()
	{	 

	}
	
	public void useSOM(String mapFileName, String searchDir, String outDir, boolean equalWeights, boolean flipColors, String[] states){
		System.setProperty("java.awt.headless", "true");
		UseMap u = new UseMap(mapFileName, searchDir, outDir, true, equalWeights, states);
		if(flipColors) {
			u.setColorA(Color.blue);
			u.setColorB(Color.red);
		}

		u.searchSystem();
	}
	
	public void useSOMpairwise(String mapFileName, String searchDir, String outDir, boolean equalWeights, boolean flipColors){
		System.setProperty("java.awt.headless", "true");
		UseMap u = new UseMap(mapFileName, searchDir, outDir, true, equalWeights, null);
		if(flipColors) {
			u.setColorA(Color.blue);
			u.setColorB(Color.red);
		}
		u.searchSystemDouble();
	}
	
//	public void corSOM(String mapFileName, String searchDir, String outDirectory){
//		System.setProperty("java.awt.headless", "true");
//		ArrayList<String> folders = new ArrayList<String>();
//		folders.add(searchDir);
//		Correlations c = new Correlations(mapFileName, folders, outDirectory);
//		c.searchSystem();
//	}
	public void SetScan(String mapFile, String searchDir, String outDir, String[] stateFiles)
	{
		SetScanner ss = new SetScanner(mapFile, searchDir, outDir, stateFiles);
		ss.searchSystem();
	}
	public static void main(String[] args)
	{
		SOMAnalysis s = new SOMAnalysis();
		Boolean colorFlip = false;
		if(args.length>2)
		{
			String[] states = new String[args.length-3];
			for(int i = 0; i< states.length; i++)
			{
				states[i] = args[i+3];
				//System.out.println(states[i]);
			}
			System.out.println("set scan");
			s.SetScan(args[0], args[1], args[2], states);
//			s.useSOM(args[0], args[1], args[2], true, colorFlip, states);
		}
		else
			s.useSOM(args[0], args[1], args[2], true, colorFlip, null);
	
	}
	
}

