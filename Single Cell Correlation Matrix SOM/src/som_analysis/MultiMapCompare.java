package som_analysis;

import java.util.ArrayList;

public class MultiMapCompare 
{
	public int shells;
	public String[] mapFiles;
	public String geneOfInterest;
	ArrayList<SingleReader> soms;
	public String[] genes;
	public int[][] dists; //[gene index][file]  = distance
	
	public MultiMapCompare(String g, int s, String[] f)
	{
		shells = s;
		geneOfInterest = g;
		mapFiles = f;
		
		soms = new ArrayList<SingleReader>();
		for(int i = 0; i < mapFiles.length; i++)
		{
			SingleReader r = new SingleReader(false);
			r.readSOM(mapFiles[i]);
			r.mapMaker();
			r.geneCheck(geneOfInterest, shells);
			if(r.notFound)
				System.err.print(g + " not found in file " + mapFiles[i]);
			soms.add(r);
		}
	}
	
	public ArrayList<String[][]> indexGenes()
	{
		ArrayList<String[][]> geneIndexes = new ArrayList<String[][]>();
		for(int i = 0; i<soms.size(); i++)
		{
			geneIndexes.add(soms.get(i).indexGenesByDistance()); //returns String[i][0] = name String[i][1] = distance
		}
		return geneIndexes;
	}
	
	public void unifyGeneDistances(ArrayList<String[][]> gAndDs)
	{
		ArrayList<String> geno = new ArrayList<String>();
		ArrayList<Double> disto = new ArrayList<Double>();
		ArrayList<String> inBothList = new ArrayList<String>();
		int inBoth = 0;
		//adds all genes on the first file
		for(int i = 0; i<gAndDs.get(0).length; i++)
		{
			geno.add(gAndDs.get(0)[i][0]);
			//disto.add(Integer.getInteger(gAndDs.get(0)[i][1]));
		}
		
		//checks and adds non-duplicate genes for the rest of the input files
		for(int i = 1; i<gAndDs.size(); i++)
		{
			for(int j = 0; j<gAndDs.get(i).length; j++)
			{
				String checkGene = gAndDs.get(i)[j][0];
				//System.out.println(checkGene);
				Boolean addIt = true;
				if(checkGene.length()>0) 
				{
					for(int k = 0; k < geno.size(); k++)
					{
						int x = gAndDs.get(0).length-1;
						if(checkGene.equalsIgnoreCase(geno.get(k))&&k<x)
						{
							addIt = false;
							inBothList.add(checkGene);
							//only works b/c geno.get(k) for setting the checked gene is the same index as gsAndDs.get(0)[k][1], because of the initial matrix filling loop above
							//will break if more than 2 maps are used, will need to fix this. Sort of hacked together for testing purposes
							//will also break if the equals is satisfied while k > gAndDs size, which can happen, but not sure why
							
							double dist1 = Double.parseDouble(gAndDs.get(i)[j][1]);
							double dist2 = Double.parseDouble(gAndDs.get(0)[k][1]);
							disto.add(dist1+dist2);
			
							inBoth++;
						}
					}
				}
				if(addIt)
					geno.add(checkGene);
			}
		}
		for(int i = 1; i<disto.size(); i++)
		{
			double avg = disto.get(i)/(double)gAndDs.size();
			disto.set(i, avg);
		}
		System.out.println(inBoth + " genes in files within " + (shells-1) + " shells of " + geneOfInterest);
		for(int i = 0; i < inBothList.size(); i++)
		{
			boolean numbies = false;
			double dist = disto.get(i);
//			double weightedDist = (double)Math.pow(1.7, -dist);
//			weightedDist*=1000;
			if(numbies)
				System.out.println(inBothList.get(i) + "\t" + dist);// + "\t" + weightedDist);
			else
				System.out.println(inBothList.get(i));
		
		}
		boolean singles = true;
		if(singles)
		{
			System.out.println(geno.size() + " genes observed in only 1 dataset within "  + shells + " shells");
			for (int i = 0; i<geno.size(); i++)
			{
				System.out.print(geno.get(i) + ", ");
			}
		}
	}
	public static void main(String[] args)
	{
		String gene = args[0];
		int shl = Integer.parseInt(args[1]);
		String[] maps = new String[args.length-2];
		for(int i =0; i< maps.length;i++)
		{
			maps[i] = args[i+2];
		}
		
		
		MultiMapCompare mm = new MultiMapCompare(gene, shl, maps);
		
		ArrayList<String[][]> geneIndexes = mm.indexGenes();
		
		mm.unifyGeneDistances(geneIndexes);
		/**Note to self: probably the best way to generate Venn diagrams for more than 2 input SOMs will be to generate pairwise "inBoth" gene lists, then sub-compare those smaller lists
		 * then for avg distances, go top down, check only those avg distance included in all, then all-1...all-n until you reach all pairwise dists
		 * */
		
		
	}
}
