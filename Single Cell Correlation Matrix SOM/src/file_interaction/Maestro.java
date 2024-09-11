package file_interaction;

//Is this class used at all?

public class Maestro 
{
	public static void main(String[] s)
	{
		CorrelationMatrixReader r = new CorrelationMatrixReader(s[0], s[1]);
		double[][] matrix = r.getMatrix();
		String[] genes = r.getGenes();
		for(int i = 0; i<matrix.length; i++)
		{
			System.out.print(genes[i] +"\t");
			for(int j = 0; j < matrix[i].length; j++)
			{
				System.out.print(matrix[i][j] + "\t");
			}
			System.out.println();
		}
	}
}
