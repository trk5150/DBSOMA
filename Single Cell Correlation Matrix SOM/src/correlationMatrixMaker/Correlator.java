package correlationMatrixMaker;

import java.util.Arrays;

public class Correlator 
{
	Correlator()
	{
		
	}
	public double[][] buildCorrelationMatrix(double[][] parsedCountTranspose) //counts is in the transposed form of [genes][cells]
	{
		int genes = parsedCountTranspose.length;
		double[][] cors = new double[genes][genes];
		
		for(int i = 0; i<genes; i++)
		{
			for (int j = i; j<genes; j++)
			{
				if(i == j)
					cors[i][j] = 1;
				else
				{
					double pears = pearson(parsedCountTranspose[i], parsedCountTranspose[j]);
					cors[i][j] = pears;
					cors[j][i] = pears;
				}
			}
			System.out.println("count: " + i);
		}
		
		return cors;
	}
	public double pearson(double[] parsedCountTranspose, double[] parsedCountTranspose2)
	{
		double[] nnn = parsedCountTranspose;
		double[] ddd = parsedCountTranspose2;
		double totaldd = 0;
		double totalnn = 0;
		for(int i = 0; i<nnn.length; i++)
		{
			totalnn += nnn[i];
		}
		for(int i = 0; i<ddd.length; i++)
		{
			totaldd+=ddd[i];
		}
		for(int i = 0; i<nnn.length; i++) {nnn[i]/=totalnn;}
		for(int i = 0; i<ddd.length; i++) {ddd[1]/=totaldd;}
		
		double sx = 0.0;
	    double sy = 0.0;
	    double sxx = 0.0;
	    double syy = 0.0;
	    double sxy = 0.0;

	    int n = nnn.length;
	    for(int i = 0; i < n; i++) 
	    {
	      double x = nnn[i];
	      double y = ddd[i];

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
	    //System.out.println(rr);
	    return rr;
	}
}
