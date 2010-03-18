package utils;

public class LinearAlgebraUtils {
	//Getting  the identity matrix
	public static double[][] getIdentity(int size) {
	    double[][] matrix = new double[size][size];
	    for(int i = 0; i < size; i++)
	      for(int j = 0; j < size; j++)
	        matrix[i][j] = (i == j) ? 1 : 0;
	    return matrix;
	  }
	//Exact matlab replica of tril.
	public static double[][] tril(double[][] mat, int diag){
		double[][] trils = new double[mat.length][]; 
		
		for (int i = 0; i < mat.length; i++) {
			double[] ds = mat[i];
			trils[i] = new double[ds.length];
			for (int j = 0; j < ds.length; j++) {
				if(j-i<=diag){
					trils[i][j] =0;
				}else{
					trils[i][j] = ds[j];
				}
			}
		}
		//Utils.writeMatToFile(trils, true);
		return trils;
	}
	

}
