package utils;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.ArrayRealVector;

public class MathUtils {
	
	public static ArrayRealVector[] zeros2D(int d1, int d2){
		ArrayRealVector[] W = new ArrayRealVector[d1];
		for(int i = 0; i<d1; i++){
			W[i]=new ArrayRealVector(d2);
			W[i].set(0);
		}
		return W;
	}
	
	public static ArrayRealVector zeros1D(int size){
		ArrayRealVector W = new ArrayRealVector(size);
		W.set(0);
		return W;
	}
	public static ArrayRealVector[] random2DVector(int d1, int d2){
		ArrayRealVector[] W = new ArrayRealVector[d1];
		for(int i = 0; i<d1; i++){
			double[] temp = new double[d2];
			for(int j = 0 ;j<d2; j++){
				temp[j]= Math.random();
			}
			W[i] = new ArrayRealVector(temp); 
		}
		return W;
	}
	public static Array2DRowRealMatrix random2DMatrix(int d1, int d2){
		Array2DRowRealMatrix W = new Array2DRowRealMatrix(d1, d2);
		for(int i = 0; i<d1; i++){
			for(int j = 0 ;j<d2; j++){
				W.setEntry(i, j, Math.random());
			}
		}
		return W;
	}
	public static ArrayRealVector random1D(int size){
		ArrayRealVector W = new ArrayRealVector(size);
		for(int i = 0; i<size; i++)
			W.setEntry(i, Math.random());
		return W;
	}
	/**
	 * Add the logarithm 
	 */
	public static double logadd(double l1, double l2){
		return 0;
		
	}
}
