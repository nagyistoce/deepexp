package utils;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.TreeMap;

import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.MatrixUtils;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;

import datastructures.Pair;

public class Utils {
	//Exact behaviour as find in matlab...
	public static<T extends Number> void writeMatToFile(T[][] data, boolean append){
		try {
			BufferedWriter bfr = new BufferedWriter(new FileWriter("TempFile.txt", append));
			System.out.println("Writing te array to TempFile.txt");
			bfr.write("Writing te array to TempFile.txt\n");
			for (int j = 0; j < data.length; j++) {
				T[] ts = data[j];
				for (int i = 0; i < ts.length; i++) {
					bfr.write(ts[i]+", ");
				}bfr.newLine();
				
				
			}
			bfr.close();
		}
		catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	public static double [][] hadamardProd(double[][]  dat1 , double[][] dat2){
		double [][] datnew = new double[dat1.length][dat1[0].length]; 
		for (int i = 0; i < dat1.length; i++) {
			double[] ds1 = dat1[i];
			double[] ds2 = dat2[i];
			for (int j = 0; j < ds1.length; j++) {
				datnew[i][j] = ds1[j]*ds2[j];
				
			}
			
		}
		return datnew;
	}
	/**
	 * Ovveride support for the hadamard product for the Apache math library data structures
	 * @param rm1 - matrix 1 
	 * @param rm2 - matrix 2
	 * 
	 */
	public static RealMatrix hadamardProd(RealMatrix rm1, RealMatrix rm2){
		RealMatrix rm3 = new Array2DRowRealMatrix(rm1.getRowDimension(), rm1.getColumnDimension());
		for(int i =0;i<rm1.getColumnDimension();i++){
			rm3.setColumnVector(i, rm1.getColumnVector(i).ebeMultiply(rm2.getColumnVector(i)));
		}
		return rm3;
	}
	
	
	public static void writeMatToFile(double[][] data, boolean append){
		try {
			
			BufferedWriter bfr = new BufferedWriter(new FileWriter("TempFile.txt", append));
			System.out.println("Writing te array to TempFile.txt");
			bfr.write("Writing te array to TempFile.txt\n");
			for (int j = 0; j < data.length; j++) {
				double[] ts = data[j];
				for (int i = 0; i < ts.length; i++) {
					bfr.write(ts[i]+", ");
				}bfr.newLine();
				
				
			}
			bfr.close();
		}
		catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	public static<T extends Number> void writeArrToFile(T[] data, boolean append){
		try {
			BufferedWriter bfr = new BufferedWriter(new FileWriter("TempFile.txt", append));
			System.out.println("Writing te array to TempFile.txt");
			bfr.write("Writing te array to TempFile.txt\n");
			for (int i = 0; i < data.length; i++) {
				bfr.write(data[i]+", ");
			}bfr.newLine();
				
				
			bfr.close();
		}
		catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	public static int[] findPos(double[] vector){
		// find the indices of the matrix>0
		int[] data = new int[vector.length];
		int count =0;
		for (int j = 0; j < data.length; j++) {
			if(vector[j]>0){
				data[count] = j;
				count++;
			}
		}
		return (int[])ArrayUtils.subarray(data, 0, count);
		
		
	}
	/**
	 * Matlab clone of find(a==b). Means find the indices of vector a where the 
	 * value b occured in a. T must override equals for this to work in the intuitive sense 
	 * */
	public static<T extends Comparable<T> > int[] find(T[] data, T target){
		ArrayList<Integer> indices = new ArrayList<Integer>();
		for (int i = 0; i < data.length; i++) {
			if(data[i].equals(target))
				indices.add(i);
		}
		int count = 0;
		int[] primInd = new int[indices.size()];
		for (Iterator<Integer> iterator = indices.iterator(); iterator.hasNext();) {
			Integer integer = iterator.next();
			primInd[count] = integer;
			count++;
		}
		return primInd;
		
	}
	/**
	 * Matlab replica of unique. Gives you unique numbers(Set) interpretation
	 * of a bag of numbers some or all of which may be repeated.
	 * */
	public static<T extends Number> T[] unique(T[] component, T[] data){
		//T[] component1 = new T[component.length];
		//Was unable to create a generic array and return it as such.
		if(data==null) throw new NullPointerException();
		
		HashSet<T> set = new HashSet<T>();
		for (T t : component) {
			set.add(t);
		}
		int count = 0;
		for (Iterator<T> iterator = set.iterator(); iterator.hasNext();) {
			T t =  iterator.next();
			data[count]= t;
			count++;
		}
		//data reference is changed... Why is the cast safe?
		data = (T[])ArrayUtils.subarray(data, 0, count);//This would work 
		return data;
	}
	/**
	 * A Useful utility to find the submatrix of a Number
	 * Can we send a array of objects of type T back?
	 * The cast to T[] will never work. Right now I can only think 
	 * of sending
	 * */
	public static <T extends Number> T[] subMatrix(T data[], int indices[], T[] subMatrix){
		
		if(subMatrix==null || data ==null) throw new NullPointerException();
		//Array.n
		//Class<?> c = Class.forName(data[0].getClass());
        //Object o = Array.newInstance(c, n);

		ArrayList<T> sub = new ArrayList<T>();
		for (int i = 0; i < indices.length; i++) {
			sub.add(data[indices[i]]);
		}
		int count = 0;
		for (Iterator iterator = sub.iterator(); iterator.hasNext();) {
			T t = (T) iterator.next();
			subMatrix[count] = t;
			count++;
		}
		subMatrix = (T[])ArrayUtils.subarray(subMatrix, 0, count);//This would work 
		return subMatrix;
		
	}
	/**
	 * Return the first nth statistic statistics of an array and there indices.
	 * The idea is to return the elements and there indices in the sorted order
	 * @param <T>
	 * @return
	 */
	public static<T> TreeMap<T, Integer> nLargest(T[] data, int n){
		TreeMap<T, Integer> tm = new TreeMap<T, Integer>();
		for(int i = 0;i<data.length;i++)
			tm.put(data[i], i);
		int size = data.length;
		
		T[] st= (T[])tm.keySet().toArray();
		
		for (int i=0;i<size-n;i++) {
			tm.remove(st[i]);
		}
		return tm;
		
	}
	
	
	
	///////////////////////////////////////////////////MAX MIN FUNCTIONS///////////////////////////////////////////////
	/**
	 * Return a the maximum Index of the array
	 * @param <T>
	 * @param array
	 * @return
	 */
	
	public static<T extends Comparable<T> > int maxInd(T[] array){
		T max = array[0];
		int maxInd=0;
		for (int i = 0; i < array.length; i++) {
			if(array[i].compareTo(max)>0){
				max = array[i];
				maxInd = i;
			}
		}
		return maxInd;
	}
	/**
	 * Return the max value and its index. If the element is repeated then the first occurance is returned
	 * @param rv
	 * @return
	 */
	public static Pair<Integer,Double> max(RealVector rv){
		Double maxOne = new Double(Double.MIN_VALUE);
		Integer maxInd = new Integer(-1);
		for(int i = 0; i< rv.getDimension() ; i++ ){
			if(maxOne<rv.getEntry(i)){
					maxOne = rv.getEntry(i);
					maxInd = i; 
			}
			
		}
		return new Pair<Integer,Double>(maxInd,maxOne);
	}
	/**
	 * The function will find the min value and its index and return in
	 * @param RealVector- The vector whose min needs to be found 
	 * @return The min value and its index
	 */
	public static Pair<Integer,Double> min(RealVector rv){
		Double minOne = new Double(Double.MAX_VALUE);
		Integer minInd = new Integer(-1);
		for(int i = 0; i< rv.getDimension() ; i++ ){
			if(minOne > rv.getEntry(i)){
					minOne = rv.getEntry(i);
					minInd = i; 
			}
			
		}
		return new Pair<Integer,Double>(minInd,minOne);
	}
	/**
	 * The idea of this function is to return the max value along a specific dimension.
	 * If X = [[2, 8, 4];   
     *        [7, 3, 9]], then max(X,0) is [7 8 9],
     * So it is max of each column.        
     * max(X,1) gives [8,9]       
	 * @param mtr
	 * @param rowOrCol
	 * @return The function returns the indices as well as the max values
	 */
	public static Pair<RealVector,RealVector> max(RealMatrix mtr, int rowOrCol){
		RealVector maxVec = new ArrayRealVector(mtr.getColumnDimension());
		RealVector maxInd = new ArrayRealVector(mtr.getColumnDimension());
		if(rowOrCol==0){
			//find max of all columns
			for(int i = 0; i < mtr.getColumnDimension();i++){
				 RealVector rv = mtr.getColumnVector(i);
				 Pair<Integer,Double> p = max(rv);
				 maxVec.setEntry(i, p.getData2());
				 maxInd.setEntry(i, p.getData1());
			}
		}else{
			for(int i = 0; i < mtr.getRowDimension();i++){
				 RealVector cv = mtr.getColumnVector(i);
				 Pair<Integer,Double> p = max(cv);
				 maxVec.setEntry(i, p.getData2());
				 maxInd.setEntry(i, p.getData1());
			}
		}
		return new Pair<RealVector,RealVector>(maxInd,maxVec);
		
	}
	/**
	 * The idea of this function is to return the min value along a specific dimension.
	 * If X = [[2, 8, 4];   
     *        [7, 3, 9]], then max(X,0) gives [2 3 4] as values and [1,2,1] as indices
     * So it is max of each column.        
     * max(X,1) gives [8,9] and[2,3]      
	 * @param mtr
	 * @param rowOrCol
	 * @return The function returns the indices as well as the min values
	 */
	public static Pair<RealVector,RealVector> min(RealMatrix mtr, int rowOrCol){
		RealVector minVec = new ArrayRealVector(mtr.getColumnDimension());
		RealVector minInd = new ArrayRealVector(mtr.getColumnDimension());
		if(rowOrCol==0){
			//find max of all columns
			for(int i = 0; i < mtr.getColumnDimension();i++){
				 RealVector rv = mtr.getColumnVector(i);
				 Pair<Integer,Double> p = min(rv);
				 minVec.setEntry(i, p.getData2());
				 minInd.setEntry(i, p.getData1());
			}
		}else{
			for(int i = 0; i < mtr.getRowDimension();i++){
				 RealVector cv = mtr.getColumnVector(i);
				 Pair<Integer,Double> p = min(cv);
				 minVec.setEntry(i, p.getData2());
				 minInd.setEntry(i, p.getData1());
			}
		}
		return new Pair<RealVector,RealVector>(minVec,minInd);
		
	}
	/* 
	 * Does not work....
	*/
    public static<T> T[] toArray(Collection<?> col){
    	Object[] dat = new Object[col.size()];
    	int count =0;
    	for (Iterator<T> iterator = (Iterator<T>) col.iterator(); iterator.hasNext();) {
			T object =  iterator.next();
			dat[count] = (T)object;
			count++;
		}
    	return (T[])dat;
    }
	//1 to n for any data type double float int....
	public static<T> T[] oneToN(T[] component){
		//T[] component = new T[size];
		for (int i = 0; i < component.length; i++) {
			component[i] =(T)((Integer)(i+1));
		}
		return component;
	}
	/**
	 * Make a N element vector with zero to N-1 Natural numbers 
	 * @param <T> T is the generic type which extends Number.
	 * @param component is the pre sized array
	 * @return
	 */
	public static<T extends Number> T[] zeroToN(T[] component){
		//T[] component = new T[size];
		for (int i = 0; i < component.length; i++) {
			component[i] =(T)((Integer)(i));
		}
		return component;
	}
	/**
	 * returns a vector representing the row or the columns, Same semantics as matlab find
	 * @param mtr The matrix to be summed
	 * @param rowOrCol if 1 then all the rows are added up else all the columns are added up. 
	 * @return the summation along the dimension.
	 */
	public static ArrayRealVector sum(RealMatrix mtr, int rowOrCol ){
		if(rowOrCol==1){
			//sum of all rows into 1 vector
			int numRows = mtr.getRowDimension();
			ArrayRealVector vect = new ArrayRealVector(mtr.getColumnDimension());	
			while(numRows > 0){
				numRows--;
				vect = (ArrayRealVector)vect.add(mtr.getRow(numRows));
			}
			return vect;
			
		}else{
			//sum of all columns into a vector
			int numCols = mtr.getColumnDimension();
			ArrayRealVector vect = new ArrayRealVector(mtr.getRowDimension());	
			while(numCols > 0){
				numCols--;
				vect = (ArrayRealVector)vect.add(mtr.getColumn(numCols));
			}
			return vect;
		}
	}
	
	/**
	 * This function helps in replicating rows*cols tilings of a vector v 
	 * */
	public static RealMatrix repmat(RealMatrix v, int rows, int cols){
		int tcols = v.getColumnDimension()*cols;
		int trows = v.getRowDimension()*rows;
		RealMatrix rm = MatrixUtils.createRealMatrix(trows,tcols);
		for(int i = 0;i<rows;i++){
			for(int j =0;j<cols;j++){
				rm.setSubMatrix(v.getData(), i*v.getRowDimension(), j*v.getColumnDimension());
			}
			
		}
		return rm;
	}
	
	public static double[] sum(double[][] mtr, int rowOrCol ){
		return sum(new Array2DRowRealMatrix(mtr),  rowOrCol ).getData();
	}
	
	public static RealMatrix mapToAbs(RealMatrix m){
		
		RealMatrix m1 = new Array2DRowRealMatrix(m.getRowDimension(), m.getColumnDimension());
		
		for(int i = 0;i<m.getRowDimension();i++){
			for(int j = 0 ; j<m.getColumnDimension();j++){
				if(m.getEntry(i, j)<0)
					m1.setEntry(i, j, -m.getEntry(i, j));
				else
					m1.setEntry(i, j, m.getEntry(i,j));
			}
		}
		return m1;
		
	}
	/**
	 * The matrix vv will be normalized. 
	 * @param vv - The matrix to e normalized
	 * @param rowsOrColumns if this is 0 then Ultimately we end up normalizing each row by dividing it by its 1st norm. 
	 * Implementation: The columns of the N*k matrix(vv^2) will be  added to make a 
	 * N*1 vector. Which will then divide the k Columns. 
	 * 
	 * If this is 1 then k*1 vector will divide the N rows.  Ultimately we end up normalizing each 
	 * column by its first norm. 
	 * 
	 * @return the new normalized matrix
	 * 
	 * Additionally we need to take care of the fact that any row is completely zero... 
	 * We need to leave it untouched that is keep it as such...
	 * Reemember a zero sum is possible if and only if all the elements are zero...
	 */
	public static RealMatrix normalize(RealMatrix vv, int rowsOrColumns){
		RealMatrix mm = new Array2DRowRealMatrix(vv.getRowDimension(),vv.getColumnDimension());
		if(rowsOrColumns==0){
			//sum up the columns and N*K matrix becomes a N*1 vector. 
		    int p =  vv.getColumnDimension();
		    RealVector den =  Utils.sum(Utils.hadamardProd(vv, vv), 0).mapPowToSelf(0.5);	
		   // Utils.writeMatToFile(Utils.hadamardProd(vv, vv).getData(), false);
		   // Utils.writeArrToFile(ArrayUtils.toObject(den.getData()), false);
		    Utils.writeMatToFile(Utils.hadamardProd(vv, vv).getData(), false);
		    for (int i = 0; i < den.getDimension(); i++) {
				if(den.getEntry(i)==0.0)
					den.setEntry(i, 1.0);
				// the zero sums will be if all elements of vv^2 were zero. In that case dividing it by 1 wont 
				//make a diff
			}
		    //int[] zeroInd = Utils.find((Double[])ArrayUtils.toObject(den.toArray()), new Double(0));
		    
		    //Now divide the column by the vector
		    for(int i=0;i<p;i++){
		    	mm.setColumnVector(i, vv.getColumnVector(i).ebeDivide(den));
		    }
		    
		    
		}else{
			int p =  vv.getRowDimension();
		    
		    RealVector den =  Utils.sum(Utils.hadamardProd(vv, vv), 1).mapPowToSelf(0.5);
		    //Now divide the column by the vector
		    for (int i = 0; i < den.getDimension(); i++) {
				if(den.getEntry(i)==0.0)
					den.setEntry(i, 1.0);
				// the zero sums will be if all elements of vv^2 were zero. In that case dividing it by 1 wont 
				//make a diff
			}
		    for(int i=0;i<p;i++){
		    	vv.setRowVector(i, vv.getRowVector(i).ebeDivide(den));
		    }
		}
		return mm;
	}
	/**
	 * Exact matlab replica of [x,y] = find(a>0). Where a is a matrix and [x,y]
	 * will be a pair of indices. 
	 * **/ 
	public static ArrayList<Pair<Integer,Integer>> findPos(double[][] matrix){
		ArrayList<Pair<Integer,Integer>> data = new ArrayList<Pair<Integer,Integer>>();
		
		for (int i = 0; i < matrix.length; i++){
			for (int j = 0; j < matrix[i].length; j++){
				if(matrix[i][j]>0){
					data.add(new Pair<Integer,Integer>(i,j));
				}
			}
		}
		return data;
	}
	//Set the specific indices to a value...
	public static<T extends Number> void set(T[] data,int[] indices, T val){
		for (int i = 0; i < indices.length; i++) {
			if(indices[i]>= data.length)
				throw new ArrayIndexOutOfBoundsException();
			data[indices[i]] = val;
		}
	}
	/**
	 * Return a 0 matrix of size i* ncomp 
	 * @param i
	 * @param ncomp
	 * @return
	 */
	public static double[][] zeros(int i, int ncomp) {
		//new Array2D
		double[][] arr = new double[i][ncomp];
		for (int j2 = 0; j2 < arr.length; j2++) {
			double[] js = arr[j2];
			for (int j = 0; j < js.length; j++) {
				js[j] = 0;
			}
		}return arr;
		
	}
	/**
	 * Invert each element of the matrix 1/Aij . This will change the matrix. 
	 * @param m
	 */
	public static void invertEBESelf(RealMatrix m){
		if(m!=null){
			for(int i =0 ; i<m.getRowDimension(); i++){
				for(int j =0;j<m.getColumnDimension();j++){
					m.setEntry(i, j, 1./m.getEntry(i, j));
				}
			}
			
		}else{
			System.err.println("Null matrix detected. Possible error");
			
		}
		
	}
	/**
	 * Invert each element and return a new matrix
	 * @param m
	 * @return
	 */
	public static RealMatrix invertEBE(RealMatrix m){
		if(m!=null){
			RealMatrix mNew = new Array2DRowRealMatrix(m.getRowDimension(), m.getColumnDimension());
			for(int i =0 ; i<m.getRowDimension(); i++){
				mNew.setRowVector(i, m.getRowVector(i).mapInv());
			}
			return mNew;
		}else{
			System.err.println("Null matrix detected. Possible error");
			return null;
		}
		
	}
	/**
	 * Map exp(Aij) for each element Aij of vector A. Mainly to save space.
	 * @param m matrix to be 
	 */
	public static void exponentEBESelf(RealMatrix m){
		if(m!=null){
			
			for(int i =0 ; i<m.getRowDimension(); i++){
				for(int j =0;j<m.getColumnDimension();j++){
					m.setEntry(i, j,Math.exp(m.getEntry(i, j)));
				}
			}
		}else{
			System.err.println("Null matrix detected. Possible error");
		}
	}
	/**
	 * Map exp(Aij) for each element Aij of vector A. Return a new matrix
	 * @param m matrix to be 
	 */
	public static RealMatrix exponentEBE(RealMatrix m){
		if(m!=null){
			RealMatrix mNew = new Array2DRowRealMatrix(m.getRowDimension(), m.getColumnDimension());
			for(int i =0 ; i<m.getRowDimension(); i++){
				
				mNew.setRowVector(i, m.getRowVector(i).mapExp());
			}
			return mNew;
		}else{
			System.err.println("Null matrix detected. Possible error");
			return null;
		}
	}
	/**
	 * Map (Aij+d) for each element Aij of vector A. Store result back in A. Mainly to save space
	 * @param m matrix to be 
	 */
	public static void addEBESelf(RealMatrix m, double d){
		if(m!=null){
			for(int i =0 ; i<m.getRowDimension(); i++){
				for(int j =0;j<m.getColumnDimension();j++){
					m.setEntry(i, j,m.getEntry(i, j)+d);
				}
			}
		}else{
			System.err.println("Null matrix detected. Possible error");
		}
	}
	/**
	 * Map (A[i]+r) for each row A[i] of vector A. Return a new matrix.
	 * @param m matrix to which the  
	 * @param r vector to be added
	 */
	public static RealMatrix addVect(RealMatrix m, RealVector r){
		if(m!=null){
			RealMatrix mNew = new Array2DRowRealMatrix(m.getRowDimension(), m.getColumnDimension());
			for(int i =0 ; i<m.getRowDimension(); i++){
				mNew.setRowVector(i, m.getRowVector(i).add(r));
			}
			return mNew;
		}else{
			System.err.println("Null matrix detected. Possible error");
			return null;
		}
	}
	/**
	 * Map Aij*d for each element Aij of vector A to self. Return a new matrix. Mainly to save space
	 * @param m matrix to be 
	 */
	public static void multiplyEBESelf(RealMatrix m, double d){
		if(m!=null){
			for(int i =0 ; i<m.getRowDimension(); i++){
				for(int j =0;j<m.getColumnDimension();j++){
					m.setEntry(i, j, m.getEntry(i, j)*d);
				}
			}
		}else{
			System.err.println("Null matrix detected. Possible error");
		}
	}
	/**
	 * Map log(Aij) for each element Aij of vector A. Return a new matrix. Mainly to save space
	 * @param m matrix to be 
	 */
	public static void logEBESelf(RealMatrix m){
		if(m!=null){
			for(int i =0 ; i<m.getRowDimension(); i++){
				for(int j =0;j<m.getColumnDimension();j++){
					m.setEntry(i, j, Math.log(m.getEntry(i, j)));
				}
			}
		}else{
			System.err.println("Null matrix detected. Possible error");
		}
	}
	/**
	 * Add a= a+m. Mainly for space saving. does not replicate  
	 * @param addend 
	 * @param m
	 */
	public static void addMapToSelf(RealMatrix addend, RealMatrix m){
		if(addend!=null && m !=null){
			for(int i =0 ; i<addend.getRowDimension(); i++){
				for(int j =0;j<addend.getColumnDimension();j++){
					addend.setEntry(i, j, addend.getEntry(i, j)+m.getEntry(i, j));
				}
			}
		}else{
			System.err.println("Null matrix detected. Possible error");
		}
	}
	/**
	 * dividend/divisor and the result is stored in dividend. Mailnly space saving  
	 * @param dividend 
	 * @param divisor
	 */
	public static void divideMapToSelf(RealMatrix dividend, RealMatrix divisor){
		if(dividend!=null && divisor !=null){
			for(int i =0 ; i<dividend.getRowDimension(); i++){
				for(int j =0;j<dividend.getColumnDimension();j++){
					dividend.setEntry(i, j, dividend.getEntry(i, j)/divisor.getEntry(i, j));
				}
			}
		}else{
			System.err.println("Null matrix detected. Possible error");
		}
	}
	/**
	 * multiplier*multiplicand and the result is stored in multiplier. Mailnly space saving  
	 * @param multiplier 
	 * @param multiplicand
	 */
	public static void hadamardProdToSelf(RealMatrix multiplier, RealMatrix multiplicand){
		if(multiplier!=null && multiplicand !=null){
			for(int i =0 ; i<multiplicand.getRowDimension(); i++){
				for(int j =0;j<multiplicand.getColumnDimension();j++){
					multiplier.setEntry(i, j, multiplicand.getEntry(i, j)*multiplicand.getEntry(i, j));
				}
			}
		}else{
			System.err.println("Null matrix detected. Possible error");
		}
	}
}
