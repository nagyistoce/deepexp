package readers;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;

import analysis.DisplayImage;

import utils.Utils;

import datastructures.Pair;

public class FaceDataReader {
	private static BufferedReader br = null;
	private static BufferedReader br2 = null;
	
	public FaceDataReader(String file){
		try {
			br = new BufferedReader(new FileReader(file+"/faces/faces.txt"));
			br2 = new BufferedReader(new FileReader(file+"/faces/backgrounds.txt"));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
	
	public Pair<RealMatrix, RealMatrix> generateSplitData(Pair<Integer,Integer> ratio) {
		//Split the data into some ratio.
		
		//trainingData
		String face;
		try {
			face = br.readLine();
			String bg = br2.readLine();
			Double [] data = null;
			List<Double[]> faces = new ArrayList<Double[]>();
			while(face!=null){
				String [] facearr = face.split("\t");
				data = new Double[facearr.length];
				int i=0;
				for (String pixel : facearr) {
					data[i++]=Double.valueOf(pixel);
				}
				faces.add(data);
				face= br.readLine();
			}
			data = null;
			List<Double[]> bgs = new ArrayList<Double[]>();
			
			while(bg!=null){
				String [] bgarray = bg.split("\t");
				data = new Double[bgarray.length];
				int i=0;
				for (String pixel : bgarray) {
					data[i++]=Double.valueOf(pixel);
				}
				bgs.add(data);
				bg= br2.readLine();
			}
			RealMatrix allData = new Array2DRowRealMatrix(faces.size()+bgs.size(),data.length);
			RealMatrix trainingData = new Array2DRowRealMatrix((ratio.getData1()*(faces.size()+bgs.size()))/10,data.length);
			RealMatrix testingData = new Array2DRowRealMatrix((ratio.getData2()*(faces.size()+bgs.size()))/10,data.length);
			//we have to split the data by the ratio
			int i =0;
			int trainingSize1 = ratio.getData1()*faces.size()/10;// of faces for training. 
			int trainingSize2 = ratio.getData1()*bgs.size()/10;// of bg for training. 
			for (Iterator<Double[]> iterator = faces.iterator(); iterator.hasNext();) {
				allData.setRow(i++, ArrayUtils.toPrimitive(iterator.next()));
			}
			for (Iterator<Double[]> iterator = bgs.iterator(); iterator.hasNext();) {
				allData.setRow(i++, ArrayUtils.toPrimitive(iterator.next()));
			}
			//set the training faces
			trainingData.setSubMatrix(allData.getSubMatrix(0, trainingSize1, 0, allData.getColumnDimension()-1).getData()
					, 0, 0);
			//set the bg faces
			trainingData.setSubMatrix(allData.getSubMatrix(faces.size(), faces.size()+trainingSize2, 0, allData.getColumnDimension()-1).getData()
					, trainingSize1-1, 0);
			//Utils.writeMatToFile(trainingData.getData(), false);
			testingData.setSubMatrix(allData.getSubMatrix(trainingSize1,faces.size()-1,0,allData.getColumnDimension()-1).getData(),
					0, 0);
			testingData.setSubMatrix(allData.getSubMatrix(faces.size()+trainingSize2,faces.size()+bgs.size()-1,0,allData.getColumnDimension()-1).getData(),
					faces.size()-trainingSize1, 0);
			DisplayImage.renderImageFromArray(trainingData.getRowVector(0), 32, 32, true);
			DisplayImage.renderImageFromArray(trainingData.getRowVector(trainingSize1-2), 32, 32, false);
			DisplayImage.renderImageFromArray(trainingData.getRowVector(trainingSize1-1), 32, 32, false);
			DisplayImage.renderImageFromArray(trainingData.getRowVector(trainingSize1), 32, 32, false);
			DisplayImage.renderImageFromArray(trainingData.getRowVector(trainingSize1+1), 32, 32, false);
			
			DisplayImage.renderImageFromArray(testingData.getRowVector(0), 32, 32, false);
			DisplayImage.renderImageFromArray(testingData.getRowVector(79), 32, 32, false);
			DisplayImage.renderImageFromArray(testingData.getRowVector(80), 32, 32, false);
			DisplayImage.renderImageFromArray(testingData.getRowVector(81), 32, 32, false);
			DisplayImage.renderImageFromArray(testingData.getRowVector(200), 32, 32, false);
			
			br.close();
			br2.close();
		
		return new Pair<RealMatrix,RealMatrix>(trainingData, testingData);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}
}
