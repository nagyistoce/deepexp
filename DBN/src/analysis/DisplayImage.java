package analysis;

import java.awt.Color;

import javax.swing.JFrame;
import javax.swing.JPanel;

import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.GrayPaintScale;
import org.jfree.chart.renderer.PaintScale;
import org.jfree.chart.renderer.xy.XYBlockRenderer;
import org.jfree.data.DomainOrder;
import org.jfree.data.general.DatasetChangeListener;
import org.jfree.data.general.DatasetGroup;
import org.jfree.data.xy.MatrixSeries;
import org.jfree.data.xy.MatrixSeriesCollection;
import org.jfree.data.xy.XYZDataset;
import org.jfree.ui.RefineryUtilities;

import utils.Utils;

/*
 * This class is going to encapsulate some static methods to display images.
 * It uses JfreeChart under the hood. Use this class to display Images.
 * Feature Req:
 * 1) Display basic image
 * 2) Refresh images
 * 3) Display multiple images in a grid fashion...
 * */

public class DisplayImage {
	
	private static JFrame nj = new JFrame();
	private static JPanel chartPanel =null;
	/**
	 * This method will probably create a XYZ data set out of the RealMatrix
	 * @param m
	 */
	public static void renderImage(RealMatrix m){
		XYZDataset mat = convertToXYZ(m);
		new ChartPanel(createChart(mat,m.getRowDimension(),m.getColumnDimension(), 0.0));
		
		
	}
	/**
	 * 
	 * @param v
	 * @param rowDim how the array will be reshaped...
	 * @param colDim how
	 */
	public static void renderImageFromArray(RealVector v, int rowDim, int colDim, boolean sameFrame){
		v.mapDivideToSelf(v.getNorm());
		if(v.getDimension() != rowDim*colDim)
			throw new IllegalArgumentException("Cant resize "+v.getDimension()+" element(s) to"+rowDim+"*"+colDim);
		XYZDataset mat = convertToXYZ(v, rowDim, colDim);
		//Initialize the main window
		
		//initialize the XYchart 
		//if(chartPanel==null)
		chartPanel = new ChartPanel(createChart(mat,rowDim,colDim,  Utils.max(v).getData2()));
		
	    //What size?
		
		chartPanel.setPreferredSize(new java.awt.Dimension(300, 300));// Sets up the size of the main window.
	    nj.setContentPane(chartPanel);
	    nj.pack();
        RefineryUtilities.centerFrameOnScreen(nj);
        nj.setVisible(true);
			
	}
	/**
	 * 3*5
	 * Convert a
	 * @param m
	 * @return
	 */
	private static XYZDataset convertToXYZ(RealVector m, int rowDim, int colDim){
		MatrixSeries mat = new MatrixSeries("DefaultData",rowDim, colDim);
		MatrixSeriesCollection m1 = new MatrixSeriesCollection(mat);
		for(int i = 0; i < rowDim; i++ ){
			for(int j = 0; j < colDim; j++){
				mat.update(i, j, m.getEntry(j+((colDim)*i)));
			}
		}
		return m1;
	}
	/**
	 * Make a XYZ matrix out of a  RealMatrix... Aha.
	 * 
	 * */
	private static XYZDataset convertToXYZ(RealMatrix m){
		MatrixSeries mat = new MatrixSeries("DefaultData",m.getRowDimension(),m.getColumnDimension());
		MatrixSeriesCollection m1 = new MatrixSeriesCollection(mat);
		for(int i=0;i<m.getRowDimension();i++){
			for(int j=0; j < m.getColumnDimension(); j++){
				mat.update(i, j, m.getEntry(i, j));
			}
		}
		return m1;
	}
	/**
     * Creates a sample chart.
     * 
     * @param dataset  the dataset.
     * 
     * @return A sample chart.
     */
    private static JFreeChart createChart(XYZDataset dataset, int rowSize, int colSize, double scaleVal) {
        NumberAxis xAxis = new NumberAxis("X");
        xAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
        xAxis.setLowerMargin(0.0);
        xAxis.setUpperMargin(100.0);
        xAxis.setRange(0,rowSize);
        NumberAxis yAxis = new NumberAxis("Y");
        yAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
        yAxis.setLowerMargin(0.0);
        yAxis.setUpperMargin(100.0);
        yAxis.setRange(0, colSize);
        
        XYBlockRenderer renderer = new XYBlockRenderer();//A Chart with (x,y) values...
        PaintScale scale = new GrayPaintScale(0.0, scaleVal);// the values of the pixels are scaled using this...?
        renderer.setPaintScale(scale);
        XYPlot plot = new XYPlot(dataset, xAxis, yAxis, renderer);
        plot.setBackgroundPaint(Color.lightGray);
        plot.setDomainGridlinesVisible(false);
        plot.setRangeGridlinePaint(Color.white);
        JFreeChart chart = new JFreeChart("Default Image", plot);
        chart.removeLegend();
        chart.setBackgroundPaint(Color.white);
        return chart;
    }
}
