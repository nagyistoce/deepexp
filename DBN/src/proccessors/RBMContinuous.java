package proccessors;
/**
 * Contains the truncated exponentials to sample data
 */
import utils.MathUtils;
import utils.Utils;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;
/**
 * Two things I need to make sure architecturally for these machines. Make them as purely processing units.
 * 1) Not to have side effects on the data sent into the system.
 * 2)    
 * @author sid
 *
 */
public class RBMContinuous implements RBM {
	
	// RBM SPECIFIC PROPERTIES#
	int numVis, numHid;
	int debug=0;
	// weights initializad as random units
	RealMatrix W = null;//new double[numVis][numHid];
	//biases of visible and hidden layer
	RealVector visBias = null;// srandom.randn(numVis,)*0.1 # 
	RealVector hidBias = null;
	RealMatrix dw = null;
	RealVector dvisBias = null;
	RealVector dhidBias = null;
	// Algorithm specific properties. Will change when algo runs
	RealMatrix data = null;
	RealMatrix hidDistribution = null;
	RealMatrix hidSample= null;
	RealMatrix hidEnergy = null;
	
	RealMatrix visDistribution = null;
	RealMatrix visSample = null;
	RealMatrix visEnergy = null;
	
	
	public RBMContinuous(int numVis, int numHid, int N, RealMatrix WStart){
		this.numHid = numHid;
		this.numVis = numVis;
		
		RealVector[] t = MathUtils.random2DVector(numVis, numHid);
		if(WStart == null){
			W = new Array2DRowRealMatrix(numVis, numHid);//new double[numVis][numHid];
			int row = 0; 
			for (RealVector arrayRealVector : t) {
				W.setRowVector(row, arrayRealVector);
				row++;
			}
		}else
			this.W = WStart;
		//biases of visible and hidden layer
		visBias = new ArrayRealVector(numVis);
		hidBias = new ArrayRealVector(numHid);
		dw = new Array2DRowRealMatrix(numVis,numHid);
		dvisBias = MathUtils.zeros1D(numVis);
		dhidBias = MathUtils.zeros1D(numHid);
		
		data = null;
		
		this.hidSample = new Array2DRowRealMatrix(N, this.numHid);
		this.hidDistribution =  new Array2DRowRealMatrix(N, this.numHid);
		
		this.visDistribution = new Array2DRowRealMatrix(N,this.numVis);
		this.visSample = new Array2DRowRealMatrix(N, this.numVis);
	
	}
	/**
	 * P(v|h)
	 * @param initialData If the data is presented for the first time it is used to form the hidden distribution....
	 * @return
	 */
	public RealMatrix evalHiddenDistribution(RealMatrix initialData){
		
		if(initialData != null){
			this.visSample = initialData;// Also Make sure that initial Data  is not changed due to visSample...
		}
			
		if(this.debug == 1) 
			System.out.println("here 1");
		
		// Discussion with bengio reveals that we implement the signs correctly in the energy and the exponentials...
	    // as per tech report here http://books.nips.cc/papers/files/nips19/NIPS2006_0739.pdf
		
		this.hidEnergy  = Utils.addVect(this.visSample.multiply(this.W), this.hidBias);
		//hidDistribution = 1/(1-exp(-hidEnergy)) - 1/hidEnergy;
		this.hidDistribution = calcExpectationEBE(this.hidEnergy);
		
		return this.hidDistribution;
		
	}
	/**Sample drawn from P(h|v)
	 * Sample the hidden layer.
	 * @return the hidden Layer sample
	 */
	public RealMatrix sampleHiddenLayer(){
		int N= this.hidDistribution.getRowDimension();
			
		//Uniformly distributed sampling points
		Array2DRowRealMatrix s = MathUtils.random2DMatrix(N, this.numHid);
		// We carry out : log(1-s*(1-exp(this.hidEnergy)))/this.hidEnergy;
		
		this.hidSample = sampleFactorialDist(this.hidEnergy);
		return this.hidSample;
		
		
	}
	/**
	 * Evaluate the distribution P(v|h)
	 * @return the distribution. Distribution is factorial.
	 */
	
	public RealMatrix evalVisDistribution(){
		
		this.visEnergy  = Utils.addVect(this.hidSample.multiply(this.W.transpose()), this.visBias );
		//this.visDistribution = 1/(1-exp(-this.visEnergy)) - 1/this.visEnergy;	
		this.visDistribution = calcExpectationEBE(this.visEnergy);
	
		return this.visDistribution;
					
	}
	
	
	//sampling needs to be done from the distribution calculated
	//by above method.
	public RealMatrix sampleVisLayer(){
		int N= this.visDistribution.getRowDimension();
		//Uniformly distributed sampling points
		Array2DRowRealMatrix s = MathUtils.random2DMatrix(N, this.numVis);
		// We carry out : log(1-s*(1-exp(this.visEnergy)))/this.visEnergy;
		this.visSample =  sampleFactorialDist(this.visEnergy);;
		
		return this.visSample;
	
	}
	/**
	 * Calculate the expectation or the distribution
	 */
	 public RealMatrix calcExpectationEBE(RealMatrix m){
		if(m!=null){
			RealMatrix mNew = new Array2DRowRealMatrix(m.getRowDimension(), m.getColumnDimension());
			for(int i =0 ; i<m.getRowDimension(); i++){
				for(int j =0;j<m.getColumnDimension();j++){
					double e = m.getEntry(i, j);
					if(Math.abs(e)<0.01)
						mNew.setEntry(i, j, 0.5 + e*(1./12. - Math.pow(e, 2)/720.));
					else{
						mNew.setEntry(i, j, 1/(1-Math.exp(-e)) - 1/e);
					}
				}
			}
			return mNew;
		}else{
			System.err.println("Null matrix detected. Possible error");
			return null;
		}
	}
	 /**
	  * In this implementation the sampling is done by 
	  * */
	public RealMatrix sampleFactorialDist(RealMatrix m){
		//log(1-s*(1-exp(m)))/m;
		if(m!=null){
			RealMatrix mNew = new Array2DRowRealMatrix(m.getRowDimension(), m.getColumnDimension());
			for(int i =0 ; i<m.getRowDimension(); i++){
				for(int j =0;j<m.getColumnDimension();j++){
					double e = m.getEntry(i, j);
					double s = Math.random();
					 if( Math.abs(e) <= 1e-5 ){
				            mNew.setEntry(i, j, s + ( s*(1 - s)/2 )*e);
					 } else{
						    mNew.setEntry(i, j, Math.log( 1-s + s*Math.exp(e))/e);
					 }
				}
			}
			return mNew;
		}else{
			System.err.println("Null matrix detected. Possible error");
			return null;
		}
	}
	/**
	 *  
	 * @param numGibbsSteps
	 */
	
	public void cdSteps(int numGibbsSteps){
		// do the Gibbs steps
		for(int i=0; i<numGibbsSteps;i++){
			this.evalHiddenDistribution(null);
			this.sampleHiddenLayer();
			this.evalVisDistribution();
			this.sampleVisLayer();
			this.evalHiddenDistribution(null);
			this.sampleHiddenLayer();
		}
	}
	/**
	 * Train the underlying machine using CD.
	 * @param ho
	 * @param cdSteps
	 * @param decay
	 * @param momentum
	 * @param epsilon
	 */
	
	public void  trainRBM(RealMatrix initialData, int cdSteps, double decay, double momentum, double epsilon){
		
		RealMatrix ho = this.evalHiddenDistribution(initialData);
		this.data = initialData;
		RealMatrix VoHo=  initialData.transpose().multiply(ho);
		cdSteps(cdSteps);//Do the CD steps..
		
		
		RealMatrix Vinf = this.visSample ;
		RealMatrix Hinf = this.hidDistribution;
		
		RealMatrix VinfHinf = Vinf.transpose().multiply(Hinf);
		
		//if isnan(this.hidDistribution).any() or isnan(this.hidSample).any() or isnan(this.visDistribution).any() or isnan(this.visSample).any():
		int N = this.data.getRowDimension();
		//update the weights
		this.dw = this.dw.scalarMultiply(momentum).add(
				VoHo.subtract(VinfHinf).scalarMultiply(1/(double)N).
					subtract(this.W.scalarMultiply(decay)).scalarMultiply(epsilon));// regularizer based penalisation
	
		this.W = this.W.add(this.dw);
		
		//update visible layer biases
		RealVector dataTerm = Utils.sum(this.data, 1);//this.data.sum(axis=0)
		RealVector modelTerm = Utils.sum(this.visSample, 1);//this.visSample.sum(axis=0)
		
		this.dvisBias = this.dvisBias.mapMultiplyToSelf(epsilon).add(
				dataTerm.subtract(modelTerm).mapDivideToSelf(N).mapMultiplyToSelf(epsilon)); 
		
		this.visBias= this.visBias.add(this.dvisBias);
		
		
		//update hidden biases
		dataTerm = Utils.sum(ho,1);
		modelTerm = Utils.sum(this.hidSample,1);
		this.dhidBias = this.dhidBias.mapMultiplyToSelf(epsilon).add(
				dataTerm.subtract(modelTerm).mapDivideToSelf(N).mapMultiplyToSelf(epsilon)); 
			  
		
		this.hidBias = this.hidBias.add(this.dhidBias);
	}
	public void updateRBM(){
		
	}
}
