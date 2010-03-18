package proccessors;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.ArrayRealVector;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;

import utils.MathUtils;
import utils.Utils;

public class RBMBinomial implements RBM {
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
	
	
	public RBMBinomial(int numVis, int numHid, int N, RealMatrix WStart){
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

	public RealMatrix evalHiddenDistribution(RealMatrix initialData){
		
		if(initialData != null){
			this.visSample = initialData;// Also Make sure that initial Data  is not changed due to visSample...
		}
			
		if(this.debug == 1) 
			System.out.println("here 1");
		
		// Discussion with bengio reveals that we implement the signs correctly in the energy and the exponentials...
	    // as per tech report here http://books.nips.cc/papers/files/nips19/NIPS2006_0739.pdf
		
		this.hidEnergy  = Utils.addVect(this.visSample.multiply(this.W), this.hidBias);
		// probs = 1./(1. + exp(-self.bh - dot(v, self.w)))
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
		this.hidSample = sampleFactorialDist(this.hidDistribution);
		return this.hidSample;
		
		
	}
	/**
	 * Evaluate the distribution P(v|h)
	 * @return the distribution. Distribution is factorial.
	 */
	
	public RealMatrix evalVisDistribution(){
		
		this.visEnergy  = Utils.addVect(this.hidSample.multiply(this.W.transpose()), this.visBias );
		this.visDistribution = calcExpectationEBE(this.visEnergy);
	
		return this.visDistribution;
					
	}
	
	
	//sampling needs to be done from the distribution calculated
	//by above method.
	public RealMatrix sampleVisLayer(){
		//Uniformly distributed sampling points
		this.visSample =  sampleFactorialDist(this.visDistribution);
		
		return this.visSample;
	
	}
	@Override
	public RealMatrix calcExpectationEBE(RealMatrix m) {
		//probs = 1./(1. + exp(-self.bh - dot(v, self.w)))
		if(m!=null){
			RealMatrix mNew = new Array2DRowRealMatrix(m.getRowDimension(), m.getColumnDimension());
			for(int i =0 ; i<m.getRowDimension(); i++){
				for(int j =0;j<m.getColumnDimension();j++){
					double e = m.getEntry(i, j);
					 //probs = 1./(1. + exp(-self.bh - dot(v, self.w)))
					mNew.setEntry(i, j, 1./(1+Math.exp(e)));
					
				}
			}
			return mNew;
		}else{
			System.err.println("Null matrix detected. Possible error");
			return null;
		}
	}

	@Override
	public RealMatrix sampleFactorialDist(RealMatrix m) {
		//log(1-s*(1-exp(m)))/m;
		if(m!=null){
			RealMatrix mNew = new Array2DRowRealMatrix(m.getRowDimension(), m.getColumnDimension());
			for(int i =0 ; i<m.getRowDimension(); i++){
				for(int j =0;j<m.getColumnDimension();j++){
					double e = m.getEntry(i, j);
					double s = Math.random();
					mNew.setEntry(i, j, e>s?0:1);
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

}
