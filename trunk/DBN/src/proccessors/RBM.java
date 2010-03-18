package proccessors;

import org.apache.commons.math.linear.RealMatrix;

//The design of an RBM can be seen as a data processing block
// It is similar to what MDP's 
public interface RBM {
	public RealMatrix calcExpectationEBE(RealMatrix m);
	public RealMatrix sampleFactorialDist(RealMatrix m);
	public RealMatrix evalHiddenDistribution(RealMatrix m);
	public RealMatrix sampleHiddenLayer();
	public RealMatrix evalVisDistribution();
	public RealMatrix sampleVisLayer();
}
