package testbed;

import org.apache.commons.math.linear.*;

import datastructures.Pair;

import proccessors.*;
import readers.FaceDataReader;
import utils.Utils;

public class RBMTest  {
	
	RBM rm = null;
	
	public static void main(String[] args){
		RBMTest tst = new RBMTest();
		FaceDataReader f = new FaceDataReader("/home/sid/datasets/");
		f.generateSplitData(new Pair<Integer,Integer>(8,2));
		//tst.testRBM();
	}
	public void testRBM(){
		RealMatrix m;
		/*double[][] W = {{1,0,1,0},{0,1,0,1}};
		rm = new RBMBinomial(2, 4, 4,new Array2DRowRealMatrix(W));
		double[][] a = {{0,0},{1,0},{0,1},{1,1}};
		m = new Array2DRowRealMatrix(a);
		rm.evalHiddenDistribution(m);
		
		for(int i=0;i<0;i++){
			rm.evalHiddenDistribution(null);
			
			rm.sampleHiddenLayer();
			//System.out.println(m);
			rm.evalVisDistribution();
			//System.out.println(m);
			m = rm.sampleVisLayer();
			System.out.println(m);
			//System.out.println(m);
			System.out.println("========");
		
		}*/
		//generate data
		double[][] data = new double[1000][4];
		for(int n =0;n<1000;n++){
	        double r = Math.random();
	        double[] p1 =  {1,0,1,0};
	        double[] p2 =  {0,1,0,1};
	        if(r>0.666){
	        	data[n] =p1;
	        }
	        else if(r>0.333) {
	        	data[n] = p2;
	        }
		}
		m = new Array2DRowRealMatrix(data);
		rm = new RBMContinuous(4, 4, 1000, null);
		double err = Double.MAX_VALUE;
		long l1 = System.currentTimeMillis();
		for(int i=0;i<1500 ;i++){
			((RBMContinuous)rm).trainRBM(m, 1, 0.001, 0, 0.3);
			if(i%100==0){
				RealVector v  = Utils.sum(rm.sampleVisLayer(), 1);
				RealVector v1 = Utils.sum(m, 1);
				err = v.subtract(v1).mapDivide(1000.).getL1Norm();
				System.out.println("Err: "+err);
				
				System.out.println("Sample from RBM: "+rm.sampleVisLayer().getRowVector(0));
				System.out.println("Sample from data: "+m.getRowVector(0));
			}
		}
		System.out.println("Time taken :"+(System.currentTimeMillis()-l1));
		

	}
/*
 * def test_rbm_sample_h():
    I, J = 2, 4

    bm = rbm.RBM(I,J)
    bm.w[0,:] = [1,0,1,0]
    bm.w[1,:] = [0,1,0,1]
    bm.w *= 2e4
    bm.bv *= 0
    bm.bh *= 0

    # test 1
    v = scipy.array([[0,0],[1,0],[0,1],[1,1.]])
    h = []
    for n in range(100):
        prob, sample = bm.sample_h(v)
        h.append(sample)

    # check inferred probabilities
    expected_probs = scipy.array([[0.5, 0.5, 0.5, 0.5],
                                  [1.0, 0.5, 1.0, 0.5],
                                  [0.5, 1.0, 0.5, 1.0],
                                  [1.0, 1.0, 1.0, 1.0]])
    assert_array_almost_equal(prob, expected_probs, 8)

    # check sampled units
    h = scipy.array(h)
    for n in range(4):
        distr = h[:,n,:].mean(axis=0)
        assert_array_almost_equal(distr, expected_probs[n,:], 1)

    # test 2, with bias
    bm.bh -= 1e4
    h = []
    for n in range(100):
        prob, sample = bm.sample_h(v)
        h.append(sample)

    # check inferred probabilities
    expected_probs = scipy.array([[0., 0., 0., 0.],
                                  [1.0, 0., 1.0, 0.],
                                  [0., 1.0, 0., 1.0],
                                  [1.0, 1.0, 1.0, 1.0]])
    assert_array_almost_equal(prob, expected_probs, 8)

    # check sampled units
    h = scipy.array(h)
    for n in range(4):
        distr = h[:,n,:].mean(axis=0)
        assert_array_almost_equal(distr, expected_probs[n,:], 1)
		 * */
}
