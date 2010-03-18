package datastructures;

public class Triple<T,W,V> {
	T data1;
	W data2;
	V data3;
	public Triple(T data1, W data2, V data3){
		this.data1 = data1;
		this.data2 = data2;
		this.data3 = data3;
	}
	public T getData1(){
		return data1;
	}
	public W getData2(){
		return data2;
	}
	public V getData3(){
		return data3;
	}

}
