package datastructures;

public class Pair<T,W> {
	T data1;
	W data2;
	public Pair(T data1, W data2){
		this.data1 = data1;
		this.data2 = data2;
	}
	public T getData1(){
		return data1;
	}
	public W getData2(){
		return data2;
	}
	public  String toString(){
		return data1+"::"+data2;
	}
	
}
