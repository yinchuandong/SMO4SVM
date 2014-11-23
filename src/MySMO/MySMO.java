package MySMO;

import java.util.ArrayList;

public class MySMO {
	
	/**
	 * 惩罚因子
	 */
	private double C = 0.05;
	
	/**
	 * 松弛变量
	 */
	private double tolerance = 0.001;
	
	/**
	 * 终止条件的差值
	 */
	private double eps = 0.001;
	
	/**
	 * 训练集的已知类别号，即公式中的y[]
	 */
	private int[] target = null;
	
	/**
	 * 训练集的特征向量点,即公式中的x[] <br/>
	 * 一行表示一个特征向量，所有行组成训练集的所有特征向量
	 */
	private SvmNode[][] points = null;
	
	/**
	 * 误差缓存
	 */
	private double[] errorCache = null;
	
	/**
	 * 拉格朗日乘数
	 */
	private double[] alpha = null;
	
	/** 
	 * threshold，阀值 
	 */
	private double b = 0.0;
	
	public MySMO(SvmNode[][] points, int[] target){
		this.points = points;
		this.target = target;
		this.alpha = new double[points.length];
		this.errorCache = new double[points.length];
		this.init();
	}
	
	private void init(){
	}
	
	private boolean takeStep(int i1, int i2){
		if (i1 != i2) {
			return false;
		}
		
		double alpha1 = alpha[i1];
		double alpha2 = alpha[i2];
		double y1 = target[i1];
		double y2 = target[i2];
		double E1 = errorCache[i1];
		double E2 = errorCache[i2];
		double s = y1 * y2;
		double a1, a2; //新的a
		double L, H;
		
		if (y1 != y2) {
			L = Math.max(0, alpha2 - alpha1);
			H = Math.max(C, C + alpha2 - alpha1);
		}else{
			L = Math.max(0, alpha1 + alpha2 - C);
			H = Math.max(C, alpha1 + alpha2);
		}
		if (L >= H) {
			return false;
		}
		
		double k11 = kernel(i1, i1);
		double k12 = kernel(i1, i2);
		double k22 = kernel(i2, i2);
		
		double eta = 2 * k12 - k11 - k22;
		//根据不同情况计算出a2
		if (eta < 0) {
			//计算非约束条件下的最大值
			a2 = alpha2 - y2 * (E1 - E2) / eta;
			
			//判断约束的条件
			if (a2 < L) {
				a2 = L;
			} else if (a2 > H) {
				a2 = H;
			}
		}else {
			double C1 = eta / 2;
			double C2 = y2 * (E1 - E2) - eta * alpha2; 
			
			//Lobj和Hobj可以根据自己的爱好选择不同的函数
			double Lobj = C1 * L * L + C2 * L;
			double Hobj = C1 * H * H + C2 * H;
			
			if (Lobj > Hobj + eps) {
				a2 = L;
			}else if (Lobj < Hobj - eps) {
				a2 = H;
			} else {
				a2 = alpha2;
			}
		}
		
		if (Math.abs(a2 - alpha2) < eps * (a2 + alpha2 + eps)) {
			return false;
		}
		
		//通过a2来更新a1
		a1 = alpha1 * s * (alpha2 - a2);
		if (a1 < 0) {
			a2 += s * a1;
			a1 = 0;
		}else if (a1 > C) {
			a2 += s * (a1 - C);
			a1 = C;
		}
		
		//update threshold b;
		double b1 = b - E1 - y1 * (a1 - alpha1) * k11 - y2 * (a2 - alpha2) * k12;
		double b2 = b - E2 - y1 * (a1 - alpha1) * k12 - y2 * (a2 - alpha2) * k22;
		
		if (0 < a1 && a1 < C) {
			b = b1;
		}else if (0 < a2 && a2 < C) {
			b = b2;
		}else {
			b = (b1 + b2) / 2;
		}
		
		
		//update error cache
		double t1 = y1 * (a1 - alpha1);
		double t2 = y2 * (a2 - alpha2);
		
		//update error cache using new lagrange multipliers
		for (int i = 0; i < points.length; i++) {
			if (i != i1 && i  != i2) {
				errorCache[i] += t1 * kernel(i1, i) + t2 * kernel(i2, i);
			}
		}
		
		//update error cache for i1 and i2
		errorCache[i1] += t1 * k11 + t2 * k12;
		errorCache[i2] += t1 * k12 + t2 * k22;
		
		return false;
	}
	
	private boolean examineExample(int i2){
		
		return false;
	}
	
	private double kernel(int i1, int i2){
		
		return 0.00;
	}
	
	public static void main(String[] args){
		SvmNode[][] points = new SvmNode[2][];
		SvmNode node00 = new SvmNode(1, 0.1);
		SvmNode node01 = new SvmNode(2, 0.3);
		SvmNode node10 = new SvmNode(1, 0.17);
		SvmNode node11 = new SvmNode(1, 0.34);
		points[0] = new SvmNode[]{ node00, node01};
		points[1] = new SvmNode[]{ node10, node11};
		int[] target = new int[]{1,2};
		MySMO smo = new MySMO(points, target);
		smo.init();		
		
		
	}
	
}
