package MySMO;

import java.util.Random;

import sun.security.pkcs11.Secmod.DbMode;

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
	private int[] y = null;
	
	/**
	 * 训练集的特征向量点,即公式中的x[] <br/>
	 * 一行表示一个特征向量，所有行组成训练集的所有特征向量
	 */
	private SvmNode[][] x = null;
	
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
	
	/**
	 * rbf kernel for exp(-gamma*|u-v|^2), 默认为0.5，也可设为1/num
	 */
	private double gamma = 0.5;
	
	/**
	 * 对points点积的缓存
	 */
	private double[][] dotDache = null;
	
	/**
	 * 所有向量的数目
	 */
	private int N = 0;
	
	private Random random = null;
	
	public MySMO(SvmNode[][] x, int[] y){
		this.x = x;
		this.y = y;
		this.N = x.length;
		this.alpha = new double[N];
		this.errorCache = new double[N];
		this.dotDache = new double[N][N];
		this.random = new Random();
		this.init();
	}
	
	private void init(){
		//初始化点积dotCache
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				this.dotDache[i][j] = dot(x[i], x[j]);
			}
		}
	}
	
	private boolean takeStep(int i1, int i2){
		if (i1 == i2) {
			return false;
		}
		
		double alpha1 = alpha[i1];
		double alpha2 = alpha[i2];
		double y1 = y[i1];
		double y2 = y[i2];
		double E1 = 0;
		double E2 = 0;
		double s = y1 * y2;
		double a1, a2; //新的a
		double L, H;
		
		if (0 < alpha1 && alpha1 < C) {
			E1 = errorCache[i1];
		}else{
			E1 = learnFunc(i1) - y1;
		}
		
		if (0 < alpha2 && alpha2 < C) {
			E2 = errorCache[i2];
		}else{
			E2 = learnFunc(i2) - y2;
		}
		
		if (y1 != y2) {
			L = Math.max(0, alpha2 - alpha1);
			H = Math.min(C, C + alpha2 - alpha1);
		}else{
			L = Math.max(0, alpha1 + alpha2 - C);
			H = Math.min(C, alpha1 + alpha2);
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
		a1 = alpha1 + s * (alpha2 - a2);
		
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
		
		//the other way to update b
//		double b1 = b + E1 + y1 * (a1 - alpha1) * k11 + y2 * (a2 - alpha2) * k12;
//		double b2 = b + E2 + y1 * (a1 - alpha1) * k12 + y2 * (a2 - alpha2) * k22;
		
		double bNew = 0;
		double deltaB = 0;
		if (0 < a1 && a1 < C) {
			bNew = b1;
		}else if (0 < a2 && a2 < C) {
			bNew = b2;
		}else {
			bNew = (b1 + b2) / 2;
		}
		deltaB = bNew - this.b; //b的增量
		
		//update error cache
		double t1 = y1 * (a1 - alpha1);
		double t2 = y2 * (a2 - alpha2);
		
		//update error cache using new lagrange multipliers
		for (int i = 0; i < N; i++) {
			if (0 < alpha[i] && alpha[i] < C) { // condition in i != i1 && i != i2
				errorCache[i] += t1 * kernel(i1, i) + t2 * kernel(i2, i) - deltaB;
			}
		}
		
		//update error cache for i1 and i2
		errorCache[i1] += t1 * k11 + t2 * k12;
		errorCache[i2] += t1 * k12 + t2 * k22;
		
//		errorCache[i1] = 0.0;
//		errorCache[i2] = 0.0;
		
		//store a1, a2 in alpha array
		alpha[i1] = a1;
		alpha[i2] = a2;
		
		
		return true;
	}
	
	/**
	 * 检查最好的样本，并进行takeStep计算
	 * @param i1
	 * @return
	 */
	private boolean examineExample(int i1){
		double y1 = y[i1];
		double alpha1 = alpha[i1];
		double E1 = 0;
		
		if (0 < alpha1 && alpha1 < C) {
			E1 = errorCache[i1];
		}else{
			E1 = learnFunc(i1) - y1;
		}
		
		double r1 = y1 * E1;
		if ((r1 < -tolerance && alpha1 < C) || (r1 > tolerance && alpha1 > 0)) {

			//选择 E1 - E2 差最大的两点
			int i2 = this.findMax(E1);
			if (i2 >= 0) {
				if (takeStep(i1, i2)) {
					return true;
				}
			}
			
			//先选择 0 < alpha < C的点
			int k0 = randomSelect(i1);
			for (int k = k0; k < N + k0; k++) {
				i2 = k % N;
				if (0 < alpha[i2] && alpha[i2] < C) {
					if (takeStep(i1, i2)) {
						return true;
					}
				}
			}
			
			//如果不符合，再遍历全部点
			k0 = randomSelect(i1);
			for (int k = k0; k < N + k0; k++) {
				i2 = k % N;
				if (takeStep(i1, i2)) {
					return true;
				}
			}
			
		}
		
		return false;
	}
	
	
	private void train(){
		System.out.println("begin train");
		
		int maxIter = 500;
		int iterCount = 0;
		int numChanged = 0;
		boolean examineAll = true;
		
		while((iterCount < maxIter) && (numChanged > 0 || examineAll)){
			numChanged = 0;
			
			if (examineAll) {
				for (int i = 0; i < N; i++) {
					if (examineExample(i)) {
						numChanged ++;
					}
				}
			}else{
				for (int i = 0; i < N; i++) {
					if (alpha[i] != 0 && alpha[i] != C) {
						if (examineExample(i)) {
							numChanged ++;
						}
					}
				}
			}
			
			iterCount ++;
			if (examineAll) {
				examineAll = false;
			}else if (numChanged == 0) {
				examineAll = true;
			}
		}
		
		System.out.println("end of train");
	}
	
	
	/**
	 * 找到|E1 - E2|差最大的点的下标
	 * @param E1
	 * @return
	 */
	private int findMax(double E1){
		int i2 = -1;
		double tmax = 0.0;
		for (int k = 0; k < N; k++) {
			if (0 < alpha[k] && alpha[k] < C) {
				double E2 = errorCache[k];
				double tmp = Math.abs(E2 - E1);
				if (tmp > tmax) {
					tmax = tmp;
					i2 = k;
				}
			}
		}
		
		return i2;
	}
	
	/**
	 * 随机选择i2，但要求i1 != i2
	 * @param i1
	 * @return
	 */
	private int randomSelect(int i1){
		int i2 = 0;
		do {
			i2 = random.nextInt(N);
		} while (i1 == i2);
		return i2;
	}
	
	/**
	 * 核函数 exp(-gamma*|u-v|^2)
	 * @param i1
	 * @param i2
	 * @return
	 */
	private double kernel(int i1, int i2){
		double result = 0.0;
		result = Math.exp(-gamma * (dotDache[i1][i1] + dotDache[i2][i2] - 2 * dotDache[i1][i2]));
		
		return result;
	}
	
	/**
	 * 对两个向量进行点积，需要x,y向量都按照index升序排序
	 * @param x
	 * @param y
	 * @return
	 */
	private double dot(SvmNode[] x, SvmNode[] y){
		double sum = 0.0;
		int xLen = x.length;
		int yLen = y.length;
		int i = 0;
		int j = 0;
		
		while(i < xLen && j < yLen){
			if (x[i].getIndex() == y[j].getIndex()) {
				sum += x[i].getValue() * y[j].getValue();
				i++;
				j++;
			}else{
				if (x[i].getIndex() > y[j].getIndex()) {
					j++;
				}else{
					i++;
				}
			}
		}
		return sum;
	}
	
	/**
	 * 学习函数u，算误差的时候要用
	 * @param k
	 * @return
	 */
	private double learnFunc(int k){
		double sum = 0.0;
		for (int i = 0; i < N; i++) {
			sum += alpha[i]*y[i]*kernel(i, k);
		}
		sum += this.b;
		return sum;
	}
	
	public static void main(String[] args){
		
		long start = System.currentTimeMillis();
		
		SvmData data = FileUtil.loadTrainFile("heart_scale");
		MySMO smo = new MySMO(data.getX(), data.getY());
		smo.train();
		
		long end = System.currentTimeMillis();
		double delay = (double)(end - start) / 1000.00;
		System.out.println("耗时：" + delay + "s");
	}
	
}
