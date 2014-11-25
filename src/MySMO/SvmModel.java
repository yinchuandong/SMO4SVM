package MySMO;

public class SvmModel {
	
	private double[] alpha;
	private int[] y;
	private double b;
	
	public SvmModel(double[] alpha, int[] y, double b){
		this.alpha = alpha;
		this.y = y;
		this.b = b;
	}

	public double[] getAlpha() {
		return alpha;
	}

	public void setAlpha(double[] alpha) {
		this.alpha = alpha;
	}

	public int[] getY() {
		return y;
	}

	public void setY(int[] y) {
		this.y = y;
	}

	public double getB() {
		return b;
	}

	public void setB(double b) {
		this.b = b;
	}
	
	
}
