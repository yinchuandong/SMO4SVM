package MySMO;

public class SvmDataLine {

	private SvmNode[] x;
	private int y;
	
	public SvmDataLine(SvmNode[] x, int y){
		this.x = x;
		this.y = y;
	}

	public SvmNode[] getX() {
		return x;
	}

	public void setX(SvmNode[] x) {
		this.x = x;
	}

	public int getY() {
		return y;
	}

	public void setY(int y) {
		this.y = y;
	}


	
	
}
