/**
 * 
 */
package jpp.rungekutta;

import java.math.BigDecimal;

/**
 * 
 * @author poa32kc
 *
 */
public class RungeKutta {
	BigDecimal hello = new BigDecimal("2.9");
	private double G = 2.959122 * Math.pow(10, -4.0);
	private double m1 = 1;
	private double m2 = 3.23 * Math.pow(10, -7.0);

	/**
	 * Calculates the delta value from the given position, velocity and interval
	 * 
	 * @param position
	 * @param velocity
	 * @param interval
	 * @return
	 */
	public double[] dglAcceleration(double[] position, double[] velocity,
			double interval) {
		double mu = this.G * (this.m1 + this.m2);
		double[] position_next = substractVector(position,
				scalarProduct(velocity, interval));
		double[] acceleration = scalarProduct(position_next,
				(mu * -1.0) / Math.pow(calculateNorm(position_next), 3.0));
		return acceleration;
	}

	/**
	 * 
	 * @param position
	 * @param velocity
	 * @param interval
	 * @return
	 */
	public double[] rungeKuttaVelocity(double[] position, double[] velocity,
			double interval) {
		double[] f1 = scalarProduct(
				dglAcceleration(position, velocity, 0), interval);
		double[] f2_param = addVector(velocity, scalarProduct(f1, 1.0 / 4.0));
		double[] f2 = scalarProduct(
				dglAcceleration(position, f2_param, (interval / 4.0)),
				interval);
		double[] f3_param = addVector(velocity, scalarProduct(f2, 1.0 / 8.0),
				scalarProduct(f1, 1.0 / 8.0));
		double[] f3 = scalarProduct(
				dglAcceleration(position, f3_param, (interval / 4.0)),
				interval);
		double[] f4_param = addVector(velocity, scalarProduct(f2, -1.0 / 2.0),
				f3);
		double[] f4 = scalarProduct(
				dglAcceleration(position, f4_param, (interval / 2.0)),
				interval);
		double[] f5_param = addVector(velocity, scalarProduct(f1, 3.0 / 16.0),
				scalarProduct(f4, 9.0 / 16.0));
		double[] f5 = scalarProduct(
				dglAcceleration(position, f5_param, (interval * 3.0
						/ 4.0)), interval);
		double[] f6_param = addVector(velocity, scalarProduct(f1, -3.0 / 7.0),
				scalarProduct(f2, 2.0 / 7.0), scalarProduct(f3, 12.0 / 7.0),
				scalarProduct(f4, -12.0 / 7.0), scalarProduct(f5, 8.0 / 7.0));
		double[] f6 = scalarProduct(
				dglAcceleration(position, f6_param, interval),
				interval);
		double[] dY = calculateDY(f1, f3, f4, f5, f6);
		double[] result = addVector(velocity, dY);
		return result;
	}

	public double[] calculateDY(double[]... f) {
		double[] result = new double[f[0].length];
		double[] scalarParams = { 7.0, 32.0, 12.0, 32.0, 7.0 };
		int i = 0;
		for (double[] fVector : f) {
			result = addVector(result, scalarProduct(fVector, scalarParams[i]));
			i = i + 1;
		}
		return scalarProduct(result, 1.0 / 90.0);
	}

	/**
	 * 
	 * @param position
	 * @param velocity
	 * @param interval
	 * @return
	 */
	public double[] rungeKuttaPosition(double[] position, double[] velocity,
			double interval) {
		return null;
	}

	/**
	 * 
	 * @param m1
	 * @param m2
	 * @param G
	 * @return double mu
	 */
	public double calculateMu(double m1, double m2, double G) {
		return -G * m1 * m2;
	}

	/**
	 * 
	 * @param v
	 * @return double[] vectors added.
	 */
	public double[] addVector(double[]... v) {
		double[] result = new double[v[0].length];
		for (int i = 0; i < result.length; i++) {
			for (double[] vector : v) {
				result[i] = result[i] + vector[i];
			}
		}
		return result;
	}

	/**
	 * 
	 * @param v1
	 * @return
	 */
	public double calculateNorm(double[] v1) {
		double result = 0.0;
		for (int i = 0; i < v1.length; i++) {
			result = result + Math.pow(v1[i], 2);
		}
		return Math.sqrt(result);
	}

	/**
	 * 
	 * @param v1
	 * @param v2
	 * @return
	 */
	public double[] substractVector(double[] v1, double[] v2) {
		return addVector(v1, scalarProduct(v2, -1));
	}

	/**
	 * 
	 * @param v
	 * @param s
	 * @return double[] scalar product
	 */
	public double[] scalarProduct(double[] v, double s) {
		double[] result = new double[v.length];
		for (int i = 0; i < v.length; i++) {
			result[i] = v[i] * s;
		}
		return result;
	}

	/**
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		double[] v1 = { 1, 2, 3 };
		double[] v2 = { 3, 4, 5 };
		double[] v3 = { 3, 4, 10 };
		double s = 100;
		RungeKutta rk = new RungeKutta();
		double[] res = rk.addVector(v1, v2, v3);
		double[] sk = rk.scalarProduct(v1, s);
		double a = 3.0 / 4.0;
		System.out.println(a);
		for (int i = 0; i < res.length; i++) {
			System.out.println(res[i]);
		}
	}

}
