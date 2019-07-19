public class NPi {
	private static final long D = 14; // # of digits of working precision
	private static final long M = (long) Math.pow(16, D);
	private static final long shift = 4 * D;
	private static final long mask = M - 1;

	/* Calculates the nth digit of PI in hexadecimal using the BBP algorithm.
	 * Adapted from: https://web.archive.org/web/20150627225748/http://en.literateprograms.org/Pi_with_the_BBP_formula_%28Python%29
	 * With help from: https://math.stackexchange.com/a/1696997
	 * @param n The nth place to calculate
	*/
	public static void BBPPi(int n) 
	{
		n--;
		long x = (4 * S(1, n) - 2 * S(4, n) - S(5, n) - S(6, n)) & mask;
		System.out.format("BBP: %14x\n", x);
	}

	/* Calculates PI to n digits using the Gauss-Legendre algorithm.
	 * The implementation is limited by the 64 bit double precision (53*log(2) ~ 15 digits),
	 * values of n > 15 are rejected as a result of this; A better approach could use a 
	 * bit array or other variable-length data structure.
	 * This algorithm was implemented as a test, I left it in for kicks.
	 * @param n The nth place to calculate
	*/
	public static void GaussLegendrePi(int n)
	{
		if (n > 15) {
			System.out.println("Gauss-Legendre cannot be used with n > 15");
			return;
		}

		// Initialize
		double a = 1;
		double b = 1.0 / Math.sqrt(2);
		double t = 1.0 / 4.0;
		int p = 1;
		double i = Math.pow(10, -n);

		/* Iterate
		 * Note: The loop quickly becomes unstoppable as i (and n by proxy) increases.
		 * Any value of n > 15 will not produce a result because of precision limitations.
		*/
		while (Math.abs(a - b) > i) {
			double an = a;
			double bn = b;
			double tn = t;
			int pn = p;
			a = (an + bn) / 2.0;
			b = Math.sqrt(an * bn);
			t = tn - pn * Math.pow(an - a, 2);
			p = 2 * pn;
		}

		// Print result
		double pi = Math.pow(a + b, 2) / (4 * t);
		System.out.format("Gauss-Legendre: %f\n", pi);
	}

	/* Helper to calculate partial sums used in the BBP forumla.
	 *
	 * @param j The specific term's 'J' variable
	 * @param n The nth place to calculate
	 * @return The final sum
	*/
	private static long S(int j, int n) 
	{
		/* Left sum:
		 * sum_(k=0)^n (16^(n - k) mod (8 k + j))/(8 k + j)
		*/ 
		long sum = 0;
		long k = 0;
		while(k <= n) {
			long r = 8 * k + j;
			sum = (sum + idiv(((long) (Math.pow(16, n - k) % r) << shift), (double) r)) & mask;
			k++;
		}

		/* Right sum (a convergent series):
		 * sum_(k=n + 1)^inf 16^(n - k)/(8 k + j)
		*/
		long total = 0;
		k = n + 1;
		while(true) {
			double xp = Math.round(Math.pow(16, n - k) * M);
			long nTotal = total + idiv(xp, (8 * k + j));
			if (total == nTotal) {
				break;
			} else {
				total = nTotal;
			}
			k++;
		}

		return sum + total;
	}
	
	/* Helper that calculates the fractional part of a number.
	 * Implementation matches the one described in the Wolfram Language: http://mathworld.wolfram.com/FractionalPart.html
	 * @param x A real number
	 * @return The fractional part of x
	*/
	private static double frac(double x)
	{
		if (x >= 0) {
			return x - Math.floor(x);
		} else {
			return x - Math.ceil(x);
		}
	}

	/* Helper that calculates the integer division between two numbers.
	 * @param a The first number
	 * @param b The second number
	 * @return Result of the integer division a / b
	*/
	private static long idiv(double a, double b)
	{
		double x = a / b;
		if (x >= 0) {
			return (long) Math.floor(x);
		} else {
			return (long) Math.ceil(x);
		}
	}

	public static void main(String[] args) 
	{
		if (args.length < 1) {
			System.out.println("First argument must be a number.");
		} else {
			int n = 0;
			try {
				n = Integer.parseInt(args[0]);
			} catch (Exception ex) {
				System.err.println("Invalid number.");
				System.exit(0);
			}
			System.out.format("n = %d\n", n);
			BBPPi(n);
			GaussLegendrePi(n);
		}
	}
}