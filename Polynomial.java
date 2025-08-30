import java.io.*;
import java.math.BigInteger;
import java.util.*;
import java.util.regex.*;

/**
 * PolynomialSolver
 * - Reads JSON from stdin (format in problem statement)
 * - Each top-level numeric key (like "1", "2", "6") is used as x-coordinate (int)
 * - Each entry contains "base" and "value" where "value" is encoded in that base -> y
 * - Uses first k points (sorted by key) to solve for k coefficients (degree = k-1)
 * - Solves exactly using rational numbers (BigInteger numerator/denominator)
 *
 * Usage:
 *   javac PolynomialSolver.java
 *   java PolynomialSolver < input.json
 */
public class Polynomial {

    // Fraction class using BigInteger
    static class Fraction {
        BigInteger num; // numerator
        BigInteger den; // denominator, always positive

        Fraction(BigInteger n, BigInteger d) {
            if (d.signum() == 0) throw new ArithmeticException("Denominator zero");
            if (d.signum() < 0) { n = n.negate(); d = d.negate(); }
            BigInteger g = n.gcd(d);
            if (!g.equals(BigInteger.ONE)) {
                n = n.divide(g);
                d = d.divide(g);
            }
            this.num = n;
            this.den = d;
        }
        Fraction(long n) { this(BigInteger.valueOf(n), BigInteger.ONE); }
        Fraction(BigInteger n) { this(n, BigInteger.ONE); }

        boolean isZero() { return num.equals(BigInteger.ZERO); }

        Fraction add(Fraction o) {
            BigInteger n = num.multiply(o.den).add(o.num.multiply(den));
            BigInteger d = den.multiply(o.den);
            return new Fraction(n, d);
        }
        Fraction sub(Fraction o) {
            BigInteger n = num.multiply(o.den).subtract(o.num.multiply(den));
            BigInteger d = den.multiply(o.den);
            return new Fraction(n, d);
        }
        Fraction mul(Fraction o) {
            return new Fraction(num.multiply(o.num), den.multiply(o.den));
        }
        Fraction div(Fraction o) {
            if (o.isZero()) throw new ArithmeticException("Divide by zero fraction");
            return new Fraction(num.multiply(o.den), den.multiply(o.num));
        }
        public String toString() {
            if (den.equals(BigInteger.ONE)) return num.toString();
            return num.toString() + "/" + den.toString();
        }
    }

    // parse base-converted string to BigInteger
    static BigInteger parseValueInBase(String s, int base) {
        s = s.trim();
        BigInteger res = BigInteger.ZERO;
        BigInteger b = BigInteger.valueOf(base);
        for (int i = 0; i < s.length(); ++i) {
            char c = s.charAt(i);
            int digit = charToDigit(c);
            if (digit < 0 || digit >= base)
                throw new IllegalArgumentException("Digit '" + c + "' not valid for base " + base);
            res = res.multiply(b).add(BigInteger.valueOf(digit));
        }
        return res;
    }

    static int charToDigit(char c) {
        if (c >= '0' && c <= '9') return c - '0';
        if (c >= 'a' && c <= 'z') return 10 + (c - 'a');
        if (c >= 'A' && c <= 'Z') return 10 + (c - 'A');
        return -1;
    }

    public static void main(String[] args) throws Exception {
        // Read input from stdin
        StringBuilder sb = new StringBuilder();
        BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
        String line;
        while ((line = br.readLine()) != null) {
            sb.append(line).append("\n");
        }
        String json = sb.toString();

        // Extract k
        Pattern pk = Pattern.compile("\"k\"\\s*:\\s*(\\d+)");
        Matcher mk = pk.matcher(json);
        if (!mk.find()) {
            System.err.println("Couldn't find k in input.");
            return;
        }
        int k = Integer.parseInt(mk.group(1));
        if (k <= 0) {
            System.err.println("k must be > 0");
            return;
        }

        // find numeric top-level entries: "1": { "base": "10", "value":"4" }
        Pattern pEntry = Pattern.compile("\"(\\d+)\"\\s*:\\s*\\{[^}]*?\"base\"\\s*:\\s*\"([^\"]+)\"[^}]*?\"value\"\\s*:\\s*\"([^\"]+)\"[^}]*?\\}", Pattern.DOTALL);
        Matcher me = pEntry.matcher(json);

        TreeMap<Integer, BigInteger> points = new TreeMap<>();
        while (me.find()) {
            int key = Integer.parseInt(me.group(1));
            String baseStr = me.group(2).trim();
            String valStr = me.group(3).trim();
            int base = Integer.parseInt(baseStr);
            BigInteger y = parseValueInBase(valStr, base);
            points.put(key, y);
        }

        if (points.size() < k) {
            System.err.println("Not enough points provided: need k=" + k + " but found " + points.size());
            return;
        }

        // pick first k points (sorted by key)
        int idx = 0;
        int[] xs = new int[k];
        BigInteger[] ys = new BigInteger[k];
        for (Map.Entry<Integer, BigInteger> e : points.entrySet()) {
            if (idx >= k) break;
            xs[idx] = e.getKey();
            ys[idx] = e.getValue();
            idx++;
        }

        // Build augmented matrix of Fractions (k rows, k+1 columns)
        Fraction[][] mat = new Fraction[k][k + 1];
        for (int i = 0; i < k; ++i) {
            BigInteger x = BigInteger.valueOf(xs[i]);
            BigInteger xp = BigInteger.ONE; // x^0 initially
            for (int j = 0; j < k; ++j) {
                mat[i][j] = new Fraction(xp); // xp / 1
                xp = xp.multiply(x); // next power
            }
            mat[i][k] = new Fraction(ys[i]); // right-hand side
        }

        // Gaussian elimination (exact fractions)
        for (int col = 0; col < k; ++col) {
            // find pivot row
            int pivot = -1;
            for (int r = col; r < k; ++r) {
                if (!mat[r][col].isZero()) { pivot = r; break; }
            }
            if (pivot == -1) {
                System.err.println("Matrix singular or no unique solution (pivot at column " + col + " is zero).");
                return;
            }
            // swap rows if needed
            if (pivot != col) {
                Fraction[] tmp = mat[pivot];
                mat[pivot] = mat[col];
                mat[col] = tmp;
            }

            // normalize pivot row
            Fraction pivotVal = mat[col][col];
            for (int j = col; j <= k; ++j) {
                mat[col][j] = mat[col][j].div(pivotVal);
            }

            // eliminate other rows
            for (int r = 0; r < k; ++r) {
                if (r == col) continue;
                Fraction factor = mat[r][col];
                if (factor.isZero()) continue;
                for (int j = col; j <= k; ++j) {
                    Fraction prod = factor.mul(mat[col][j]);
                    mat[r][j] = mat[r][j].sub(prod);
                }
            }
        }

        // Read solution coefficients (a0..a_{k-1})
        System.out.println("Coefficients (a0 ... a" + (k - 1) + "):");
        for (int i = 0; i < k; ++i) {
            System.out.print(mat[i][k].toString());
            if (i < k - 1) System.out.print(" ");
        }
        System.out.println();

        // also print polynomial in human readable form
        System.out.println("\nPolynomial:");
        StringBuilder poly = new StringBuilder();
        for (int i = 0; i < k; ++i) {
            String coeff = mat[i][k].toString();
            if (i == 0) {
                poly.append(coeff);
            } else {
                poly.append(" + ").append(coeff).append("*x");
                if (i > 1) poly.append("^").append(i);
            }
        }
        System.out.println(poly.toString());
    }
}
