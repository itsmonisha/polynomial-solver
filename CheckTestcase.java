import java.io.*;
import java.math.BigInteger;
import java.util.*;
import java.util.regex.*;

/*
 CheckTestcase.java
 - Reads JSON from stdin (same format as TestCase1/2)
 - Builds polynomial from first k points (exact rational arithmetic)
 - Evaluates polynomial at ALL provided points and prints mismatches.
 Usage:
   javac CheckTestcase.java
   java CheckTestcase < TestCase1.json
   java CheckTestcase < TestCase2.json
*/

public class CheckTestcase {

    // Exact rational fraction using BigInteger
    static class Fraction {
        BigInteger num, den;
        Fraction(BigInteger n, BigInteger d) {
            if (d.signum() == 0) throw new ArithmeticException("Zero denom");
            if (d.signum() < 0) { n = n.negate(); d = d.negate(); }
            BigInteger g = n.gcd(d);
            if (!g.equals(BigInteger.ONE)) { n = n.divide(g); d = d.divide(g); }
            num = n; den = d;
        }
        Fraction(BigInteger n) { this(n, BigInteger.ONE); }
        Fraction(long n) { this(BigInteger.valueOf(n), BigInteger.ONE); }
        boolean isZero(){ return num.equals(BigInteger.ZERO); }
        Fraction add(Fraction o){ return new Fraction(num.multiply(o.den).add(o.num.multiply(den)), den.multiply(o.den)); }
        Fraction sub(Fraction o){ return new Fraction(num.multiply(o.den).subtract(o.num.multiply(den)), den.multiply(o.den)); }
        Fraction mul(Fraction o){ return new Fraction(num.multiply(o.num), den.multiply(o.den)); }
        Fraction div(Fraction o){ if (o.isZero()) throw new ArithmeticException("Div by zero"); return new Fraction(num.multiply(o.den), den.multiply(o.num)); }
        public String toString(){ if (den.equals(BigInteger.ONE)) return num.toString(); return num + "/" + den; }
    }

    static int charToDigit(char c) {
        if (c >= '0' && c <= '9') return c - '0';
        if (c >= 'a' && c <= 'z') return 10 + (c - 'a');
        if (c >= 'A' && c <= 'Z') return 10 + (c - 'A');
        return -1;
    }

    // parse a string "value" in base `base` into BigInteger
    static BigInteger parseValueInBase(String s, int base){
        s = s.trim();
        BigInteger res = BigInteger.ZERO;
        BigInteger b = BigInteger.valueOf(base);
        for (int i = 0; i < s.length(); ++i){
            int d = charToDigit(s.charAt(i));
            if (d < 0 || d >= base) throw new IllegalArgumentException("Bad digit " + s.charAt(i) + " for base " + base);
            res = res.multiply(b).add(BigInteger.valueOf(d));
        }
        return res;
    }

    public static void main(String[] args) throws Exception {
        // Read entire stdin
        StringBuilder sb = new StringBuilder();
        BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
        String line;
        while ((line = br.readLine()) != null) { sb.append(line).append("\n"); }
        String json = sb.toString();

        // Extract k
        Pattern pk = Pattern.compile("\"k\"\\s*:\\s*(\\d+)");
        Matcher mk = pk.matcher(json);
        if (!mk.find()){ System.err.println("k not found in JSON"); return; }
        int k = Integer.parseInt(mk.group(1));
        if (k <= 0){ System.err.println("invalid k"); return; }

        // Parse entries "1": { "base": "...", "value":"..." }
        Pattern pEntry = Pattern.compile("\"(\\d+)\"\\s*:\\s*\\{[^}]*?\"base\"\\s*:\\s*\"([^\"]+)\"[^}]*?\"value\"\\s*:\\s*\"([^\"]+)\"[^}]*?\\}", Pattern.DOTALL);
        Matcher me = pEntry.matcher(json);

        TreeMap<Integer, BigInteger> points = new TreeMap<>();
        while (me.find()){
            int key = Integer.parseInt(me.group(1));
            int base = Integer.parseInt(me.group(2).trim());
            String val = me.group(3).trim();
            BigInteger y = parseValueInBase(val, base);
            points.put(key, y);
        }

        if (points.size() < k){
            System.err.println("Not enough points provided: have " + points.size() + ", need k=" + k);
            return;
        }

        // SELECT first k points (sorted by key). If you want different selection, you can change this.
        int idx = 0;
        int[] xs = new int[k];
        BigInteger[] ys = new BigInteger[k];
        for (Map.Entry<Integer, BigInteger> e : points.entrySet()){
            if (idx >= k) break;
            xs[idx] = e.getKey();
            ys[idx] = e.getValue();
            idx++;
        }

        // Build augmented matrix (k x (k+1)) of Fractions
        Fraction[][] mat = new Fraction[k][k+1];
        for (int i = 0; i < k; ++i){
            BigInteger x = BigInteger.valueOf(xs[i]);
            BigInteger xp = BigInteger.ONE;
            for (int j = 0; j < k; ++j){
                mat[i][j] = new Fraction(xp);
                xp = xp.multiply(x);
            }
            mat[i][k] = new Fraction(ys[i]);
        }

        // Gaussian elimination (exact)
        for (int col = 0; col < k; ++col){
            int pivot = -1;
            for (int r = col; r < k; ++r) if (!mat[r][col].isZero()){ pivot = r; break; }
            if (pivot == -1){
                System.err.println("Singular matrix at column " + col + ". Try a different subset of k points.");
                return;
            }
            if (pivot != col){
                Fraction[] tmp = mat[pivot];
                mat[pivot] = mat[col];
                mat[col] = tmp;
            }
            Fraction piv = mat[col][col];
            for (int j = col; j <= k; ++j) mat[col][j] = mat[col][j].div(piv);
            for (int r = 0; r < k; ++r){
                if (r == col) continue;
                Fraction factor = mat[r][col];
                if (factor.isZero()) continue;
                for (int j = col; j <= k; ++j){
                    mat[r][j] = mat[r][j].sub(factor.mul(mat[col][j]));
                }
            }
        }

        // Extract coefficients a0..a_{k-1} as Fractions
        Fraction[] coeffs = new Fraction[k];
        for (int i = 0; i < k; ++i) coeffs[i] = mat[i][k];

        // Evaluate polynomial at all provided points; collect wrong ones
        ArrayList<Integer> wrongXs = new ArrayList<>();
        for (Map.Entry<Integer, BigInteger> e : points.entrySet()){
            int xInt = e.getKey();
            BigInteger yGiven = e.getValue();
            BigInteger xp = BigInteger.ONE;
            Fraction sum = new Fraction(BigInteger.ZERO);
            for (int j = 0; j < k; ++j){
                Fraction term = coeffs[j].mul(new Fraction(xp));
                sum = sum.add(term);
                xp = xp.multiply(BigInteger.valueOf(xInt));
            }
            // check sum == yGiven (as fraction)
            if (!sum.num.equals(yGiven.multiply(sum.den))){
                wrongXs.add(xInt);
            }
        }

        // Print results in a clear format you can paste in the form
        if (wrongXs.isEmpty()){
            System.out.println("No wrong data points. All provided points are consistent with the polynomial built from the first k points.");
        } else {
            System.out.print("Wrong Data Set Points (x values): ");
            for (int i = 0; i < wrongXs.size(); ++i){
                if (i > 0) System.out.print(", ");
                System.out.print(wrongXs.get(i));
            }
            System.out.println();
        }
    }
}
