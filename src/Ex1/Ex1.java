package Ex1;

import java.util.Arrays;

/**
 * Introduction to Computer Science 2026, Ariel University,
 * Ex1: arrays, static functions and JUnit
 * https://docs.google.com/document/d/1GcNQht9rsVVSt153Y8pFPqXJVju56CY4/edit?usp=sharing&ouid=113711744349547563645&rtpof=true&sd=true
 * <p>
 * This class represents a set of static methods on a polynomial functions - represented as an array of doubles.
 * The array {0.1, 0, -3, 0.2} represents the following polynomial function: 0.2x^3-3x^2+0.1
 * This is the main Class you should implement (see "add your code below")
 *
 * @author boaz.benmoshe
 */
public class Ex1 {
    /**
     * Epsilon value for numerical computation, it serves as a "close enough" threshold.
     */
    public static final double EPS = 0.001; // the epsilon to be used for the root approximation.
    /**
     * The zero polynomial function is represented as an array with a single (0) entry.
     */
    public static final double[] ZERO = {0};
    /**
     * Computes the f(x) value of the polynomial function at x.
     *
     * @param poly - polynomial function
     * @param x
     * @return f(x) - the polynomial function value at x.
     */
    public int test = 9;

    public static double f(double[] poly, double x) {
        double ans = 0;
        for (int i = 0; i < poly.length; i++) {
            double c = Math.pow(x, i);
            ans += c * poly[i];
        }
        return ans;
    }

    /**
     * Given a polynomial function (p), a range [x1,x2] and an epsilon eps.
     * This function computes an x value (x1<=x<=x2) for which |p(x)| < eps,
     * assuming p(x1)*p(x2) <= 0.
     * This function should be implemented recursively.
     *
     * @param p   - the polynomial function
     * @param x1  - minimal value of the range
     * @param x2  - maximal value of the range
     * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
     * @return an x value (x1<=x<=x2) for which |p(x)| < eps.
     */
    public static double root_rec(double[] p, double x1, double x2, double eps) {
        double f1 = f(p, x1);
        double x12 = (x1 + x2) / 2;
        double f12 = f(p, x12);
        if (Math.abs(f12) < eps) {
            return x12;
        }
        if (f12 * f1 <= 0) {
            return root_rec(p, x1, x12, eps);
        } else {
            return root_rec(p, x12, x2, eps);
        }
    }

    /**
     * This function creates a polynomial that exactly fits a given set of 2D points.
     * It supports:
     * 2 points → returns a linear polynomial (degree 1)
     * 3 points → returns a quadratic polynomial (degree 2)
     *
     * @param xx double[]
     * @param yy double[]
     * @return an array of doubles representing the coefficients of the polynom.
     */
    public static double[] PolynomFromPoints(double[] xx, double[] yy) {
        double[] ans = null;
        int lx = xx.length;
        int ly = yy.length;
        if (xx != null && yy != null && lx == ly && lx > 1 && lx < 4) {
            if (lx == 2) {
                double x1 = xx[0], y1 = yy[0];
                double x2 = xx[1], y2 = yy[1];
                double a = (y2 - y1) / (x2 - x1);
                double b = y1 - a * x1;
                ans = new double[]{b, a};
            } else {
                double x1 = xx[0], y1 = yy[0];
                double x2 = xx[1], y2 = yy[1];
                double x3 = xx[2], y3 = yy[2];
                double denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
                if (Math.abs(denom) < EPS) return null;
                double a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
                double b = (x3 * x3 * (y1 - y2) + x2 * x2 * (y3 - y1) + x1 * x1 * (y2 - y3)) / denom;
                double c = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;
                ans = new double[]{c, b, a};
            }
            ans = removeUnnecessaryZeros(ans);
        }
        return ans;
    }

    /**
     * This function checks whether two polynomials represent the same mathematical polynomial function.
     * It compares their coefficients after removing unnecessary trailing zeros.
     *
     * @param p1 first polynomial function double[]
     * @param p2 second polynomial function double[]
     * @return true iff p1 represents the same polynomial function as p2.
     */
    public static boolean equals(double[] p1, double[] p2) {
        if (p1 == null || p2 == null) return false;

        p1 = removeUnnecessaryZeros(p1);
        p2 = removeUnnecessaryZeros(p2);
        int mek1 = Mekadem(p1);
        int mek2 = Mekadem(p2);
        if (mek1 != mek2 || p1.length != p2.length) return false;
        int maxmek = Math.max(mek1, mek2);
        for (int i = 0; i <= maxmek; i++) {
            if (Math.abs(p1[i] - p2[i]) > EPS) {
                return false;
            }
        }
        return true;
    }

    /**
     * This function determines the degree of a polynomial based on its coefficients.
     * The degree is the highest index in the array whose value is not zero.
     *
     * @param p double[]
     * @return The degree of the polynomial
     */
    public static int Mekadem(double[] p) {
        for (int i = p.length - 1; i >= 0; i--) {
            if (Math.abs(p[i]) > EPS) {
                return i;
            }
        }
        return 0;
    }

    /**
     * This function removes trailing zeros from a polynomial’s coefficient array.
     * Trailing zeros do not change the value of the polynomial.
     *
     * @param arr double[]
     * @return The function returns a new array with those extra zeros removed.
     */
    public static double[] removeUnnecessaryZeros(double[] arr) {
        if (polinomZero(arr)) {
            return new double[]{0};
        }
        int count_zeros = 0;
        for (int i = arr.length - 1; i >= 0; i--) {
            if (arr[i] == 0) {
                count_zeros++;
            } else {
                break;
            }
        }
        double[] ans = new double[arr.length - count_zeros];
        for (int i = 0; i < ans.length; i++) {
            ans[i] = arr[i];
        }
        return ans;
    }

    /**
     * This function checks whether a given polynomial is the zero polynomial.
     * A zero polynomial is a polynomial where all coefficients are equal to zero (or very close to zero within a small tolerance)
     *
     * @param p1 double[]
     *           The array of polynomial coefficients. Index represents the power of x.
     * @return true or false
     */
    public static boolean polinomZero(double[] p1) {
        if (p1 == null || p1.length == 0) {
            return false;
        }
        for (int i = 0; i < p1.length; i++) {
            if (Math.abs(p1[i]) > EPS) {
                return false;
            }
        }
        return true;
    }

    /**
     * Computes a String representing the polynomial function.
     * For example the array {2,0,3.1,-1.2} will be presented as the following String  "-1.2x^3 +3.1x^2 +2.0"
     * @param poly the polynomial function represented as an array of doubles
     * @return String representing the polynomial function:
     */
    public static String poly(double[] poly) {
        if (poly == null || poly.length == 0) return "0";
        String ans = "";
        boolean first = true;
        for (int i = poly.length - 1; i >= 0; i--) {
            double coef = poly[i];
            if (Math.abs(coef) < EPS) continue;
            if (!first) {
                if (coef >= 0) ans += " +";
                else ans += " -";
            } else {
                if (coef < 0) ans += "-";
                first = false;
            }
            double c = Math.abs(coef);
            if (i == 0) {
                ans += c;
            } else if (i == 1) {
                ans += c + "x";
            } else {
                if (c != 1) {
                    ans += c;
                }
                ans += "x^" + i;
            }
        }
        if (ans.equals("")) return "0";
        return ans;
    }


    /**
     * Given two polynomial functions (p1,p2), a range [x1,x2] and an epsilon eps. This function computes an x value (x1<=x<=x2)
     * for which |p1(x) -p2(x)| < eps, assuming (p1(x1)-p2(x1)) * (p1(x2)-p2(x2)) <= 0.
     *
     * @param p1  - first polynomial function
     * @param p2  - second polynomial function
     * @param x1  - minimal value of the range
     * @param x2  - maximal value of the range
     * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
     * @return an x value (x1<=x<=x2) for which |p1(x) - p2(x)| < eps.
     */
    public static double sameValue(double[] p1, double[] p2, double x1, double x2, double eps) {
        if (p1 == null || p2 == null || eps <= 0 || x1 > x2) return -1;
        double ans = x1;
        double y1 = f(p1, x1) - f(p2, x1);
        double y2 = f(p1, x2) - f(p2, x2);
        if (y1 * y2 > 0) return -1;
        while (Math.abs(x2 - x1) > eps) {
            double xm = (x1 + x2) / 2.0;
            double ym = f(p1, xm) - f(p2, xm);
            if (Math.abs(ym) < eps) return xm;
            if (y1 * ym < 0) {
                x2 = xm;
                y2 = ym;
            } else {
                x1 = xm;
                y1 = ym;
            }
        }
        ans = (x1 + x2) / 2.0;
        return ans;
    }

    /**
     * Given a polynomial function (p), a range [x1,x2] and an integer with the number (n) of sample points.
     * This function computes an approximation of the length of the function between f(x1) and f(x2)
     * using n inner sample points and computing the segment-path between them.
     * assuming x1 < x2.
     * This function should be implemented iteratively (none recursive).
     *
     * @param p                - the polynomial function
     * @param x1               - minimal value of the range
     * @param x2               - maximal value of the range
     * @param numberOfSegments - (A positive integer value (1,2,...).
     * @return the length approximation of the function between f(x1) and f(x2).
     */
    public static double length(double[] p, double x1, double x2, int numberOfSegments) {
        double ans = 0;
        if (p == null || numberOfSegments <= 0 || x1 >= x2) return -1;
        double segment = (x2 - x1) / numberOfSegments;
        double x = x1;
        double y = f(p, x);
        for (int i = 1; i <= numberOfSegments; i++) {
            double xCurr = x1 + i * segment;
            double yCurr = f(p, xCurr);
            ans += Math.sqrt(Math.pow(xCurr - x, 2) + Math.pow(yCurr - y, 2));
            x = xCurr;
            y = yCurr;
        }
        return ans;
    }

    /**
     * This function estimates the area between two curves represented by the polynomials p1 and p2, over the interval from x1 to x2.
     *
     * @param p1                - first polynomial function
     * @param p2                - second polynomial function
     * @param x1                - minimal value of the range
     * @param x2                - maximal value of the range
     * @param numberOfTrapezoid - a natural number representing the number of Trapezoids between x1 and x2.
     * @return the approximated area between the two polynomial functions within the [x1,x2] range.
     */
    public static double area(double[] p1, double[] p2, double x1, double x2, int numberOfTrapezoid) {
        double ans = 0;
        if (p1 == null || p2 == null || numberOfTrapezoid <= 0 || x2 <= x1) return -1;
        double width = (x2 - x1) / numberOfTrapezoid;
        for (int i = 0; i < numberOfTrapezoid; i++) {
            double xA = x1 + i * width;
            double xB = x1 + (i + 1) * width;
            double yA = Math.abs(f(p1, xA) - f(p2, xA));
            double yB = Math.abs(f(p1, xB) - f(p2, xB));
            ans += width * (yA + yB) / 2.0;
        }
        return ans;
    }

    /**
     * This function computes the array representation of a polynomial function from a String
     * representation. Note:given a polynomial function represented as a double array,
     * getPolynomFromString(poly(p)) should return an array equals to p.
     *
     * @param p - a String representing polynomial function.
     * @return Polynomial coefficients where index = exponent
     */
    public static double[] getPolynomFromString(String p) {
        p = normalize(p);
        if (p.equals("")) return new double[]{0};

        int max = extractMaxPower(p);
        double[] ans = new double[max + 1];

        int[] iRef = {0};
        while (iRef[0] < p.length()) {
            double coef = parseCoef(p, iRef);
            int pow = parsePow(p, iRef);
            ans[pow] += coef;
        }
        ans = removeUnnecessaryZeros(ans);
        return ans;
    }

    /**
     * This function prepares a polynomial string for parsing or further processing.
     * It removes spaces and ensures that the expression starts with either a + or - sign.
     * @param s String
     * @return A string representing a polynomial expression
     */
    public static String normalize(String s) {
        if (s == null) return "";
        s = s.replace(" ", "");
        if (s.length() == 0) return "";
        if (s.charAt(0) != '-') s = "+" + s;
        return s;
    }

    /**
     * This function scans a polynomial string and finds the highest power of x that appears in the expression.
     * It helps determine the degree of the polynomial before parsing its coefficients.
     * @param s String
     * @return A polynomial expression in text form
     */
    public static int extractMaxPower(String s) {
        int maxPow = 0;
        for (int i = 0; i < s.length(); i++) {
            if (s.charAt(i) == '^') {
                int pow = 0;
                int j = i + 1;
                while (j < s.length() && Character.isDigit(s.charAt(j))) {
                    pow = pow * 10 + (s.charAt(j) - '0');
                    j++;
                }
                if (pow > maxPow) maxPow = pow;
            } else if (s.charAt(i) == 'x') {
                if (maxPow < 1) maxPow = 1;
            }
        }
        return maxPow;
    }

    /**
     * This function reads and extracts the power (exponent) of the variable x from a polynomial string while parsing it.
     * @param s string
     * @param iRef int[]
     * @return x without a power (e.g., "x" or "3x") and Power explicitly defined (e.g., "x^7")
     */
    public static int parsePow(String s, int[] iRef) {
        int i = iRef[0];
        int pow = 0;
        if (i < s.length() && s.charAt(i) == 'x') {
            pow = 1;
            i++;
            if (i < s.length() && s.charAt(i) == '^') {
                i++;
                pow = 0;
                while (i < s.length() && Character.isDigit(s.charAt(i))) {
                    pow = pow * 10 + (s.charAt(i) - '0');
                    i++;
                }
            }
        }
        iRef[0] = i;
        return pow;
    }

    /**
     * This function extracts the coefficient of a polynomial term from a string while parsing.
     * @param s String
     * @param iRef int[]
     * @return The numeric coefficient for the term
     */
    public static double parseCoef(String s, int[] iRef) {
        int i = iRef[0];
        char sign = s.charAt(i++);
        StringBuilder num = new StringBuilder();
        while (i < s.length() && (Character.isDigit(s.charAt(i)) || s.charAt(i) == '.')) {
            num.append(s.charAt(i));
            i++;
        }
        double coef = (num.length() > 0) ? Double.parseDouble(num.toString()) : 1.0;
        if (sign == '-') coef = -coef;
        iRef[0] = i;
        return coef;
    }

    /**
     * This function computes the polynomial function which is the sum of two polynomial functions (p1,p2)
     *
     * @param p1 double[]
     * @param p2 double[]
     * @return A new polynomial array representing p1 + p2
     */
    public static double[] add(double[] p1, double[] p2) {
        double[] ans = ZERO;//
        if (p1 == null || p2 == null) return ans;
        int max = Math.max(p1.length, p2.length);
        ans = new double[max];
        for (int i = 0; i < max; i++) {
            double a = i < p1.length ? p1[i] : 0;
            double b = i < p2.length ? p2[i] : 0;
            ans[i] = a + b;
        }
        ans = removeUnnecessaryZeros(ans);
        return ans;
    }

    /**
     * This function computes the polynomial function which is the multiplication of two polynoms (p1,p2)
     * It multiplies two polynomials using coefficient convolution and returns the resulting polynomial.
     * @param p1 double[]
     * @param p2 double[]
     * @return The polynomial representing p1 * p2
     */
    public static double[] mul(double[] p1, double[] p2) {
        double[] ans = ZERO;
        if (p1 == null || p2 == null) return ans;
        ans = new double[p1.length + p2.length - 1];
        for (int i = 0; i < p1.length; i++) {
            for (int j = 0; j < p2.length; j++) {
                ans[i + j] += p1[i] * p2[j];
            }
        }
        ans = removeUnnecessaryZeros(ans);
        return ans;
    }

    /**
     * This function computes the derivative of the p0 polynomial function.
     * @param po double[]
     * @return The coefficients of the derivative polynomial
     */
    public static double[] derivative(double[] po) {
        double[] ans = ZERO;
        if (po == null || po.length == 1) return ans;
        ans = new double[po.length - 1];
        for (int i = 1; i < po.length; i++) {
            ans[i - 1] = i * po[i];
        }
        return ans;
    }

    static void main() {
        double[][] d1 = {{0}, {1}, {1, 2, 0, 0}};
        double[][] d2 = {Ex1.ZERO, {1 + Ex1.EPS / 2}, {1, 2}};
        double[][] xx = {{-2 * Ex1.EPS}, {1 + Ex1.EPS * 1.2}, {1, 2, Ex1.EPS / 2}};
        double[] p1 = {1, -3, 2};
        double[] p2 = {0};
        double ans = sameValue(p1, p2, 0, 1, EPS);
        System.out.println(ans);
        double[] p5 = {0, 1};
        double len = length(p5, 0, 10, 10000);
        System.out.println(len);
        double[] p6 = {2, -3, 5};
        double[] nigzeret = derivative(p6);
        System.out.println(Arrays.toString(nigzeret));

    }
}
