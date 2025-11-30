package Ex1;

import org.junit.jupiter.api.Test;

import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.*;

/**
 * * Introduction to Computer Science 2026, Ariel University,
 * * Ex1: arrays, static functions and JUnit
 * <p>
 * This JUnit class represents a JUnit (unit testing) for Ex1-
 * It contains few testing functions for the polynomial functions as define in Ex1.
 * Note: you should add additional JUnit testing functions to this class.
 *
 * @author boaz.ben-moshe
 */

class Ex1Test {
    static final double[] P1 = {2, 0, 3, -1, 0}, P2 = {0.1, 0, 1, 0.1, 3};
    static double[] po1 = {2, 2}, po2 = {-3, 0.61, 0.2};
    ;
    static double[] po3 = {2, 1, -0.7, -0.02, 0.02};
    static double[] po4 = {-3, 0.61, 0.2};

    @Test
    /**
     * Tests that f(x) == poly(x).
     */
    void testF() {
        double fx0 = Ex1.f(po1, 0);
        double fx1 = Ex1.f(po1, 1);
        double fx2 = Ex1.f(po1, 2);
        assertEquals(fx0, 2, Ex1.EPS);
        assertEquals(fx1, 4, Ex1.EPS);
        assertEquals(fx2, 6, Ex1.EPS);
    }

    @Test
    /**
     * Tests that p1(x) + p2(x) == (p1+p2)(x)
     */
    void testF2() {
        double x = Math.PI;
        double[] po12 = Ex1.add(po1, po2);
        double f1x = Ex1.f(po1, x);
        double f2x = Ex1.f(po2, x);
        double f12x = Ex1.f(po12, x);
        assertEquals(f1x + f2x, f12x, Ex1.EPS);
    }

    @Test
    /**
     * Tests that p1+p2+ (-1*p2) == p1
     */
    void testAdd() {
        double[] p12 = Ex1.add(po1, po2);
        double[] minus1 = {-1};
        double[] pp2 = Ex1.mul(po2, minus1);
        double[] p1 = Ex1.add(p12, pp2);
        assertTrue(Ex1.equals(p1, po1));
    }

    @Test
    public void testPoly() {
        double[] p1 = {2, 0, 3.1, -1.2};
        assertEquals("-1.2x^3 +3.1x^2 +2.0", Ex1.poly(p1));
        double[] p2 = {-1.1, 2.3, 3.1};
        assertEquals("3.1x^2 +2.3x -1.1", Ex1.poly(p2));
        double[] p3 = {5};
        assertEquals("5.0", Ex1.poly(p3));
        double[] p4 = {0, -1};
        assertEquals("-1.0x", Ex1.poly(p4));
        double[] p5 = {0, 0, 4};
        assertEquals("4.0x^2", Ex1.poly(p5));
        double[] p6 = {0, 0, 0, 0};
        assertEquals("0", Ex1.poly(p6));
        double[] p7 = null;
        assertEquals("0", Ex1.poly(p7));
        double[] p8 = {};
        assertEquals("0", Ex1.poly(p8));
    }

    @Test
    void testMekadem() {
        assertEquals(3, Ex1.Mekadem(new double[]{0, 1, 0, 2}));
        assertEquals(2, Ex1.Mekadem(new double[]{1, 5, 3, 0, 0}));
        assertEquals(1, Ex1.Mekadem(new double[]{2, 1, Ex1.EPS / 2}));
        assertEquals(0, Ex1.Mekadem(new double[]{7}));
        assertEquals(0, Ex1.Mekadem(new double[]{0, 0, 0}));
        assertEquals(2, Ex1.Mekadem(new double[]{0, 0, Ex1.EPS * 2}));
    }

    @Test
    void testRemoveZeros() {

        assertArrayEquals(new double[]{3, 8, 2, 9}, Ex1.removeUnnecessaryZeros(new double[]{3, 8, 2, 9, 0, 0, 0}));
        assertArrayEquals(new double[]{2, 4, 0, 2}, Ex1.removeUnnecessaryZeros(new double[]{2, 4, 0, 2, 0, 0, 0}));
        assertArrayEquals(new double[]{0}, Ex1.removeUnnecessaryZeros(new double[]{0, 0, 0}));
        assertArrayEquals(new double[]{0, 3}, Ex1.removeUnnecessaryZeros(new double[]{0, 3}));
        assertArrayEquals(new double[]{9}, Ex1.removeUnnecessaryZeros(new double[]{9}));
        assertArrayEquals(new double[]{0}, Ex1.removeUnnecessaryZeros(new double[]{0}));
    }

    @Test
    void testpolinomZero() {
        assertTrue(Ex1.polinomZero(new double[]{0, 0, 0, 0}));
        assertFalse(Ex1.polinomZero(new double[]{4, 2, 1, 4, 8, 8, 24, 0, 0, 0, 0}));
        assertFalse(Ex1.polinomZero(new double[]{0, 0, 0, 0, 0.55, 0, 0, 0, 0, 0}));
    }

    @Test
    /**
     * Tests that p1+p2 == p2+p1
     */
    void testAdd2() {
        double[] p12 = Ex1.add(po1, po2);
        double[] p21 = Ex1.add(po2, po1);
        assertTrue(Ex1.equals(p12, p21));
    }

    @Test
    /**
     * Tests that p1+0 == p1
     */
    void testAdd3() {
        double[] p1 = Ex1.add(po1, Ex1.ZERO);
        assertTrue(Ex1.equals(p1, po1));
    }

    @Test
    /**
     * Tests that p1*0 == 0
     */
    void testMul1() {
        double[] p1 = Ex1.mul(po1, Ex1.ZERO);
        assertTrue(Ex1.equals(p1, Ex1.ZERO));
    }

    @Test
    /**
     * Tests that p1*p2 == p2*p1
     */
    void testMul2() {
        double[] p12 = Ex1.mul(po1, po2);
        double[] p21 = Ex1.mul(po2, po1);
        assertTrue(Ex1.equals(p12, p21));
    }

    @Test
    /**
     * Tests that p1(x) * p2(x) = (p1*p2)(x),
     */
    void testMulDoubleArrayDoubleArray() {
        double[] xx = {0, 1, 2, 3, 4.1, -15.2222};
        double[] p12 = Ex1.mul(po1, po2);
        for (int i = 0; i < xx.length; i = i + 1) {
            double x = xx[i];
            double f1x = Ex1.f(po1, x);
            double f2x = Ex1.f(po2, x);
            double f12x = Ex1.f(p12, x);
            assertEquals(f12x, f1x * f2x, Ex1.EPS);
        }
    }

    @Test
    /**
     * Tests a simple derivative examples - till ZERO.
     */
    void testDerivativeArrayDoubleArray() {
        double[] p = {1, 2, 3}; // 3X^2+2x+1
        double[] pt = {2, 6}; // 6x+2
        double[] dp1 = Ex1.derivative(p); // 2x + 6
        double[] dp2 = Ex1.derivative(dp1); // 2
        double[] dp3 = Ex1.derivative(dp2); // 0
        double[] dp4 = Ex1.derivative(dp3); // 0
        assertTrue(Ex1.equals(dp1, pt));
        assertTrue(Ex1.equals(Ex1.ZERO, dp3));
        assertTrue(Ex1.equals(dp4, dp3));
    }

    @Test
    /**
     * Tests the parsing of a polynom in a String like form.
     */
    public void testFromString() {
        double[] p = {-1.1, 2.3, 3.1}; // 3.1X^2+ 2.3x -1.1
        String sp2 = "3.1x^2 +2.3x -1.1";
        String sp = Ex1.poly(p);
        double[] p1 = Ex1.getPolynomFromString(sp);
        double[] p2 = Ex1.getPolynomFromString(sp2);
        boolean isSame1 = Ex1.equals(p1, p);
        boolean isSame2 = Ex1.equals(p2, p);
        System.out.println("sp:  " + sp);
        System.out.println("p1:  " + Arrays.toString(p1));
        System.out.println("poly(p1): " + Ex1.poly(p1));
        if (!isSame1) {
            fail();
        }
        if (!isSame2) {
            fail();
        }
        assertEquals(sp, Ex1.poly(p1));
    }

    @Test
    public void testNormalize() {
        assertEquals("", Ex1.normalize(null));
        assertEquals("", Ex1.normalize(""));
        assertEquals("+3.1x^2+2.3x-1.1", Ex1.normalize("3.1x^2 +2.3x -1.1"));
        assertEquals("-2x", Ex1.normalize("-2x"));
    }

    @Test
    public void testExtractMaxPower() {
        assertEquals(2, Ex1.extractMaxPower("+3.1x^2+2.3x-1.1"));
        assertEquals(1, Ex1.extractMaxPower("+5x-3"));
        assertEquals(0, Ex1.extractMaxPower("+7"));
    }

    @Test
    public void testParseCoef() {
        int[] idx;

        idx = new int[]{0};
        assertEquals(3.1, Ex1.parseCoef("+3.1x", idx), Ex1.EPS);

        idx = new int[]{0};
        assertEquals(-2.5, Ex1.parseCoef("-2.5x", idx), Ex1.EPS);

        idx = new int[]{0};
        assertEquals(1, Ex1.parseCoef("+x", idx), Ex1.EPS);

        idx = new int[]{0};
        assertEquals(-1, Ex1.parseCoef("-x", idx), Ex1.EPS);
    }

    @Test
    /**
     * Tests the equality of pairs of arrays.
     */
    public void testEquals() {
        double[][] d1 = {{0}, {1}, {1, 2, 0, 0}};
        double[][] d2 = {Ex1.ZERO, {1 + Ex1.EPS / 2}, {1, 2}};
        double[][] xx = {{-2 * Ex1.EPS}, {1 + Ex1.EPS * 1.2}, {1, 2, Ex1.EPS / 2}};
        for (int i = 0; i < d1.length; i = i + 1) {
            assertTrue(Ex1.equals(d1[i], d2[i]));
        }
        for (int i = 0; i < d1.length; i = i + 1) {
            assertFalse(Ex1.equals(d1[i], xx[i]));
        }
    }

    @Test
    /**
     * Tests is the sameValue function is symmetric.
     */
    public void testSameValue2() {
        double x1 = -4, x2 = 0;
        double rs1 = Ex1.sameValue(po1, po2, x1, x2, Ex1.EPS);
        double rs2 = Ex1.sameValue(po2, po1, x1, x2, Ex1.EPS);
        assertEquals(rs1, rs2, Ex1.EPS);
    }

    @Test
    /**
     * Test the area function - it should be symmetric.
     */
    public void testArea() {
        double x1 = -4, x2 = 0;
        double a1 = Ex1.area(po1, po2, x1, x2, 100);
        double a2 = Ex1.area(po2, po1, x1, x2, 100);
        assertEquals(a1, a2, Ex1.EPS);
    }

    @Test
    /**
     * Test the area f1(x)=0, f2(x)=x;
     */
    public void testArea2() {
        double[] po_a = Ex1.ZERO;
        double[] po_b = {0, 1};
        double x1 = -1;
        double x2 = 2;
        double a1 = Ex1.area(po_a, po_b, x1, x2, 1);
        double a2 = Ex1.area(po_a, po_b, x1, x2, 2);
        double a3 = Ex1.area(po_a, po_b, x1, x2, 3);
        double a100 = Ex1.area(po_a, po_b, x1, x2, 100);
        double area = 2.5;
        assertEquals(a1, area, Ex1.EPS);
        assertEquals(a2, area, Ex1.EPS);
        assertEquals(a3, area, Ex1.EPS);
        assertEquals(a100, area, Ex1.EPS);
    }

    @Test
    public void testPolynomFromPoints() {
        double[] xx1 = {0, 2};
        double[] yy1 = {3, 7};
        double[] expected1 = {3, 2};
        double[] res1 = Ex1.PolynomFromPoints(xx1, yy1);
        assertTrue(Ex1.equals(expected1, res1));
        double[] xx2 = {0, 1, 2};
        double[] yy2 = {2, 4, 10};
        double[] expected2 = {2, 0, 2};
        double[] res2 = Ex1.PolynomFromPoints(xx2, yy2);
        assertTrue(Ex1.equals(expected2, res2));
        double[] xx3 = {0, 1, 2, 3};
        double[] yy3 = {1, 2, 3, 4};
        double[] res3 = Ex1.PolynomFromPoints(xx3, yy3);
        assertNull(res3);
        double[] xx4 = {0, 1};
        double[] yy4 = {4};
        double[] res4 = Ex1.PolynomFromPoints(xx4, yy4);
        assertNull(res4);
    }

    @Test
    /**
     * Test the area function.
     */
    public void testArea3() {
        double[] po_a = {2, 1, -0.7, -0.02, 0.02};
        double[] po_b = {6, 0.1, -0.2};
        double x1 = Ex1.sameValue(po_a, po_b, -10, -5, Ex1.EPS);
        double a1 = Ex1.area(po_a, po_b, x1, 6, 8);
        double area = 58.5658;
        assertEquals(a1, area, Ex1.EPS);
    }
}
