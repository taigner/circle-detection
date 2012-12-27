import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;

import java.awt.*;
import java.awt.image.BufferedImage;

/**
 * Hough transform for detecting circles in an image
 *
 * @author Tobias Aigner
 */
public class Hough_Transform implements PlugInFilter {
    private static final int RESULTING_CIRCLE_LINE_WIDTH = 2;
    private static final Color RESULTING_CIRCLES_COLOR = Color.RED;

    private static final String FOUND_CIRCLE_MESSAGE = "Found circle (x, y, radius): (%s, %s, %s)";
    // area which is searched for circles, 1.0 = 100%
    private static final double SEARCH_AREA_PERCENT = 1.0;

    private ImageProcessor original;

    private int[][][] accumulator;

    private int radiusMin;
    private int radiusMax;
    private int numberOfSearchResults;

    private int x_midpoint;
    private int y_midpoint;
    private int x_offset;
    private int y_offset;

    private void houghTransform(Rectangle boundingBox) {
        prepareWith(boundingBox);
        ImageProcessor edges = buildEdgeImageFrom(original);

        for (int u = radiusMax; u < x_midpoint - radiusMax; u++)
            for (int v = radiusMax; v < y_midpoint - radiusMax; v++)
                if (edges.getPixel(u + x_offset, v + y_offset) != 0)
                    for (int p = 0; p < (radiusMax - radiusMin); p++)
                        bresenhamCircleInAccumulator(u, v, p);
    }

    private void prepareWith(Rectangle boundingBox) {
        x_midpoint = (int) (boundingBox.width * SEARCH_AREA_PERCENT);
        y_midpoint = (int) (boundingBox.height * SEARCH_AREA_PERCENT);
        x_offset = (int) boundingBox.getMinX();
        y_offset = (int) boundingBox.getMinY();

        accumulator = new int[x_midpoint][y_midpoint][radiusMax - radiusMin];
    }

    private ImageProcessor buildEdgeImageFrom(ImageProcessor original) {
        ImageProcessor duplicate = original.duplicate();
        duplicate.smooth();

        CannyEdgeDetector edgeDetector = buildCannyEdgeDetector(duplicate.getBufferedImage());
        edgeDetector.process();

        return new ColorProcessor(edgeDetector.getEdgesImage()).convertToByte(false);

//        // edge detection with built-in imagej methods
//        ImageProcessor edges = original.duplicate().convertToByte(true);
//        edges.smooth();
//        // open operation
//        edges.erode();
//        edges.dilate();
//        edges.findEdges();
//        edges.autoThreshold();
//        return edges;
    }

    private CannyEdgeDetector buildCannyEdgeDetector(BufferedImage image) {
        CannyEdgeDetector edgeDetector = new CannyEdgeDetector();
        edgeDetector.setLowThreshold(0.8f);
        edgeDetector.setHighThreshold(2.5f);
        edgeDetector.setSourceImage(image);
        return edgeDetector;
    }

    private void findAndDrawCircles(int maximumNumberOfCircles) {
        int lastMax;

        int[][] results = new int[maximumNumberOfCircles][3];
        int[] circle = new int[3];

        // skip 1% from edges
        final int resultArea = (int) (x_midpoint * 0.01);

        for (int maxCircles = 0; maxCircles < maximumNumberOfCircles; maxCircles++) {
            lastMax = -1;
            for (int u = resultArea; u < (x_midpoint - resultArea); u++) {
                for (int v = resultArea; v < (y_midpoint - resultArea); v++) {
                    for (int p = 0; p < (radiusMax - radiusMin); p++) {
                        if (accumulator[u][v][p] > lastMax) {
                            circle[0] = u;
                            circle[1] = v;
                            circle[2] = p;
                            lastMax = accumulator[u][v][p];
                        }
                    }
                }
            }

            System.arraycopy(circle, 0, results[maxCircles], 0, 3);
            clearFoundAreaInAccumulator(results[maxCircles]);
        }

        drawFoundCircles(results);
    }

    private void drawFoundCircles(int[][] results) {
        for (int[] result : results) {
            System.out.println(String.format(FOUND_CIRCLE_MESSAGE, result[0] + x_offset, result[1] + y_offset, result[2] + radiusMin));

            drawCircle(new ImagePlus("results", original).getProcessor(), result[0] + x_offset, result[1] + y_offset, result[2] + radiusMin);
        }
    }

    private void clearFoundAreaInAccumulator(int[] results) {
        int x = results[0];
        int y = results[1];
        int radius = results[2] + radiusMax + 3;

        for (int u = (x - radius); u < (x + radius); u++)
            for (int v = (y - radius); v < (y + radius); v++)
                for (int p = 0; p < (radiusMax - radiusMin); p++)
                    if (u > 0 && v > 0 && p > 0)
                        accumulator[u][v][p] = 0;
    }

    private void incrementInAccumulator(int x, int y, int p) {
        if (x < 0 || y < 0 || p < 0)
            return;

        accumulator[x][y][p]++;
    }

    /**
     * Bresenham algorithm for drawing circles inside accumulator.
     * For more details see http://en.wikipedia.org/wiki/Midpoint_circle_algorithm
     *
     * @param x0     Midpoint for x
     * @param y0     Midpoint for y
     * @param radius Radius of the circle
     */
    void bresenhamCircleInAccumulator(int x0, int y0, int radius) {
        int p = radius;
        radius += radiusMin;

        int f = 1 - radius;
        int ddF_x = 0;
        int ddF_y = -2 * radius;
        int x = 0;
        int y = radius;

        incrementInAccumulator(x0, y0 + radius, p);
        incrementInAccumulator(x0, y0 - radius, p);
        incrementInAccumulator(x0 + radius, y0, p);
        incrementInAccumulator(x0 - radius, y0, p);

        while (x < y) {
            if (f >= 0) {
                y--;
                ddF_y += 2;
                f += ddF_y;
            }

            x++;
            ddF_x += 2;
            f += ddF_x + 1;

            incrementInAccumulator(x0 + x, y0 + y, p);
            incrementInAccumulator(x0 - x, y0 + y, p);
            incrementInAccumulator(x0 + x, y0 - y, p);
            incrementInAccumulator(x0 - x, y0 - y, p);
            incrementInAccumulator(x0 + y, y0 + x, p);
            incrementInAccumulator(x0 - y, y0 + x, p);
            incrementInAccumulator(x0 + y, y0 - x, p);
            incrementInAccumulator(x0 - y, y0 - x, p);
        }
    }

    @Override
    public int setup(String arg, ImagePlus imp) {
        if (arg.equals("about")) {
            showAbout();
            return DONE;
        }

        return DOES_8G | DOES_16 | DOES_32 | DOES_RGB;
    }

    public void showAbout() {
        IJ.showMessage("HoughTransform", "detecting circles in an image");
    }

    @Override
    public void run(ImageProcessor imageProcessor) {
        original = imageProcessor;

        if (displayConfigurationDialog()) {
            houghTransform(original.getRoi());
            findAndDrawCircles(numberOfSearchResults);
        }
    }

    private boolean displayConfigurationDialog() {
        GenericDialog dialog = new GenericDialog("Hough Transform");

        dialog.addNumericField("Minimum Radius", 10, 0);
        dialog.addNumericField("Maximum Radius", 20, 0);
        dialog.addNumericField("Maximum Circles", 2, 0);

        dialog.showDialog();
        if (dialog.wasCanceled())
            return false;

        radiusMin = (int) dialog.getNextNumber();
        radiusMax = (int) dialog.getNextNumber();
        numberOfSearchResults = (int) dialog.getNextNumber();

        return true;
    }

    private void drawCircle(ImageProcessor IP, double x0, double y0, double r) {
        IP.setLineWidth(RESULTING_CIRCLE_LINE_WIDTH);
        IP.setColor(RESULTING_CIRCLES_COLOR);

        int n = 400;
        double dtheta = (Math.PI * 2) / n;

        int ix = (int) Math.round(x0 + r);
        int iy = (int) Math.round(y0);

        IP.moveTo(ix, iy);

        for (int i = 1; i <= n; i++) {
            double theta = i * dtheta;
            double x = x0 + r * Math.cos(theta);
            double y = y0 + r * Math.sin(theta);
            ix = (int) Math.round(x);
            iy = (int) Math.round(y);
            IP.lineTo(ix, iy);
        }
    }

    public static void main(String[] args) {
        Class<?> clazz = Hough_Transform.class;
        String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
        System.setProperty("plugins.dir", url.substring(5, url.length() - clazz.getName().length() - 6));

        new ImageJ();
    }
}
