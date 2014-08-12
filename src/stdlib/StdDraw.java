package stdlib;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.FileDialog;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.MediaTracker;
import java.awt.RenderingHints;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.Arc2D;
import java.awt.geom.Ellipse2D;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.awt.image.DirectColorModel;
import java.awt.image.WritableRaster;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.LinkedList;
import java.util.TreeSet;
import javax.imageio.ImageIO;
import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.KeyStroke;

/*************************************************************************
 *  Compilation:  javac StdDraw.java
 *  Execution:    java StdDraw
 *
 *  Standard drawing library. This class provides a basic capability for
 *  creating drawings with your programs. It uses a simple graphics model that
 *  allows you to create drawings consisting of points, lines, and curves
 *  in a window on your computer and to save the drawings to a file.
 *
 *  Todo
 *  ----
 *    -  Add support for gradient fill, etc.
 *
 *  Remarks
 *  -------
 *    -  don't use AffineTransform for rescaling since it inverts
 *       images and strings
 *    -  careful using setFont in inner loop within an animation -
 *       it can cause flicker
 *
 *************************************************************************/

/**
 *  <i>Standard draw</i>. This class provides a basic capability for
 *  creating drawings with your programs. It uses a simple graphics model that
 *  allows you to create drawings consisting of points, lines, and curves
 *  in a window on your computer and to save the drawings to a file.
 *  <p>
 *  For additional documentation, see <a href="http://introcs.cs.princeton.edu/15inout">Section 1.5</a> of
 *  <i>Introduction to Programming in Java: An Interdisciplinary Approach</i> by Robert Sedgewick and Kevin Wayne.
 */
public final class StdDraw implements ActionListener, MouseListener, MouseMotionListener, KeyListener {

	// pre-defined colors
	public static final Color BLACK      = Color.BLACK;
	public static final Color BLUE       = Color.BLUE;
	public static final Color CYAN       = Color.CYAN;
	public static final Color DARK_GRAY  = Color.DARK_GRAY;
	public static final Color GRAY       = Color.GRAY;
	public static final Color GREEN      = Color.GREEN;
	public static final Color LIGHT_GRAY = Color.LIGHT_GRAY;
	public static final Color MAGENTA    = Color.MAGENTA;
	public static final Color ORANGE     = Color.ORANGE;
	public static final Color PINK       = Color.PINK;
	public static final Color RED        = Color.RED;
	public static final Color WHITE      = Color.WHITE;
	public static final Color YELLOW     = Color.YELLOW;

	/**
	 * Shade of blue used in Introduction to Programming in Java.
	 * It is Pantone 300U. The RGB values are approximately (9, 90, 266).
	 */
	public static final Color BOOK_BLUE       = new Color(  9,  90, 166);
	public static final Color BOOK_LIGHT_BLUE = new Color(103, 198, 243);

	/**
	 * Shade of red used in Algorithms 4th edition.
	 * It is Pantone 1805U. The RGB values are approximately (150, 35, 31).
	 */
	public static final Color BOOK_RED = new Color(150, 35, 31);

	// default colors
	private static final Color DEFAULT_PEN_COLOR   = StdDraw.BLACK;
	private static final Color DEFAULT_CLEAR_COLOR = StdDraw.WHITE;

	// current pen color
	private static Color penColor;

	// default canvas size is DEFAULT_SIZE-by-DEFAULT_SIZE
	private static final int DEFAULT_SIZE = 512;
	private static int width  = StdDraw.DEFAULT_SIZE;
	private static int height = StdDraw.DEFAULT_SIZE;

	// default pen radius
	private static final double DEFAULT_PEN_RADIUS = 0.002;

	// current pen radius
	private static double penRadius;

	// show we draw immediately or wait until next show?
	private static boolean defer = false;

	// boundary of drawing canvas, 5% border
	private static final double BORDER = 0.05;
	private static final double DEFAULT_XMIN = 0.0;
	private static final double DEFAULT_XMAX = 1.0;
	private static final double DEFAULT_YMIN = 0.0;
	private static final double DEFAULT_YMAX = 1.0;
	private static double xmin, ymin, xmax, ymax;

	// for synchronization
	private static Object mouseLock = new Object();
	private static Object keyLock = new Object();

	// default font
	private static final Font DEFAULT_FONT = new Font("SansSerif", Font.PLAIN, 16);

	// current font
	private static Font font;

	// double buffered graphics
	private static BufferedImage offscreenImage, onscreenImage;
	private static Graphics2D offscreen, onscreen;

	// singleton for callbacks: avoids generation of extra .class files
	private static StdDraw std = new StdDraw();

	// the frame for drawing to the screen
	private static JFrame frame;

	// mouse state
	private static boolean mousePressed = false;
	private static double mouseX = 0;
	private static double mouseY = 0;

	// queue of typed key characters
	private static LinkedList<Character> keysTyped = new LinkedList<Character>();

	// set of key codes currently pressed down
	private static TreeSet<Integer> keysDown = new TreeSet<Integer>();


	// not instantiable
	private StdDraw() { }


	// static initializer
	static { StdDraw.init(); }

	/**
	 * Set the window size to the default size 512-by-512 pixels.
	 */
	public static void setCanvasSize() {
		StdDraw.setCanvasSize(StdDraw.DEFAULT_SIZE, StdDraw.DEFAULT_SIZE);
	}

	/**
	 * Set the window size to w-by-h pixels.
	 *
	 * @param w the width as a number of pixels
	 * @param h the height as a number of pixels
	 * @throws a RunTimeException if the width or height is 0 or negative
	 */
	public static void setCanvasSize(final int w, final int h) {
		if (w < 1 || h < 1) throw new RuntimeException("width and height must be positive");
		StdDraw.width = w;
		StdDraw.height = h;
		StdDraw.init();
	}

	// init
	private static void init() {
		if (StdDraw.frame != null) StdDraw.frame.setVisible(false);
		StdDraw.frame = new JFrame();
		StdDraw.offscreenImage = new BufferedImage(StdDraw.width, StdDraw.height, BufferedImage.TYPE_INT_ARGB);
		StdDraw.onscreenImage  = new BufferedImage(StdDraw.width, StdDraw.height, BufferedImage.TYPE_INT_ARGB);
		StdDraw.offscreen = StdDraw.offscreenImage.createGraphics();
		StdDraw.onscreen  = StdDraw.onscreenImage.createGraphics();
		StdDraw.setXscale();
		StdDraw.setYscale();
		StdDraw.offscreen.setColor(StdDraw.DEFAULT_CLEAR_COLOR);
		StdDraw.offscreen.fillRect(0, 0, StdDraw.width, StdDraw.height);
		StdDraw.setPenColor();
		StdDraw.setPenRadius();
		StdDraw.setFont();
		StdDraw.clear();

		// add antialiasing
		final RenderingHints hints = new RenderingHints(RenderingHints.KEY_ANTIALIASING,
				RenderingHints.VALUE_ANTIALIAS_ON);
		hints.put(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
		StdDraw.offscreen.addRenderingHints(hints);

		// frame stuff
		final ImageIcon icon = new ImageIcon(StdDraw.onscreenImage);
		final JLabel draw = new JLabel(icon);

		draw.addMouseListener(StdDraw.std);
		draw.addMouseMotionListener(StdDraw.std);

		StdDraw.frame.setContentPane(draw);
		StdDraw.frame.addKeyListener(StdDraw.std);    // JLabel cannot get keyboard focus
		StdDraw.frame.setResizable(false);
		StdDraw.frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);            // closes all windows
		// frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);      // closes only current window
		StdDraw.frame.setTitle("Standard Draw");
		StdDraw.frame.setJMenuBar(StdDraw.createMenuBar());
		StdDraw.frame.pack();
		StdDraw.frame.requestFocusInWindow();
		StdDraw.frame.setVisible(true);
	}

	// create the menu bar (changed to private)
	private static JMenuBar createMenuBar() {
		final JMenuBar menuBar = new JMenuBar();
		final JMenu menu = new JMenu("File");
		menuBar.add(menu);
		final JMenuItem menuItem1 = new JMenuItem(" Save...   ");
		menuItem1.addActionListener(StdDraw.std);
		menuItem1.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_S,
				Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
		menu.add(menuItem1);
		return menuBar;
	}


	/*************************************************************************
	 *  User and screen coordinate systems
	 *************************************************************************/

	/**
	 * Set the x-scale to be the default (between 0.0 and 1.0).
	 */
	public static void setXscale() { StdDraw.setXscale(StdDraw.DEFAULT_XMIN, StdDraw.DEFAULT_XMAX); }

	/**
	 * Set the y-scale to be the default (between 0.0 and 1.0).
	 */
	public static void setYscale() { StdDraw.setYscale(StdDraw.DEFAULT_YMIN, StdDraw.DEFAULT_YMAX); }

	/**
	 * Set the x-scale (a 10% border is added to the values)
	 * @param min the minimum value of the x-scale
	 * @param max the maximum value of the x-scale
	 */
	public static void setXscale(final double min, final double max) {
		final double size = max - min;
		StdDraw.xmin = min - StdDraw.BORDER * size;
		StdDraw.xmax = max + StdDraw.BORDER * size;
	}

	/**
	 * Set the y-scale (a 10% border is added to the values).
	 * @param min the minimum value of the y-scale
	 * @param max the maximum value of the y-scale
	 */
	public static void setYscale(final double min, final double max) {
		final double size = max - min;
		StdDraw.ymin = min - StdDraw.BORDER * size;
		StdDraw.ymax = max + StdDraw.BORDER * size;
	}

	/**
	 * Set the x-scale and y-scale (a 10% border is added to the values)
	 * @param min the minimum value of the x- and y-scales
	 * @param max the maximum value of the x- and y-scales
	 */
	public static void setScale(final double min, final double max) {
		StdDraw.setXscale(min, max);
		StdDraw.setYscale(min, max);
	}

	// helper functions that scale from user coordinates to screen coordinates and back
	private static double  scaleX(final double x) { return StdDraw.width  * (x - StdDraw.xmin) / (StdDraw.xmax - StdDraw.xmin); }
	private static double  scaleY(final double y) { return StdDraw.height * (StdDraw.ymax - y) / (StdDraw.ymax - StdDraw.ymin); }
	private static double factorX(final double w) { return w * StdDraw.width  / Math.abs(StdDraw.xmax - StdDraw.xmin);  }
	private static double factorY(final double h) { return h * StdDraw.height / Math.abs(StdDraw.ymax - StdDraw.ymin);  }
	private static double   userX(final double x) { return StdDraw.xmin + x * (StdDraw.xmax - StdDraw.xmin) / StdDraw.width;    }
	private static double   userY(final double y) { return StdDraw.ymax - y * (StdDraw.ymax - StdDraw.ymin) / StdDraw.height;   }


	/**
	 * Clear the screen to the default color (white).
	 */
	public static void clear() { StdDraw.clear(StdDraw.DEFAULT_CLEAR_COLOR); }
	/**
	 * Clear the screen to the given color.
	 * @param color the Color to make the background
	 */
	public static void clear(final Color color) {
		StdDraw.offscreen.setColor(color);
		StdDraw.offscreen.fillRect(0, 0, StdDraw.width, StdDraw.height);
		StdDraw.offscreen.setColor(StdDraw.penColor);
		StdDraw.draw();
	}

	/**
	 * Get the current pen radius.
	 */
	public static double getPenRadius() { return StdDraw.penRadius; }

	/**
	 * Set the pen size to the default (.002).
	 */
	public static void setPenRadius() { StdDraw.setPenRadius(StdDraw.DEFAULT_PEN_RADIUS); }
	/**
	 * Set the radius of the pen to the given size.
	 * @param r the radius of the pen
	 * @throws RuntimeException if r is negative
	 */
	public static void setPenRadius(final double r) {
		if (r < 0) throw new RuntimeException("pen radius must be positive");
		StdDraw.penRadius = r * StdDraw.DEFAULT_SIZE;
		final BasicStroke stroke = new BasicStroke((float) StdDraw.penRadius, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);
		// BasicStroke stroke = new BasicStroke((float) penRadius);
		StdDraw.offscreen.setStroke(stroke);
	}

	/**
	 * Get the current pen color.
	 */
	public static Color getPenColor() { return StdDraw.penColor; }

	/**
	 * Set the pen color to the default color (black).
	 */
	public static void setPenColor() { StdDraw.setPenColor(StdDraw.DEFAULT_PEN_COLOR); }
	/**
	 * Set the pen color to the given color. The available pen colors are
	 * BLACK, BLUE, CYAN, DARK_GRAY, GRAY, GREEN, LIGHT_GRAY, MAGENTA,
	 * ORANGE, PINK, RED, WHITE, and YELLOW.
	 * @param color the Color to make the pen
	 */
	public static void setPenColor(final Color color) {
		StdDraw.penColor = color;
		StdDraw.offscreen.setColor(StdDraw.penColor);
	}

	/**
	 * Get the current font.
	 */
	public static Font getFont() { return StdDraw.font; }

	/**
	 * Set the font to the default font (sans serif, 16 point).
	 */
	public static void setFont() { StdDraw.setFont(StdDraw.DEFAULT_FONT); }

	/**
	 * Set the font to the given value.
	 * @param f the font to make text
	 */
	public static void setFont(final Font f) { StdDraw.font = f; }


	/*************************************************************************
	 *  Drawing geometric shapes.
	 *************************************************************************/

	/**
	 * Draw a line from (x0, y0) to (x1, y1).
	 * @param x0 the x-coordinate of the starting point
	 * @param y0 the y-coordinate of the starting point
	 * @param x1 the x-coordinate of the destination point
	 * @param y1 the y-coordinate of the destination point
	 */
	public static void line(final double x0, final double y0, final double x1, final double y1) {
		StdDraw.offscreen.draw(new Line2D.Double(StdDraw.scaleX(x0), StdDraw.scaleY(y0), StdDraw.scaleX(x1), StdDraw.scaleY(y1)));
		StdDraw.draw();
	}

	/**
	 * Draw one pixel at (x, y).
	 * @param x the x-coordinate of the pixel
	 * @param y the y-coordinate of the pixel
	 */
	private static void pixel(final double x, final double y) {
		StdDraw.offscreen.fillRect((int) Math.round(StdDraw.scaleX(x)), (int) Math.round(StdDraw.scaleY(y)), 1, 1);
	}

	/**
	 * Draw a point at (x, y).
	 * @param x the x-coordinate of the point
	 * @param y the y-coordinate of the point
	 */
	public static void point(final double x, final double y) {
		final double xs = StdDraw.scaleX(x);
		final double ys = StdDraw.scaleY(y);
		final double r = StdDraw.penRadius;
		// double ws = factorX(2*r);
		// double hs = factorY(2*r);
		// if (ws <= 1 && hs <= 1) pixel(x, y);
		if (r <= 1) StdDraw.pixel(x, y);
		else StdDraw.offscreen.fill(new Ellipse2D.Double(xs - r/2, ys - r/2, r, r));
		StdDraw.draw();
	}

	/**
	 * Draw a circle of radius r, centered on (x, y).
	 * @param x the x-coordinate of the center of the circle
	 * @param y the y-coordinate of the center of the circle
	 * @param r the radius of the circle
	 * @throws RuntimeException if the radius of the circle is negative
	 */
	public static void circle(final double x, final double y, final double r) {
		if (r < 0) throw new RuntimeException("circle radius can't be negative");
		final double xs = StdDraw.scaleX(x);
		final double ys = StdDraw.scaleY(y);
		final double ws = StdDraw.factorX(2*r);
		final double hs = StdDraw.factorY(2*r);
		if (ws <= 1 && hs <= 1) StdDraw.pixel(x, y);
		else StdDraw.offscreen.draw(new Ellipse2D.Double(xs - ws/2, ys - hs/2, ws, hs));
		StdDraw.draw();
	}

	/**
	 * Draw filled circle of radius r, centered on (x, y).
	 * @param x the x-coordinate of the center of the circle
	 * @param y the y-coordinate of the center of the circle
	 * @param r the radius of the circle
	 * @throws RuntimeException if the radius of the circle is negative
	 */
	public static void filledCircle(final double x, final double y, final double r) {
		if (r < 0) throw new RuntimeException("circle radius can't be negative");
		final double xs = StdDraw.scaleX(x);
		final double ys = StdDraw.scaleY(y);
		final double ws = StdDraw.factorX(2*r);
		final double hs = StdDraw.factorY(2*r);
		if (ws <= 1 && hs <= 1) StdDraw.pixel(x, y);
		else StdDraw.offscreen.fill(new Ellipse2D.Double(xs - ws/2, ys - hs/2, ws, hs));
		StdDraw.draw();
	}


	/**
	 * Draw an ellipse with given semimajor and semiminor axes, centered on (x, y).
	 * @param x the x-coordinate of the center of the ellipse
	 * @param y the y-coordinate of the center of the ellipse
	 * @param semiMajorAxis is the semimajor axis of the ellipse
	 * @param semiMinorAxis is the semiminor axis of the ellipse
	 * @throws RuntimeException if either of the axes are negative
	 */
	public static void ellipse(final double x, final double y, final double semiMajorAxis, final double semiMinorAxis) {
		if (semiMajorAxis < 0) throw new RuntimeException("ellipse semimajor axis can't be negative");
		if (semiMinorAxis < 0) throw new RuntimeException("ellipse semiminor axis can't be negative");
		final double xs = StdDraw.scaleX(x);
		final double ys = StdDraw.scaleY(y);
		final double ws = StdDraw.factorX(2*semiMajorAxis);
		final double hs = StdDraw.factorY(2*semiMinorAxis);
		if (ws <= 1 && hs <= 1) StdDraw.pixel(x, y);
		else StdDraw.offscreen.draw(new Ellipse2D.Double(xs - ws/2, ys - hs/2, ws, hs));
		StdDraw.draw();
	}

	/**
	 * Draw an ellipse with given semimajor and semiminor axes, centered on (x, y).
	 * @param x the x-coordinate of the center of the ellipse
	 * @param y the y-coordinate of the center of the ellipse
	 * @param semiMajorAxis is the semimajor axis of the ellipse
	 * @param semiMinorAxis is the semiminor axis of the ellipse
	 * @throws RuntimeException if either of the axes are negative
	 */
	public static void filledEllipse(final double x, final double y, final double semiMajorAxis, final double semiMinorAxis) {
		if (semiMajorAxis < 0) throw new RuntimeException("ellipse semimajor axis can't be negative");
		if (semiMinorAxis < 0) throw new RuntimeException("ellipse semiminor axis can't be negative");
		final double xs = StdDraw.scaleX(x);
		final double ys = StdDraw.scaleY(y);
		final double ws = StdDraw.factorX(2*semiMajorAxis);
		final double hs = StdDraw.factorY(2*semiMinorAxis);
		if (ws <= 1 && hs <= 1) StdDraw.pixel(x, y);
		else StdDraw.offscreen.fill(new Ellipse2D.Double(xs - ws/2, ys - hs/2, ws, hs));
		StdDraw.draw();
	}


	/**
	 * Draw an arc of radius r, centered on (x, y), from angle1 to angle2 (in degrees).
	 * @param x the x-coordinate of the center of the circle
	 * @param y the y-coordinate of the center of the circle
	 * @param r the radius of the circle
	 * @param angle1 the starting angle. 0 would mean an arc beginning at 3 o'clock.
	 * @param angle2 the angle at the end of the arc. For example, if
	 *        you want a 90 degree arc, then angle2 should be angle1 + 90.
	 * @throws RuntimeException if the radius of the circle is negative
	 */
	public static void arc(final double x, final double y, final double r, final double angle1, double angle2) {
		if (r < 0) throw new RuntimeException("arc radius can't be negative");
		while (angle2 < angle1) angle2 += 360;
		final double xs = StdDraw.scaleX(x);
		final double ys = StdDraw.scaleY(y);
		final double ws = StdDraw.factorX(2*r);
		final double hs = StdDraw.factorY(2*r);
		if (ws <= 1 && hs <= 1) StdDraw.pixel(x, y);
		else StdDraw.offscreen.draw(new Arc2D.Double(xs - ws/2, ys - hs/2, ws, hs, angle1, angle2 - angle1, Arc2D.OPEN));
		StdDraw.draw();
	}

	/**
	 * Draw a square of side length 2r, centered on (x, y).
	 * @param x the x-coordinate of the center of the square
	 * @param y the y-coordinate of the center of the square
	 * @param r radius is half the length of any side of the square
	 * @throws RuntimeException if r is negative
	 */
	public static void square(final double x, final double y, final double r) {
		if (r < 0) throw new RuntimeException("square side length can't be negative");
		final double xs = StdDraw.scaleX(x);
		final double ys = StdDraw.scaleY(y);
		final double ws = StdDraw.factorX(2*r);
		final double hs = StdDraw.factorY(2*r);
		if (ws <= 1 && hs <= 1) StdDraw.pixel(x, y);
		else StdDraw.offscreen.draw(new Rectangle2D.Double(xs - ws/2, ys - hs/2, ws, hs));
		StdDraw.draw();
	}

	/**
	 * Draw a filled square of side length 2r, centered on (x, y).
	 * @param x the x-coordinate of the center of the square
	 * @param y the y-coordinate of the center of the square
	 * @param r radius is half the length of any side of the square
	 * @throws RuntimeException if r is negative
	 */
	public static void filledSquare(final double x, final double y, final double r) {
		if (r < 0) throw new RuntimeException("square side length can't be negative");
		final double xs = StdDraw.scaleX(x);
		final double ys = StdDraw.scaleY(y);
		final double ws = StdDraw.factorX(2*r);
		final double hs = StdDraw.factorY(2*r);
		if (ws <= 1 && hs <= 1) StdDraw.pixel(x, y);
		else StdDraw.offscreen.fill(new Rectangle2D.Double(xs - ws/2, ys - hs/2, ws, hs));
		StdDraw.draw();
	}


	/**
	 * Draw a rectangle of given half width and half height, centered on (x, y).
	 * @param x the x-coordinate of the center of the rectangle
	 * @param y the y-coordinate of the center of the rectangle
	 * @param halfWidth is half the width of the rectangle
	 * @param halfHeight is half the height of the rectangle
	 * @throws RuntimeException if halfWidth or halfHeight is negative
	 */
	public static void rectangle(final double x, final double y, final double halfWidth, final double halfHeight) {
		if (halfWidth  < 0) throw new RuntimeException("half width can't be negative");
		if (halfHeight < 0) throw new RuntimeException("half height can't be negative");
		final double xs = StdDraw.scaleX(x);
		final double ys = StdDraw.scaleY(y);
		final double ws = StdDraw.factorX(2*halfWidth);
		final double hs = StdDraw.factorY(2*halfHeight);
		if (ws <= 1 && hs <= 1) StdDraw.pixel(x, y);
		else StdDraw.offscreen.draw(new Rectangle2D.Double(xs - ws/2, ys - hs/2, ws, hs));
		StdDraw.draw();
	}

	/**
	 * Draw a filled rectangle of given half width and half height, centered on (x, y).
	 * @param x the x-coordinate of the center of the rectangle
	 * @param y the y-coordinate of the center of the rectangle
	 * @param halfWidth is half the width of the rectangle
	 * @param halfHeight is half the height of the rectangle
	 * @throws RuntimeException if halfWidth or halfHeight is negative
	 */
	public static void filledRectangle(final double x, final double y, final double halfWidth, final double halfHeight) {
		if (halfWidth  < 0) throw new RuntimeException("half width can't be negative");
		if (halfHeight < 0) throw new RuntimeException("half height can't be negative");
		final double xs = StdDraw.scaleX(x);
		final double ys = StdDraw.scaleY(y);
		final double ws = StdDraw.factorX(2*halfWidth);
		final double hs = StdDraw.factorY(2*halfHeight);
		if (ws <= 1 && hs <= 1) StdDraw.pixel(x, y);
		else StdDraw.offscreen.fill(new Rectangle2D.Double(xs - ws/2, ys - hs/2, ws, hs));
		StdDraw.draw();
	}


	/**
	 * Draw a polygon with the given (x[i], y[i]) coordinates.
	 * @param x an array of all the x-coordindates of the polygon
	 * @param y an array of all the y-coordindates of the polygon
	 */
	public static void polygon(final double[] x, final double[] y) {
		final int N = x.length;
		final GeneralPath path = new GeneralPath();
		path.moveTo((float) StdDraw.scaleX(x[0]), (float) StdDraw.scaleY(y[0]));
		for (int i = 0; i < N; i++)
			path.lineTo((float) StdDraw.scaleX(x[i]), (float) StdDraw.scaleY(y[i]));
		path.closePath();
		StdDraw.offscreen.draw(path);
		StdDraw.draw();
	}

	/**
	 * Draw a filled polygon with the given (x[i], y[i]) coordinates.
	 * @param x an array of all the x-coordindates of the polygon
	 * @param y an array of all the y-coordindates of the polygon
	 */
	public static void filledPolygon(final double[] x, final double[] y) {
		final int N = x.length;
		final GeneralPath path = new GeneralPath();
		path.moveTo((float) StdDraw.scaleX(x[0]), (float) StdDraw.scaleY(y[0]));
		for (int i = 0; i < N; i++)
			path.lineTo((float) StdDraw.scaleX(x[i]), (float) StdDraw.scaleY(y[i]));
		path.closePath();
		StdDraw.offscreen.fill(path);
		StdDraw.draw();
	}



	/*************************************************************************
	 *  Drawing images.
	 *************************************************************************/

	// get an image from the given filename
	private static Image getImage(final String filename) {

		// to read from file
		ImageIcon icon = new ImageIcon(filename);

		// try to read from URL
		if ((icon == null) || (icon.getImageLoadStatus() != MediaTracker.COMPLETE)) {
			try {
				final URL url = new URL(filename);
				icon = new ImageIcon(url);
			} catch (final Exception e) { /* not a url */ }
		}

		// in case file is inside a .jar
		if ((icon == null) || (icon.getImageLoadStatus() != MediaTracker.COMPLETE)) {
			final URL url = StdDraw.class.getResource(filename);
			if (url == null) throw new RuntimeException("image " + filename + " not found");
			icon = new ImageIcon(url);
		}

		return icon.getImage();
	}

	/**
	 * Draw picture (gif, jpg, or png) centered on (x, y).
	 * @param x the center x-coordinate of the image
	 * @param y the center y-coordinate of the image
	 * @param s the name of the image/picture, e.g., "ball.gif"
	 * @throws RuntimeException if the image is corrupt
	 */
	public static void picture(final double x, final double y, final String s) {
		final Image image = StdDraw.getImage(s);
		final double xs = StdDraw.scaleX(x);
		final double ys = StdDraw.scaleY(y);
		final int ws = image.getWidth(null);
		final int hs = image.getHeight(null);
		if (ws < 0 || hs < 0) throw new RuntimeException("image " + s + " is corrupt");

		StdDraw.offscreen.drawImage(image, (int) Math.round(xs - ws/2.0), (int) Math.round(ys - hs/2.0), null);
		StdDraw.draw();
	}

	/**
	 * Draw picture (gif, jpg, or png) centered on (x, y),
	 * rotated given number of degrees
	 * @param x the center x-coordinate of the image
	 * @param y the center y-coordinate of the image
	 * @param s the name of the image/picture, e.g., "ball.gif"
	 * @param degrees is the number of degrees to rotate counterclockwise
	 * @throws RuntimeException if the image is corrupt
	 */
	public static void picture(final double x, final double y, final String s, final double degrees) {
		final Image image = StdDraw.getImage(s);
		final double xs = StdDraw.scaleX(x);
		final double ys = StdDraw.scaleY(y);
		final int ws = image.getWidth(null);
		final int hs = image.getHeight(null);
		if (ws < 0 || hs < 0) throw new RuntimeException("image " + s + " is corrupt");

		StdDraw.offscreen.rotate(Math.toRadians(-degrees), xs, ys);
		StdDraw.offscreen.drawImage(image, (int) Math.round(xs - ws/2.0), (int) Math.round(ys - hs/2.0), null);
		StdDraw.offscreen.rotate(Math.toRadians(+degrees), xs, ys);

		StdDraw.draw();
	}

	/**
	 * Draw picture (gif, jpg, or png) centered on (x, y), rescaled to w-by-h.
	 * @param x the center x coordinate of the image
	 * @param y the center y coordinate of the image
	 * @param s the name of the image/picture, e.g., "ball.gif"
	 * @param w the width of the image
	 * @param h the height of the image
	 * @throws RuntimeException if the width height are negative
	 * @throws RuntimeException if the image is corrupt
	 */
	public static void picture(final double x, final double y, final String s, final double w, final double h) {
		final Image image = StdDraw.getImage(s);
		final double xs = StdDraw.scaleX(x);
		final double ys = StdDraw.scaleY(y);
		if (w < 0) throw new RuntimeException("width is negative: " + w);
		if (h < 0) throw new RuntimeException("height is negative: " + h);
		final double ws = StdDraw.factorX(w);
		final double hs = StdDraw.factorY(h);
		if (ws < 0 || hs < 0) throw new RuntimeException("image " + s + " is corrupt");
		if (ws <= 1 && hs <= 1) StdDraw.pixel(x, y);
		else {
			StdDraw.offscreen.drawImage(image, (int) Math.round(xs - ws/2.0),
					(int) Math.round(ys - hs/2.0),
					(int) Math.round(ws),
					(int) Math.round(hs), null);
		}
		StdDraw.draw();
	}


	/**
	 * Draw picture (gif, jpg, or png) centered on (x, y), rotated
	 * given number of degrees, rescaled to w-by-h.
	 * @param x the center x-coordinate of the image
	 * @param y the center y-coordinate of the image
	 * @param s the name of the image/picture, e.g., "ball.gif"
	 * @param w the width of the image
	 * @param h the height of the image
	 * @param degrees is the number of degrees to rotate counterclockwise
	 * @throws RuntimeException if the image is corrupt
	 */
	public static void picture(final double x, final double y, final String s, final double w, final double h, final double degrees) {
		final Image image = StdDraw.getImage(s);
		final double xs = StdDraw.scaleX(x);
		final double ys = StdDraw.scaleY(y);
		final double ws = StdDraw.factorX(w);
		final double hs = StdDraw.factorY(h);
		if (ws < 0 || hs < 0) throw new RuntimeException("image " + s + " is corrupt");
		if (ws <= 1 && hs <= 1) StdDraw.pixel(x, y);

		StdDraw.offscreen.rotate(Math.toRadians(-degrees), xs, ys);
		StdDraw.offscreen.drawImage(image, (int) Math.round(xs - ws/2.0),
				(int) Math.round(ys - hs/2.0),
				(int) Math.round(ws),
				(int) Math.round(hs), null);
		StdDraw.offscreen.rotate(Math.toRadians(+degrees), xs, ys);

		StdDraw.draw();
	}


	/*************************************************************************
	 *  Drawing text.
	 *************************************************************************/

	/**
	 * Write the given text string in the current font, centered on (x, y).
	 * @param x the center x-coordinate of the text
	 * @param y the center y-coordinate of the text
	 * @param s the text
	 */
	public static void text(final double x, final double y, final String s) {
		StdDraw.offscreen.setFont(StdDraw.font);
		final FontMetrics metrics = StdDraw.offscreen.getFontMetrics();
		final double xs = StdDraw.scaleX(x);
		final double ys = StdDraw.scaleY(y);
		final int ws = metrics.stringWidth(s);
		final int hs = metrics.getDescent();
		StdDraw.offscreen.drawString(s, (float) (xs - ws/2.0), (float) (ys + hs));
		StdDraw.draw();
	}

	/**
	 * Write the given text string in the current font, centered on (x, y) and
	 * rotated by the specified number of degrees
	 * @param x the center x-coordinate of the text
	 * @param y the center y-coordinate of the text
	 * @param s the text
	 * @param degrees is the number of degrees to rotate counterclockwise
	 */
	public static void text(final double x, final double y, final String s, final double degrees) {
		final double xs = StdDraw.scaleX(x);
		final double ys = StdDraw.scaleY(y);
		StdDraw.offscreen.rotate(Math.toRadians(-degrees), xs, ys);
		StdDraw.text(x, y, s);
		StdDraw.offscreen.rotate(Math.toRadians(+degrees), xs, ys);
	}


	/**
	 * Write the given text string in the current font, left-aligned at (x, y).
	 * @param x the x-coordinate of the text
	 * @param y the y-coordinate of the text
	 * @param s the text
	 */
	public static void textLeft(final double x, final double y, final String s) {
		StdDraw.offscreen.setFont(StdDraw.font);
		final FontMetrics metrics = StdDraw.offscreen.getFontMetrics();
		final double xs = StdDraw.scaleX(x);
		final double ys = StdDraw.scaleY(y);
		final int hs = metrics.getDescent();
		StdDraw.offscreen.drawString(s, (float) (xs), (float) (ys + hs));
		StdDraw.show();
	}

	/**
	 * Write the given text string in the current font, right-aligned at (x, y).
	 * @param x the x-coordinate of the text
	 * @param y the y-coordinate of the text
	 * @param s the text
	 */
	public static void textRight(final double x, final double y, final String s) {
		StdDraw.offscreen.setFont(StdDraw.font);
		final FontMetrics metrics = StdDraw.offscreen.getFontMetrics();
		final double xs = StdDraw.scaleX(x);
		final double ys = StdDraw.scaleY(y);
		final int ws = metrics.stringWidth(s);
		final int hs = metrics.getDescent();
		StdDraw.offscreen.drawString(s, (float) (xs - ws), (float) (ys + hs));
		StdDraw.show();
	}



	/**
	 * Display on screen, pause for t milliseconds, and turn on
	 * <em>animation mode</em>: subsequent calls to
	 * drawing methods such as <tt>line()</tt>, <tt>circle()</tt>, and <tt>square()</tt>
	 * will not be displayed on screen until the next call to <tt>show()</tt>.
	 * This is useful for producing animations (clear the screen, draw a bunch of shapes,
	 * display on screen for a fixed amount of time, and repeat). It also speeds up
	 * drawing a huge number of shapes (call <tt>show(0)</tt> to defer drawing
	 * on screen, draw the shapes, and call <tt>show(0)</tt> to display them all
	 * on screen at once).
	 * @param t number of milliseconds
	 */
	public static void show(final int t) {
		StdDraw.defer = false;
		StdDraw.draw();
		try { Thread.sleep(t); }
		catch (final InterruptedException e) { System.out.println("Error sleeping"); }
		StdDraw.defer = true;
	}


	/**
	 * Display on-screen and turn off animation mode:
	 * subsequent calls to
	 * drawing methods such as <tt>line()</tt>, <tt>circle()</tt>, and <tt>square()</tt>
	 * will be displayed on screen when called. This is the default.
	 */
	public static void show() {
		StdDraw.defer = false;
		StdDraw.draw();
	}

	// draw onscreen if defer is false
	private static void draw() {
		if (StdDraw.defer) return;
		StdDraw.onscreen.drawImage(StdDraw.offscreenImage, 0, 0, null);
		StdDraw.frame.repaint();
	}


	/*************************************************************************
	 *  Save drawing to a file.
	 *************************************************************************/

	/**
	 * Save onscreen image to file - suffix must be png, jpg, or gif.
	 * @param filename the name of the file with one of the required suffixes
	 */
	public static void save(final String filename) {
		final File file = new File(filename);
		final String suffix = filename.substring(filename.lastIndexOf('.') + 1);

		// png files
		if (suffix.toLowerCase().equals("png")) {
			try { ImageIO.write(StdDraw.offscreenImage, suffix, file); }
			catch (final IOException e) { e.printStackTrace(); }
		}

		// need to change from ARGB to RGB for jpeg
		// reference: http://archives.java.sun.com/cgi-bin/wa?A2=ind0404&L=java2d-interest&D=0&P=2727
		else if (suffix.toLowerCase().equals("jpg")) {
			final WritableRaster raster = StdDraw.onscreenImage.getRaster();
			WritableRaster newRaster;
			newRaster = raster.createWritableChild(0, 0, StdDraw.width, StdDraw.height, 0, 0, new int[] {0, 1, 2});
			final DirectColorModel cm = (DirectColorModel) StdDraw.onscreenImage.getColorModel();
			final DirectColorModel newCM = new DirectColorModel(cm.getPixelSize(),
					cm.getRedMask(),
					cm.getGreenMask(),
					cm.getBlueMask());
			final BufferedImage rgbBuffer = new BufferedImage(newCM, newRaster, false,  null);
			try { ImageIO.write(rgbBuffer, suffix, file); }
			catch (final IOException e) { e.printStackTrace(); }
		}

		else {
			System.out.println("Invalid image file type: " + suffix);
		}
	}


	/**
	 * This method cannot be called directly.
	 */
	public void actionPerformed(final ActionEvent e) {
		final FileDialog chooser = new FileDialog(StdDraw.frame, "Use a .png or .jpg extension", FileDialog.SAVE);
		chooser.setVisible(true);
		final String filename = chooser.getFile();
		if (filename != null) {
			StdDraw.save(chooser.getDirectory() + File.separator + chooser.getFile());
		}
	}


	/*************************************************************************
	 *  Mouse interactions.
	 *************************************************************************/

	/**
	 * Is the mouse being pressed?
	 * @return true or false
	 */
	public static boolean mousePressed() {
		synchronized (StdDraw.mouseLock) {
			return StdDraw.mousePressed;
		}
	}

	/**
	 * What is the x-coordinate of the mouse?
	 * @return the value of the x-coordinate of the mouse
	 */
	public static double mouseX() {
		synchronized (StdDraw.mouseLock) {
			return StdDraw.mouseX;
		}
	}

	/**
	 * What is the y-coordinate of the mouse?
	 * @return the value of the y-coordinate of the mouse
	 */
	public static double mouseY() {
		synchronized (StdDraw.mouseLock) {
			return StdDraw.mouseY;
		}
	}


	/**
	 * This method cannot be called directly.
	 */
	public void mouseClicked(final MouseEvent e) { }

	/**
	 * This method cannot be called directly.
	 */
	public void mouseEntered(final MouseEvent e) { }

	/**
	 * This method cannot be called directly.
	 */
	public void mouseExited(final MouseEvent e) { }

	/**
	 * This method cannot be called directly.
	 */
	public void mousePressed(final MouseEvent e) {
		synchronized (StdDraw.mouseLock) {
			StdDraw.mouseX = StdDraw.userX(e.getX());
			StdDraw.mouseY = StdDraw.userY(e.getY());
			StdDraw.mousePressed = true;
		}
	}

	/**
	 * This method cannot be called directly.
	 */
	public void mouseReleased(final MouseEvent e) {
		synchronized (StdDraw.mouseLock) {
			StdDraw.mousePressed = false;
		}
	}

	/**
	 * This method cannot be called directly.
	 */
	public void mouseDragged(final MouseEvent e)  {
		synchronized (StdDraw.mouseLock) {
			StdDraw.mouseX = StdDraw.userX(e.getX());
			StdDraw.mouseY = StdDraw.userY(e.getY());
		}
	}

	/**
	 * This method cannot be called directly.
	 */
	public void mouseMoved(final MouseEvent e) {
		synchronized (StdDraw.mouseLock) {
			StdDraw.mouseX = StdDraw.userX(e.getX());
			StdDraw.mouseY = StdDraw.userY(e.getY());
		}
	}


	/*************************************************************************
	 *  Keyboard interactions.
	 *************************************************************************/

	/**
	 * Has the user typed a key?
	 * @return true if the user has typed a key, false otherwise
	 */
	public static boolean hasNextKeyTyped() {
		synchronized (StdDraw.keyLock) {
			return !StdDraw.keysTyped.isEmpty();
		}
	}

	/**
	 * What is the next key that was typed by the user? This method returns
	 * a Unicode character corresponding to the key typed (such as 'a' or 'A').
	 * It cannot identify action keys (such as F1
	 * and arrow keys) or modifier keys (such as control).
	 * @return the next Unicode key typed
	 */
	public static char nextKeyTyped() {
		synchronized (StdDraw.keyLock) {
			return StdDraw.keysTyped.removeLast();
		}
	}

	/**
	 * Is the keycode currently being pressed? This method takes as an argument
	 * the keycode (corresponding to a physical key). It can handle action keys
	 * (such as F1 and arrow keys) and modifier keys (such as shift and control).
	 * See <a href = "http://download.oracle.com/javase/6/docs/api/java/awt/event/KeyEvent.html">KeyEvent.java</a>
	 * for a description of key codes.
	 * @return true if keycode is currently being pressed, false otherwise
	 */
	public static boolean isKeyPressed(final int keycode) {
		return StdDraw.keysDown.contains(keycode);
	}


	/**
	 * This method cannot be called directly.
	 */
	public void keyTyped(final KeyEvent e) {
		synchronized (StdDraw.keyLock) {
			StdDraw.keysTyped.addFirst(e.getKeyChar());
		}
	}

	/**
	 * This method cannot be called directly.
	 */
	public void keyPressed(final KeyEvent e) {
		StdDraw.keysDown.add(e.getKeyCode());
	}

	/**
	 * This method cannot be called directly.
	 */
	public void keyReleased(final KeyEvent e) {
		StdDraw.keysDown.remove(e.getKeyCode());
	}




	/**
	 * Test client.
	 */
	public static void main(final String[] args) {
		StdDraw.square(.2, .8, .1);
		StdDraw.filledSquare(.8, .8, .2);
		StdDraw.circle(.8, .2, .2);

		StdDraw.setPenColor(StdDraw.BOOK_RED);
		StdDraw.setPenRadius(.02);
		StdDraw.arc(.8, .2, .1, 200, 45);

		// draw a blue diamond
		StdDraw.setPenRadius();
		StdDraw.setPenColor(StdDraw.BOOK_BLUE);
		final double[] x = { .1, .2, .3, .2 };
		final double[] y = { .2, .3, .2, .1 };
		StdDraw.filledPolygon(x, y);

		// text
		StdDraw.setPenColor(StdDraw.BLACK);
		StdDraw.text(0.2, 0.5, "black text");
		StdDraw.setPenColor(StdDraw.WHITE);
		StdDraw.text(0.8, 0.8, "white text");
	}

}
