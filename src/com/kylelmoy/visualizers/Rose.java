package com.kylelmoy.visualizers;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;

import stdlib.StdDraw;

import com.kylelmoy.wrm2eig.DataFile;

public class Rose {
	final static int resolution = 72;
	static int[] histogram;
	public static void main(String[] args) throws IOException {
		DataFile vectors = new DataFile(new File("data/Loop/vectors.dat"));
		histogram = new int[resolution];
		binValues(vectors);
		StdDraw.show(100);
		draw();
		//StdDraw.show(100);
		StdDraw.save("data/graphs/72bins.png");
	}
	private static void binValues(DataFile data) {
		int len = data.caseCount() * data.caseLength();
		for (int i = 0; i < len; i++) {
			double value = data.get(i);
			for (int j = 0; j < resolution; j++) {
				double cmp = (2 * Math.PI / resolution) * (j + 1);
				if (value < cmp) {
					histogram[j] ++;
					break;
				}
			}
		}
	}
	private static void draw() {
		StdDraw.setCanvasSize(1024,1024);
		StdDraw.setXscale(-128,128);
		StdDraw.setYscale(-128,128);
		StdDraw.setPenRadius(0.001);
		
		StdDraw.line(-128, 0, 128, 0);
		StdDraw.line(0, -128, 0, 128);
		int max = 0;
		for (int val : histogram) {
			max = val > max ? val : max;
		}
		Color[] color = {StdDraw.RED,StdDraw.ORANGE,StdDraw.YELLOW,StdDraw.GREEN,StdDraw.BLUE,StdDraw.MAGENTA};
		for (int i = 0; i < resolution; i++) {
			double value = ((histogram[i] / (double)max) * 2048.0);
			double angle1 = (2 * Math.PI / resolution) * i;
			double angle2 = (2 * Math.PI / resolution) * (i + 1);
			double subangle = angle2 - (Math.PI / resolution);
			double[] x = new double[3];
			double[] y = new double[3];
	        double textx = 90 * Math.cos(angle1);
	        double texty = 90 * Math.sin(angle1);
	        double valx = (100 + (i % 2 == 0 ? 10 : 0)) * Math.cos(subangle);
	        double valy = (100 + (i % 2 == 0 ? 10 : 0)) * Math.sin(subangle);
			x[0] = 0;
			y[0] = 0;
			x[1] = value * Math.cos(angle1);
			y[1] = value * Math.sin(angle1);
			x[2] = value * Math.cos(angle2);
			y[2] = value * Math.sin(angle2);
			StdDraw.setPenColor(StdDraw.LIGHT_GRAY);
	        StdDraw.line(0, 0, textx, texty);
			StdDraw.setPenColor(StdDraw.GRAY);
	        StdDraw.line(0, 0, valx, valy);
			StdDraw.setPenColor(color[i % 6]);
			StdDraw.filledPolygon(x, y);
	        DecimalFormat df = new DecimalFormat("#.##");
			StdDraw.setPenColor(StdDraw.BLACK);
			StdDraw.text(textx, texty, "" + df.format(angle1));
			StdDraw.setPenColor(StdDraw.BOOK_BLUE);
			StdDraw.text(valx, valy, "" + histogram[i]);
			//StdDraw.line(0, 0, x1, y1);
			//StdDraw.line(0, 0, x2, y2);
			//StdDraw.line(x1, y1, x2, y2);
		}
	}
}
