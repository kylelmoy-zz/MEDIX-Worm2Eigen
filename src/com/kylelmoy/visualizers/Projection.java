package com.kylelmoy.visualizers;
import java.awt.event.KeyEvent;
import java.io.File;
import java.io.IOException;

import stdlib.StdDraw;

import com.kylelmoy.wrm2eig.DataFile;
public class Projection {
	static DataFile vec;
	static DataFile projected;
	static boolean reload = false;
	static int eig = 48;
	static int n = 48;
	public static void main(String[] args) throws IOException {
		StdDraw.setCanvasSize(256,256);
		StdDraw.setXscale(-128,128);
		StdDraw.setYscale(-128,128);
		StdDraw.show(100);
		StdDraw.setPenRadius(0.005);
		vec = new DataFile(new File("data/vectors.dat"));
		projected = new DataFile(new File("data/n48/48.dat"));
		int f = 0;
		int[] interest = {4,5,6,14,24,36};
		while (true) {
			if (StdDraw.isKeyPressed(KeyEvent.VK_Q)) {
				break;
			}
			if (StdDraw.isKeyPressed(KeyEvent.VK_SPACE)) {
				reload = false;
				projected =  new DataFile(new File("data/n48/" + eig + ".dat"));
				draw(f);
			}
			if (StdDraw.isKeyPressed(KeyEvent.VK_UP)) {
				if (eig < 48) {
					reload = true;
					eig ++;
					draw(f);
				}
			}
			if (StdDraw.isKeyPressed(KeyEvent.VK_DOWN)) {
				if (eig > 0) {
					reload = true;
					eig --;
					draw(f);
				}
			}
			if (StdDraw.isKeyPressed(KeyEvent.VK_RIGHT)) {
				if (StdDraw.isKeyPressed(KeyEvent.VK_SHIFT)) {
					f += 10;
					draw(f);
				}
				else draw(f++);
			}
			if (StdDraw.isKeyPressed(KeyEvent.VK_LEFT)) {
				if (StdDraw.isKeyPressed(KeyEvent.VK_SHIFT)) {
					f -= 10;
					draw(f);
				}
				else draw(f--);
			}
			if (StdDraw.isKeyPressed(KeyEvent.VK_S)) {
				for (int i = 0; i < interest.length; i ++) {
					eig = interest[i];
					projected =  new DataFile(new File("data/n48/" + eig + ".dat"));
					draw(f);
					StdDraw.show(10);
					StdDraw.save("data/images/frame " + f + " - " + eig + ".png");
				}			}
			StdDraw.show(100);
		}
	}
	private static void draw(int frame) {
		//Draw Worm
		StdDraw.clear();
		StdDraw.setPenColor(StdDraw.BLACK);
		StdDraw.setPenRadius(0.003);
		StdDraw.setXscale(-64,64);
		StdDraw.setYscale(-64,64);
		StdDraw.text(0, 64, "Frame: " + frame);
		
		if (reload) StdDraw.setPenColor(StdDraw.RED);
		StdDraw.text(0, 56, "Eigenworms: " + eig);

		StdDraw.setPenColor(StdDraw.MAGENTA);
		StdDraw.setPenRadius(0.005);
		//Recreate from vectors
		double[] x = new double[n + 1];
		double[] y = new double[n + 1];
		double dist = 3.0;
		int offset = frame * n;
		for (int i = 0; i < n; i++) {
			double angle = vec.get(offset + i);
			x[i+1] = x[i] + (dist * Math.cos(angle));
			y[i+1] = y[i] + (dist * Math.sin(angle));
		}
		double _x = 0;
		double _y = 0;
		for (int i = 0; i < n + 1; i++) {
			_x += x[i];
			_y += y[i];
		}
		_x /= n + 1;
		_y /= n + 1;
		for (int i = 0; i < n + 1; i++) {
			x[i] -= _x;
			y[i] -= _y;
			StdDraw.point(x[i], y[i]);
		}
		
		//From projected
		StdDraw.setPenColor(StdDraw.GREEN);
		StdDraw.setPenRadius(0.005);
		x = new double[n + 1];
		y = new double[n + 1];
		for (int i = 0; i < n; i++) {
			x[i+1] = x[i] + (dist * Math.cos(projected.get(offset + i)));
			y[i+1] = y[i] + (dist * Math.sin(projected.get(offset + i)));
		}
		_x = 0;
		_y = 0;
		for (int i = 0; i < n + 1; i++) {
			_x += x[i];
			_y += y[i];
		}
		_x /= n + 1;
		_y /= n + 1;
		for (int i = 0; i < n + 1; i++) {
			x[i] -= _x;
			y[i] -= _y;
			StdDraw.point(x[i], y[i]);
		}
		//StdDraw.line(-128, 0, 128, 0);
	}
}
