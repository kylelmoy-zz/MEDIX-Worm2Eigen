package com.kylelmoy.wrm2eig;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;

public class Wrm2Eig {
	private static class ComputeJob implements Runnable {
		private final DataFile job;
		private DataFile result;
		private ComputeJob (DataFile job) {
			this.job = job;
		}
		public DataFile getResult() {
			return result;
		}
		@Override
		public void run() {
			result = calculateVectors(job);
		}
	}
	public static void main(String[] args) throws Exception {
		int t = 4;
		int d = 49;
		//Usage:
		//	wrm2eig input output
		// List of methods
		//Parse text skeleton points
		System.out.println("Parsing input...");
		//DataFile input = parseInputText(new File("data/skeleton.txt"), d);//args[0]);
		DataFile input = new DataFile(new File("data/input.dat"));
		//input.writeToFile(new File("data/input.dat"));
		//Calculate Vectors
		System.out.println("Calculating vectors...");
		long time = System.currentTimeMillis();
		Thread[] threads = new Thread[t];
		
		//Divide workload for t cores
		DataFile[] jobs = input.split(t);
		
		//Create t computers
		ComputeJob[] compute = new ComputeJob[t];
		for (int i = 0; i < t; i ++) {
			compute[i] = new ComputeJob(jobs[i]);
		}
		
		//Construct t threads and begin compute
		for (int i = 0; i < t; i ++) {
			threads[i] = new Thread(compute[i], "Job " + i);
			threads[i].start();
		}
		
		//Wait for thread completion
		for (Thread thread : threads) {
			thread.join();
		}
		
		//Join results
		DataFile vectors = compute[0].getResult();
		for (int i = 1; i < t; i ++) {
			vectors = vectors.join(compute[i].getResult());
		}
		//DataFile vectors = calculateVectors(input);
		
		//vectors.writeToFile(new File("data/multiVectors.dat"));
		//DataFile vectors = new DataFile(new File("data/newVectors.dat"));
		System.out.println("Complete: " + (System.currentTimeMillis() - time) + "ms");
		/*
		System.out.println(vectors.equals(vectors));
		System.out.println(vectors.equals(result));
		System.out.println(result.length() + " == " + vectors.length());
		System.out.println(result.caseCount() + " == " + vectors.caseCount());
		System.out.println(result.caseLength() + " == " + vectors.caseLength());
		*/
		//PCA
		//Project
		System.out.println("Done!");
	}
	
	/**
	 * Loads skeleton points from a text file in the format:
	 * 		|x1,1;y1,1|x1,2;y1,2|...|x1,n;y1,n
	 * 		...
	 * 		|xc,1;yc,1|xc,2;yc,2|...|xc,n;yc,n
	 * (format from Ron Neihaus's MATLAB video decomposition)
	 * @param file the source file
	 * @param d the number of skeleton points to sample down to
	 * @return a <code>DataFile</code> containing the loaded file data
	 * @throws FileNotFoundException if the file cannot be found
	 */
	private static DataFile parseInputText(File file, int d) throws FileNotFoundException {
		ArrayList<Integer> data = new ArrayList<Integer>();
		Scanner input = new Scanner (file);
		int caseCount = 0;
		int pointCount = 0;
		while (input.hasNext()) {
			String caseLine = input.nextLine();
			int skeletonPoints = 0;
			for (char c : caseLine.toCharArray()) {
				if (c == '|')
					skeletonPoints++;
			}
			//Ignore skeleton entries with too few points
			if (skeletonPoints < 100) {
				continue;
			}
			caseLine = caseLine.replaceAll("\\||;", " ");
			Scanner caseInput = new Scanner(caseLine);
			data.add(skeletonPoints);
			pointCount ++;
			for (int i = 0; i < skeletonPoints; i++) {
				data.add(caseInput.nextInt());
				data.add(caseInput.nextInt());
				pointCount += 2;
			}
			caseInput.close();
			caseCount ++;
		}
		input.close();
		System.out.println("Number of cases: " + caseCount);
		System.out.println("Number of points: " + pointCount);

		DataFile output = new DataFile((caseCount * d * 2), caseCount);
		//Downsample
		int offset = 0;
		while (offset < data.size()){
			int length = (int)data.get(offset++);
			double[] x = new double[length];
			double[] y = new double[length];
			int index = offset;
			for (int i = 0; i < length; i++) {
				//int index = offset + (i * 2); //Optimization, ho!
				x[i] = data.get(index);
				y[i] = data.get(index + 1);
				index += 2;
			}
			offset += length * 2;
			
			//Down Sample to d
			double s = (double)length/(double)(d-1);
			double[] dX = new double[d];
			double[] dY = new double[d];
			double c = 0;
			for (int i = 0; i < (d-1); i++) {
				dX[i] = x[(int)c];
				dY[i] = y[(int)c];
				c += s;
			}
			//Always include tail
			dX[d - 1] = x[length-1];
			dY[d - 1] = y[length-1];
			for (int i = 0; i < d; i++) {
				output.add(dX[i]);
				output.add(dY[i]);
			}
		}
		return output;
	}
	
	/**
	 * Calculates the angles between each consecutive skeleton point, and produces a <code>DataFile</code>.
	 * @param data the DataFile containing the skeleton point data
	 * @return a <code>DataFile</code> containing the vector data
	 */
	private static DataFile calculateVectors(DataFile data) {
		int n = (data.caseLength() / 2) - 1;
		DataFile output = new DataFile(data.caseCount() * n,data.caseCount());
		for (int i = 0; i < data.caseCount(); i++) {
			double[] dX = new double[n+1];
			double[] dY = new double[n+1];
			for (int j = 0; j < n+1; j ++) {
				dX[j] = data.next();
				dY[j] = data.next();
			}
			//Calculate Vectors
			double[] vector = new double[n];
			double sum = 0;
			for (int j = 0; j < n; j++) {
				double yDiff = (dY[j] - dY[j + 1]);
				double xDiff = (dX[j] - dX[j + 1]);
				double angle;
				if (xDiff == 0) {
					if (yDiff > 0) {
						angle = Math.PI/2;
					} else {
						angle = (Math.PI/2) * 3;
					}
				} else angle = Math.atan2(yDiff,xDiff);
				vector[j] = angle;
				sum += angle;
			}
			sum /= n;
			
			//Normalize (subtract mean angle, rotates to 0)
			for (int j = 0; j < n; j++) {
				vector[j] -= sum;
				output.add(vector[j]);
			}
		}
		return output;
	}
}
