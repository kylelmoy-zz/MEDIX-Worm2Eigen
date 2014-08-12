package com.kylelmoy.wrm2eig;

import java.io.DataOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Scanner;

import org.la4j.matrix.dense.Basic2DMatrix;

import Jama.Matrix;

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
	private static class EigenPair implements Comparable {
		final public double value;
		final public double[] vector;
		public EigenPair(double value, double[] vector) {
			this.value = value;
			this.vector = vector;
		}
		@Override
		public int compareTo(Object obj) {
			EigenPair that = (EigenPair)obj;
			if (this.value > that.value) return -1;
			if (this.value < that.value) return 1;
			return 0;
		}
		public String toString() {
			return ""+value;
		}
	}
	public static void main(String[] args) throws Exception {
		int t = Runtime.getRuntime().availableProcessors();
		int d = 49;
		int n = d - 1;
		int c = 48;
		//Usage:
		//	wrm2eig input output
		// List of methods
		
		
		//Parse text skeleton points
		System.out.println("Parsing input...");
		long time = System.currentTimeMillis();
		DataFile input = parseInputLoops(new File("data/skeleton.txt"), d);
		//DataFile input = new DataFile(new File("data/input.dat"));
		input.writeToFile(new File("data/input.dat"));
		System.out.println("\tComplete: " + (System.currentTimeMillis() - time) + "ms");
		
		
		//Calculate Vectors
		System.out.println("Calculating vectors...");
		time = System.currentTimeMillis();
		Thread[] threads = new Thread[t];
		//Divide workload for t cores
		DataFile[] jobs = input.split(t);
		//Create t computers
		ComputeJob[] compute = new ComputeJob[t];
		for (int i = 0; i < t; i ++) {
			compute[i] = new ComputeJob(jobs[i]);
		}
		//Construct t threads and begin compute
		System.out.println("\tStarting compute threads...");
		for (int i = 0; i < t; i ++) {
			threads[i] = new Thread(compute[i], "Job " + i);
			threads[i].start();
		}
		//Wait for thread completion
		System.out.println("\tWaiting for thread completion...");
		for (Thread thread : threads) {
			thread.join();

			System.out.println("\t\t" + thread.getName() + " complete...");
		}
		//Join results
		System.out.println("\tJoining results...");
		DataFile vectors = compute[0].getResult();
		for (int i = 1; i < t; i ++) {
			vectors = vectors.join(compute[i].getResult());
		}
		//DataFile vectors = calculateVectors(input);
		vectors.writeToFile(new File("data/vectors.dat"));
		//DataFile vectors = new DataFile(new File("data/vectors.dat"));
		System.out.println("\tComplete: " + (System.currentTimeMillis() - time) + "ms");

		//PCA
		System.out.println("Calculating principal components...");
		time = System.currentTimeMillis();
		DataFile components = calculatePrincipalComponents(vectors);
		//DataFile components = new DataFile(new File("data/components.dat"));
		components.writeToFile(new File("data/components.dat"));
		System.out.println("\tComplete: " + (System.currentTimeMillis() - time) + "ms");
		
		//Calculate amplitudes
		for (int i = 1; i <= c; i ++) {
			System.out.println(i);
			System.out.println("Calculating amplitudes...");
			DataFile amp = calculateAmplitudes(vectors,components,i);
			//amp.writeToFile(new File("data/amplitudes.dat"));
			
			
			//Project
			System.out.println("Calculating projection...");
			DataFile projected = projectData(amp, components, i);
			projected.writeToFile(new File("data/n48/" + i + ".dat"));
		}
		//Done
		System.out.println("Done!");
	}
	private static DataFile parseInputLoops(File file, int d) throws FileNotFoundException {
		ArrayList<Integer> data = new ArrayList<Integer>();
		Scanner input = new Scanner(file);
		Scanner loops = new Scanner(new File("data/isLoop.txt"));
		int caseCount = 0;
		int pointCount = 0;
		System.out.println("\tReading text...");
		while (input.hasNext()) {
			boolean isLoop = loops.nextInt() == 1 ? true : false;
			String caseLine = input.nextLine();
			if (!isLoop) continue;
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
			if (caseCount % 1000 == 0)
				System.out.println("\t\t" + caseCount);
		}
		input.close();
		System.out.println("\tNumber of cases: " + caseCount);
		System.out.println("\tNumber of points: " + pointCount);

		DataFile output = new DataFile((caseCount * d * 2), caseCount);
		//Downsample
		System.out.println("\tDown sampling skeleton to " + d + " points...");
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
		Scanner input = new Scanner(file);
		int caseCount = 0;
		int pointCount = 0;
		System.out.println("\tReading text...");
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
			if (caseCount % 1000 == 0)
				System.out.println("\t\t" + caseCount);
		}
		input.close();
		System.out.println("\tNumber of cases: " + caseCount);
		System.out.println("\tNumber of points: " + pointCount);

		DataFile output = new DataFile((caseCount * d * 2), caseCount);
		//Downsample
		System.out.println("\tDown sampling skeleton to " + d + " points...");
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
	
	/**
	 * Calculates eigenvectors from the covariance matrix of the vector data,
	 * then produces a <code>Matrix</code> of eigenvectors ordered by eigenvalue.
	 * @param data The vector data
	 * @return a <code>Matrix</code> of eigenvectors (components) ordered by greatest eigenvalue
	 */
	private static DataFile calculatePrincipalComponents(DataFile data) {
		int n = data.caseLength();
		int f = data.caseCount();
		Matrix covariance = covar(data);
		Matrix eigenVector = covariance.eig().getV();
		Matrix eigenValue = covariance.eig().getD();
		ArrayList<EigenPair> pq = new ArrayList<EigenPair>();
		for (int i = 0; i < n; i ++) {
			double value = eigenValue.get(i, i);
			double[] vector = new double[n];
			for (int j = 0; j < n; j++) {
				vector[j] = eigenVector.get(j, i);
			}
			pq.add(new EigenPair(value, vector));
		}
		Collections.sort(pq);
		double[] principalComponents = new double[n*n];
		for (int i = 0; i < n; i ++) {
			for (int j = 0; j < n; j ++) {
				principalComponents[(i*n) + j] = pq.get(i).vector[j];
			}
		}
		//Matrix pc = new Matrix(principalComponents);
		//pc.print(n, n);
		return new DataFile(principalComponents,n);
	}
	
	/**
	 * Transforms the vector data with some number of components.
	 * @param vectors The vector data
	 * @param components The components (in order of greatest eigenvalue)
	 * @param numComponents How many components to transform with
	 * @return The transformed data (amplitudes)
	 */
	private static DataFile calculateAmplitudes(DataFile vectors, DataFile components, int numComponents) {
		//Jama is such a worthless library
		Basic2DMatrix pc = new Basic2DMatrix(components.array());
		org.la4j.matrix.Matrix data = (new Basic2DMatrix(vectors.array())).transpose();
		org.la4j.matrix.Matrix feature = pc.sliceTopLeft(numComponents, pc.columns());
		//System.out.println(feature);
		//System.out.println(feature.rows() + "x" + feature.columns() + " * " + data.rows() + "x" + data.columns());
		org.la4j.matrix.Matrix transdata = feature.multiply(data);
		double[] matrix = new double[transdata.rows() * transdata.columns()];
		int c = 0;
		for (int j = 0; j < transdata.rows(); j++) {
			for (int i = 0; i < transdata.columns(); i++) {
				matrix[c++] = transdata.get(j, i);
			}
		}
		return new DataFile(matrix,transdata.columns());
	}
	
	/**
	 * Projects vector data using amplitudes and components.
	 * @param transformed The transformed vector data (amplitudes)
	 * @param components The components
	 * @param numComponents The number of components to use
	 * @return The transformed data
	 */
	private static DataFile projectData(DataFile transformed, DataFile components, int numComponents) {
		Basic2DMatrix pc = new Basic2DMatrix(components.array());
		org.la4j.matrix.Matrix feature = pc.sliceTopLeft(numComponents, pc.columns());
		org.la4j.matrix.Matrix trans = new Basic2DMatrix(transformed.array());
		org.la4j.matrix.Matrix projected = feature.transpose().multiply(trans);
		projected = projected.transpose();
		double[] matrix = new double[projected.rows() * projected.columns()];
		int c = 0;
		for (int j = 0; j < projected.rows(); j++) {
			for (int i = 0; i < projected.columns(); i++) {
				matrix[c++] = projected.get(j, i);
			}
		}
		return new DataFile(matrix,projected.rows());
	}
	
	//HELPER METHODS
	/**
	 * Calculates a covariance matrix for the given <code>DataFile</code>
	 * @param data the <code>DataFile</code> to calculate on
	 * @return the covariance <code>Matrix</code>
	 */
	private static Matrix covar(DataFile data) {
		int n = data.caseLength();
		double[][] covar = new double[n][n];
		for (int x = 0; x < n; x++) {
			for (int y = 0; y < n; y++) {
				if (covar[y][x] != 0) {
					covar[x][y] = covar[y][x];
					continue;
				}
				double xMean = mean(data, x);
				double yMean = mean(data, y);
				double sum = 0;
				for (int i = 0; i < data.caseCount(); i++) {
					sum += (data.get((i * n) + x) - xMean) * (data.get((i * n) + y) - yMean);
				}
				covar[x][y] = sum / (double)(data.caseCount()-1);
			}
		}
		return new Matrix(covar);
	}
	/**
	 * Calculates the mean of a column in a <code>DataFile</code>
	 * @param data the data to calculate on
	 * @param col the column
	 * @return the mean as a double
	 */
	private static double mean(DataFile data, int col) {
		int n = data.caseLength();
		double sum = 0;
		for (int i = 0; i < data.caseCount(); i++) {
			sum += data.get((i * n) + col);
		}
		return sum / data.caseCount();
	}
}
