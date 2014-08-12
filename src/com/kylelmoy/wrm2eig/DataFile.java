package com.kylelmoy.wrm2eig;

import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.DoubleBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.util.Arrays;

import Jama.Matrix;

/**
 * An array-backed data transfer object for conveying worm data.
 * @author Kyle Moy
 *
 */
public class DataFile {
	/**
	 * Where the data is stored.
	 */
	private transient double[] data;
	
	/**
	 * The total number data points
	 */
	private final int length;
	
	/**
	 * The number of data points in each case
	 */
	private final int caseLength;
	
	/**
	 * Keeps track of the seek position for writing
	 */
	private int writePointer;
	
	/**
	 * Keeps track of the seek position for reading
	 */
	private int readPointer;
	
	/**
	 * The number of cases
	 */
	private final int caseCount;
	
	/**
	 * Construct a new <code>DataFile</code>, and initialize it with data from a <code>File</code>.
	 * @param file the file whose data will populate the DataFile
	 * @param caseLength the length of each individual case 
	 * @param caseCount the number of cases
	 * @throws IOException if the file cannot be read
	 */
	public DataFile(File file) throws IOException {
		RandomAccessFile randomAccessFile = new RandomAccessFile(file, "r");
		FileChannel fileChannel = randomAccessFile.getChannel();
		MappedByteBuffer mappedByteBuffer = fileChannel.map(FileChannel.MapMode.READ_ONLY, 0, fileChannel.size());
		DoubleBuffer doublebuffer = mappedByteBuffer.asDoubleBuffer();
		length = (int)doublebuffer.get();
		caseLength = (int)doublebuffer.get();
		caseCount = (int)doublebuffer.get();
		data = new double[doublebuffer.remaining()];
		doublebuffer.get(data);
		randomAccessFile.close();
		if (length != data.length)
			throw new Error("Declared length, data length mismatch: " + length + " != " + data.length);
		//if (length % caseLength != 0) 
		//	throw new Error("Declared length, case length mismatch: " + length);
		
		writePointer = length;
		readPointer = 0;
	}
	
	
	/**
	 * Construct a new, empty <code>DataFile</code>.
	 * @param file the file whose data will populate the DataFile
	 * @param caseCount the number of cases
	 * @throws IllegalArgumentException if the length or case count is negative
	 */
	public DataFile(int length, int caseCount) {
		if (length < 0 || caseCount <= 0)
			throw new IllegalArgumentException();
		
		this.length = length;
		caseLength = length / caseCount;
		this.caseCount = caseCount;
		data = new double[length];
		writePointer = 0;
		readPointer = 0;
	}
	
	/**
	 * Construct a new <code>DataFile</code>, and initialize it with data from a <code>double[]</code>.
	 * @param data the data
	 * @param caseLength the length of each individual case 
	 * @throws IllegalArgumentException if the length or case length is negative, or the case length is not divisible by the case length
	 */
	public DataFile(double[] data, int caseLength) {
		if (data.length <= 0 || caseLength <= 0)
			throw new IllegalArgumentException();
		this.data = data;
		length = data.length;
		this.caseLength = caseLength;
		this.caseCount = length / caseLength;
		writePointer = 0;
		readPointer = 0;
	}
	
	/**
	 * Writes the supplied data at the writing position, and increments the writing position.
	 * @param d the data to be written
	 */
	public void add(double d) {
		synchronized (data) {
			data[writePointer++] = d;
		}
	}

	/**
	 * Writes the supplied data at the supplied position.
	 * @param index the position to write to
	 * @param d the data to be written
	 */
	public void add(int index, double d) {
		synchronized (data) {
			data[index] = d;
		}
	}
	
	/**
	 * Retrieves the data at the user-supplied index.
	 * @param index the index of the data to retrieve
	 * @return the data
	 */
	public double get(int index) {
		synchronized (data) {
			return data[index];
		}
	}
	
	/**
	 * Retrieves the next available data point.
	 * @return the next available data point
	 */
	public double next() {
		synchronized (data) {
			if (readPointer >= length)
				throw new IndexOutOfBoundsException("Read pointer: " + readPointer);
			return data[readPointer++];
		}
	}
	
	/**
	 * Retrieves all the data for a specified case
	 * @param index the case to retrieve
	 * @return an array of doubles containing the data of the case
	 */
	public double[] getCase(int index) {
		synchronized (data) {
			if ((index + 1) * caseLength > length)
				throw new IndexOutOfBoundsException("Index: " + index);
			return Arrays.copyOfRange(data, index * caseLength, (index + 1) * caseLength);
		}
	}
	
	/**
	 * Returns the specified range of cases.
	 * @param from The beginning of the slice
	 * @param to The last case in the slice
	 * @return A DataFile containing the specified range of cases
	 */
	public DataFile getSlice(int from, int to) {
		synchronized (data) {
			if (from < 0 || from > length)
				throw new IndexOutOfBoundsException("From: " + from);
			if (to < 0)
				throw new IndexOutOfBoundsException("To: " + to);
			if (to > caseCount)
				to = caseCount;
			return new DataFile(Arrays.copyOfRange(data, from * caseLength, (to * caseLength)), caseLength);
		}
	}
	
	/**
	 * Splits this <code>DataFile</code> into a number of <code>DataFile</code>s of length smaller or equal to the original length divided by the number to split.
	 * @param numSplit The number of <code>DataFile</code>s to split into
	 * @return A <code>DataFile</code> array containing the split number
	 */
	public DataFile[] split(int numSplit) {
		if (caseCount * caseLength != length)
			throw new Error("Internal counters do not match");
		DataFile[] output = new DataFile[numSplit];
		int range = (int) Math.ceil((double)caseCount / (double)numSplit);
		int from, to = 0;
		for (int i = 0; i < numSplit; i++) {
			from = to;
			to = from + range;
			output[i] = this.getSlice(from, to);
		}
		return output;
	}
	
	/**
	 * Concatenate that DataFile to the end of this DataFile.
	 * @param that the DataFile to be joined to this DataFile
	 * @return the resulting joined DataFile
	 */
	public DataFile join(DataFile that) {
		if (this.caseLength() != that.caseLength)
			throw new Error("Case length mismatch: " + this.caseLength() + " != " + that.caseLength());
		DataFile output = new DataFile(this.length() + that.length(), this.caseCount() + that.caseCount());
		for (double d : data)
			output.add(d);
		for (int i = 0; i < that.length(); i++)
			output.add(that.get(i));
		return output;
	}
	
	/**
	 * Dumps the state of this DataFile to a file on the hard disk to be loaded later
	 * @param file the <code>File</code> to write to
	 * @throws IOException if the <code>File</code> cannot be written to
	 * @see com.kylelmoy.wrm2eig.DataFile#DataFile(File)
	 */
	public void writeToFile(File file) throws IOException {
		FileOutputStream fos = new FileOutputStream(file);
		DataOutputStream dos = new DataOutputStream(fos);
		dos.writeDouble(length);
		dos.writeDouble(caseLength);
		dos.writeDouble(caseCount);
		for (double v : data)
			dos.writeDouble(v);
		dos.flush();
		dos.close();
	}
	
	/**
	 * @return a 2D array representation of the backing data array.
	 */
	public double[][] array() {
		int w = caseLength;
		int h = caseCount;
		double[][] output = new double[h][w];
		for (int y = 0; y < h; y++) {	
			for (int x = 0; x < w; x++) {
				output[y][x] = data[(y * w) + x];
			}
		}
		return output;
	}
	
	/**
	 * @return The length of this <code>Data File</code>
	 */
	public int length() {
		return length;
	}
	
	/**
	 * @return The number of cases
	 */
	public int caseCount() {
		return caseCount;
	}
	
	/**
	 * @return The case length
	 */
	public int caseLength() {
		return caseLength;
	}
	
	public boolean equals(Object obj) {
		if (!(obj instanceof DataFile))
			return false;
		DataFile that = (DataFile)obj;
		if (that.length() != this.length())
			return false;
		for (int i = 0; i < length; i ++)
			if (that.get(i) != this.get(i))
				return false;
		return true;
	}
	public static void main (String[] args) {
		double[] d = {1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14};
		DataFile input = new DataFile(d,3);
		Matrix m = new Matrix(input.array());
		m.print(input.caseLength(), input.caseCount());
	}
}