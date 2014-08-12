import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.DoubleBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;

import com.kylelmoy.wrm2eig.DataFile;


public class test {

	public static void main(String[] args) throws IOException {
		DataFile vectors = new DataFile(new File("data/vectors.dat"));
		DataFile proj = new DataFile(new File("data/projected.dat"));
		RandomAccessFile randomAccessFile = new RandomAccessFile(new File("data/48.dat"), "r");
		FileChannel fileChannel = randomAccessFile.getChannel();
		MappedByteBuffer mappedByteBuffer = fileChannel.map(FileChannel.MapMode.READ_ONLY, 0, fileChannel.size());
		DoubleBuffer doublebuffer = mappedByteBuffer.asDoubleBuffer();
		for (int i = 0; i < 10; i ++) {
			System.out.println(vectors.next() + " == " + proj.next());
		}
	}

}
