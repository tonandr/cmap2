package cmap2;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

/**
 * @author Inwoo Chung (gutomitai@gmail.com)
 */
public class CMAPLib {
	public static String DATA_DIR_PATH = "/Volumes/Time Machine Backups/"; //?
	public static String SCORES_BY_GENE_FILE_NAME = "scoresByGene";
	public static String RANKS_BY_GENE_FILE_NAME = "ranksByGene";
	public static String SCORES_BY_SIG_FILE_NAME = "scoresBySig";
	public static String RANKS_BY_SIG_FILE_NAME = "ranksBySig";
	
	private static int[] toIntArray(byte[] byteArray){
		  int times = Integer.SIZE / Byte.SIZE;
		  int[] ints = new int[byteArray.length / times];
		  for (int i = 0; i < ints.length; i++){
			byte[] rVal = new byte[4];
			int count = 0;
			
			for (int j = 3; j >= 0; j--) {
				rVal[j] = byteArray[i * times + count++];
			}
			  
		    ints[i] = ByteBuffer.wrap(rVal).getInt();
		  }
		  return ints;
	}
		 
	private static double[] toDoubleArray(byte[] byteArray){
		  int times = Double.SIZE / Byte.SIZE;
		  double[] doubles = new double[byteArray.length / times];
		  for (int i = 0; i < doubles.length; i++){
			byte[] rVal = new byte[8];
			int count = 0;
				
			for (int j = 7; j >= 0; j--) {
				rVal[j] = byteArray[i * times + count++];
			}
				  
			doubles[i] = ByteBuffer.wrap(rVal).getDouble();
		  }
		  return doubles;
	}
		 
		
	private static byte[] loadBinaryFile(String filename) {
		  try {
		    return Files.readAllBytes(sanitizedFilename(filename));
		  } catch (IOException e) {
		    StringWriter sw = new StringWriter();
		    e.printStackTrace(new PrintWriter(sw));
		    return new byte[0];
		  }
	}
		 
	private static Path sanitizedFilename(String filename) {
			Path path = null;
			
			if (filename.compareTo(SCORES_BY_GENE_FILE_NAME) == 0) {
				path = Paths.get(DATA_DIR_PATH + "score_data.bin");
			} else if (filename.compareTo(SCORES_BY_SIG_FILE_NAME) == 0) {
				path = Paths.get(DATA_DIR_PATH + "score_trans.bin");
			} else if (filename.compareTo(RANKS_BY_GENE_FILE_NAME) == 0) {
				path = Paths.get(DATA_DIR_PATH + "int32_rank_data.bin");
			} else if (filename.compareTo(RANKS_BY_SIG_FILE_NAME) == 0) {
				path = Paths.get(DATA_DIR_PATH + "int32_rank_trans.bin");
			} else {
				path = Paths.get(DATA_DIR_PATH + filename);
			}
			
			return path;
	}

	private static int saveBinaryFile(String filename, byte[] data) {
		  try {
		    Files.write(sanitizedFilename(filename), data);
		    return 0;
		  } catch (IOException e) {
		    StringWriter sw = new StringWriter();
		    e.printStackTrace(new PrintWriter(sw));
		    return -1;
		  }
	}
	
	public static int[] loadFromIntFile(String filename, long start, int length){
		if (start == 0 && length == -1) {
			return toIntArray(loadBinaryFile(filename));
		} 
		
		byte[] buf = new byte[length * 4];
		FileInputStream fis;
		long skip;
		
		try {
			fis = new FileInputStream(new File(sanitizedFilename(filename).toString()));
			skip = start * 4;
			fis.skip(skip);
			fis.read(buf);
			fis.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		int[] result = toIntArray(buf);
				
		return result;
	}

	public static double[] loadFromDoubleFile(String filename, long start, int length) {
		if (start == 0 && length == -1) {
			return toDoubleArray(loadBinaryFile(filename));
		} 
		
		byte[] buf = new byte[length * 8];
		long skip;
		
		try {
			FileInputStream fis = new FileInputStream(new File(sanitizedFilename(filename).toString()));
			skip = start * 8;
			
			fis.skip(skip);
			fis.read(buf);
			fis.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		double[] result = toDoubleArray(buf);
		
		return result;
	}
	
	public static int saveIntFile(String filename, int[] data) {		
		return saveBinaryFile(filename, toByteArray(data));
	}

	private static byte[] toByteArray(int[] data) {
		ByteBuffer byteBuffer = ByteBuffer.allocate(data.length * 4);
		byteBuffer.order(ByteOrder.LITTLE_ENDIAN);
		
		for (int v : data) {
			byteBuffer.putInt(v);
		}
		
		return byteBuffer.array();
	}

	public static int saveDoubleFile(String filename, double[] data) {
		return saveBinaryFile(filename, toByteArray(data));
	}

	private static byte[] toByteArray(double[] data) {
		ByteBuffer byteBuffer = ByteBuffer.allocate(data.length * 8);
		byteBuffer.order(ByteOrder.LITTLE_ENDIAN);
		
		for (double v : data) {
			byteBuffer.putDouble(v);
		}
		
		return byteBuffer.array();
	}
}
