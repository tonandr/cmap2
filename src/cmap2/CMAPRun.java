package cmap2;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

/**
 * @author Inwoo Chung (gutomitai@gmail.com)
 *
 */
public class CMAPRun {
	
	// Constants.
	public String GENE_LIST_FILE_NAME = "score_n10x10174.csv";
	public String QUERY_UP_SIG_FILE_NAME;
	public String QUERY_DOWN_SIG_FILE_NAME;
	public String RESULT_FILE_NAME;
	
	public static int NUM_GENES = 10174;
	public static int NUM_SIGS = 476251;
	
	// Gene list.
	private int[] geneList = new int[NUM_GENES];
	
	// Query signature data.
	private List<String[]> queryUpSigs = new ArrayList<String[]>();
	private List<String[]> queryDownSigs = new ArrayList<String[]>();
		
	// Run.
	public void run(String[] args) throws FileNotFoundException {
		
		// Load args.
		QUERY_UP_SIG_FILE_NAME = args[0];
		QUERY_DOWN_SIG_FILE_NAME = args[1];
		RESULT_FILE_NAME = args[2];

		
		// Evaluate the fast WTKS calculator.
		CMAP2Updated cmap2 = new CMAP2Updated();
		
		// Load test data.
		// Get a gene list.
		Scanner input = new Scanner(new File(GENE_LIST_FILE_NAME));
		
		// Skip a header.
		input.nextLine();
		
		for (int i = 0; i < NUM_GENES; i++) {
			String val = input.nextLine();
			String[] vals = val.split(",");
			geneList[i] = Integer.parseInt(vals[0]);
		}
		
		input.close();
		
		// Load query signature data.
		// Up.
		input = new Scanner(new File(QUERY_UP_SIG_FILE_NAME));
		
		while (input.hasNextLine()) {
			String val = input.nextLine();
			String[] vals = val.split(",");
			String[] gVals = new String[vals.length - 2];
			
			for (int i = 2; i < vals.length; i++) {
				gVals[i - 2] = vals[i]; 
			}
			
			queryUpSigs.add(gVals);
		}
		
		// Down.
		input = new Scanner(new File(QUERY_DOWN_SIG_FILE_NAME));
		
		while (input.hasNextLine()) {
			String val = input.nextLine();
			String[] vals = val.split(",");
			String[] gVals = new String[vals.length - 2];
			
			for (int i = 2; i < vals.length; i++) {
				gVals[i - 2] = vals[i]; 
			}
			
			queryDownSigs.add(gVals);
		}
		
		input.close();
				
		// Calculate WTKScont.
		int numTestedQuerySigs = queryUpSigs.size();
		double[] WTKScont = new double[numTestedQuerySigs * NUM_SIGS];
		
		cmap2.init(geneList);
				
		for (int i = 0; i < numTestedQuerySigs; i++) {
			double[] r = cmap2.getWTKScomb(queryUpSigs.get(i), queryDownSigs.get(i));
			
			for (int j = 0; j < NUM_SIGS; j++) {
				WTKScont[i * NUM_SIGS + j] = r[j];
			}
		}
				
		// Save result.
		ByteBuffer byteBuffer = ByteBuffer.allocate(WTKScont.length * 8);
		byteBuffer.order(ByteOrder.LITTLE_ENDIAN);
		
		for (double v : WTKScont) {
			byteBuffer.putDouble(v);
		}
		
		byte[] buf = byteBuffer.array();
		
		try {
			Files.write(Paths.get(RESULT_FILE_NAME), buf);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	// Main.
	public static void main(String[] args) {
		CMAPRun cmapRun = new CMAPRun();
		
		try {
			cmapRun.run(args);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
