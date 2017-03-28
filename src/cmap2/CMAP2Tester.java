package cmap2;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.Scanner;

/**
 * @author Inwoo Chung (gutomitai@gmail.com)
 */
public class CMAP2Tester {
	
	// Constants.
	public static String DATA_DIR_PATH = "/Users/gutoo/topcoder/cmap2/data/";
	public static String GENE_LIST_FILE_NAME = "score_n10x10174.csv";
	public static String QUERY_UP_SIG_FILE_NAME = "offline_query_up_n250.csv";
	public static String QUERY_DOWN_SIG_FILE_NAME = "offline_query_down_n250.csv";
	public static String GROUND_TRUTH_FILE_NAME = "offline_ground_truth_wtks_n250x476251.csv";
	
	public static int NUM_GENES = 10174;
	public static int NUM_SIGS = 476251;
	
	// Gene list.
	private int[] geneList = new int[NUM_GENES];
	
	// Query signature data.
	private List<String[]> queryUpSigs = new ArrayList<String[]>();
	private List<String[]> queryDownSigs = new ArrayList<String[]>();
	
	// WTKS ref.
	private double[][] WTKSref = new double[250][NUM_SIGS];
	
	// Number of tested query signatures.
	private int numTestedQuerySigs = 1; 
	
	// Test.
	public void test() throws FileNotFoundException {
		
		// Evaluate the fast WTKS calculator.
		CMAP2Updated cmap2 = new CMAP2Updated();
		
		// Load test data.
		// Get a gene list.
		Scanner input = new Scanner(new File(DATA_DIR_PATH + GENE_LIST_FILE_NAME));
		
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
		input = new Scanner(new File(DATA_DIR_PATH + QUERY_UP_SIG_FILE_NAME));
		
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
		input = new Scanner(new File(DATA_DIR_PATH + QUERY_DOWN_SIG_FILE_NAME));
		
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
		
		// Load ground truth.
		input = new Scanner(new File(DATA_DIR_PATH + GROUND_TRUTH_FILE_NAME));
		
		// Skip header.
		input.nextLine();
		
		for (int i = 0; i < NUM_SIGS; i++) {
			String val = input.nextLine();
			String[] vals = val.split(",");
			
			for (int j = 0; j < 250; j++) {
				WTKSref[j][i] = Double.parseDouble(vals[j + 1]);
			}
		}
		
		// Evaluate performance according to the number of tested query signatures.
		// Calculate WTKScont.
		double[][] WTKScont = new double[numTestedQuerySigs][NUM_SIGS];
		
		cmap2.init(geneList);
		
		// Calculate running time.
		Date d = new Date();
		long elapsedTime = d.getTime();
		
		for (int i = 0; i < numTestedQuerySigs; i++) {
			WTKScont[i] = cmap2.getWTKScomb(queryUpSigs.get(i), queryDownSigs.get(i));
		}
		
		elapsedTime = d.getTime() - elapsedTime;
		
		// Calculate WTKSdiff.
		double[][] WTKSdiff = new double[numTestedQuerySigs][NUM_SIGS];
		
		for (int i = 0; i < numTestedQuerySigs; i++) {
			double[] diff = new double[NUM_SIGS];
			
			for (int j = 0; j < NUM_SIGS; j++) {
				diff[i] = Math.abs(WTKScont[i][j] - WTKSref[i][j]);
			}
			
			WTKSdiff[i] = diff;
		}
		
		// Count cases greater than 0.001 for difference.
		int count = 0;
		
		for (int i = 0; i < numTestedQuerySigs; i++) {			
			for (int j = 0; j < NUM_SIGS; j++) {
				if (WTKSdiff[i][j] > 0.001) {
					count += 1;
				}
			}
		}
		
		System.out.println("Number of cases greater than 0.001 for difference: " + String.valueOf(count));
		
		// Calculate score.
		double b = 30.0 * 60.0 * 1000.0 * numTestedQuerySigs / 250.0;
		double score = 1000000.0 * b / (b + 4.0 * elapsedTime);
		
		System.out.println("Elapsed time (minute): " + String.valueOf(elapsedTime / (1000.0 * 60.0)));
		System.out.println("Score: " + String.valueOf(score));
	}
	
	// Main.
	public static void main(String[] args) {
		CMAP2Tester cmap2Tester = new CMAP2Tester();
		try {
			cmap2Tester.test();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
