package cmap2;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveAction;

/**
 * @author Inwoo Chung (gutomitai@gmail.com)
 */
public class CMAP2Updated {
	
	/** DEBUG flag. */
	public static boolean DEBUG = true;
	
	/** Fast WTKS calculation model. */
	public class FastWTKSCalModel {
		
		// Constants.
		public String SORTED_SIGNATURES_FILE_PREFIX = "sorted_sigs_"; // Zero index based.
		public String SORTED_SIGNATURES_GENE_FILE_PREFIX = "sorted_sigs_gene_"; // Zero index based.
		public int NUM_TOTAL_SORTED_SIGNATURES_FILES = 21;
		public int NUM_SORTED_SIGNATURES = 23812;
		public int NUM_REMAINED_SORTED_SIGNATURES = 11;
		
		public String SCORES_BY_GENE_FILE_NAME = "scoresByGene";
		public String RANKS_BY_GENE_FILE_NAME = "ranksByGene";
		public String SCORES_BY_SIG_FILE_NAME = "scoresBySig";
		public String RANKS_BY_SIG_FILE_NAME = "ranksBySig";
		
		public int NUM_GENES = 10174;
		public int NUM_SIGS = 476251;
		
		public int NUM_THREADS = 16;
		
		private int[] geneList;
		
		// Parallel executor to compute WTKS comb values.
		public class ParallelWTKScombCalculator extends RecursiveAction { // Java SE version?
			private boolean isParent;
			private int totalSortedSigFileNum;
			private int[] queryUpSig; // Copy.
			private int[] queryDownSig; // Copy.
			private LinkedList<Integer> sigNumsfifo; // Share.
			private double[] WTKScomb; // Share.
			private double[] sortedSigs; // Share.
			private int[] sortedSigsGene; // Share.
			
			private double[] scores = new double[NUM_GENES];
			private int[] genes = new int[NUM_GENES];
			private Map<Integer, Double> WTKSupMap = new HashMap<Integer, Double>();
			private Map<Integer, Double> WTKSdownMap = new HashMap<Integer, Double>();
			private int[] queryUpIndexes;
			private int[] queryDownIndexes;
			private double WTKSupSum = 0.0;
			private double WTKSdownSum = 0.0;
			private double dValup;
			private double dValdown;
			
			private Random rnd;
			public double sigma = 500.0; //?
			public double t0 = 1.0; //?
			public double c = 1.0; //?
			public int numIter = 100; //?
			
			public ParallelWTKScombCalculator(boolean isParent
					, int totalSortedSigFileNum
					, int[] queryUpSig
					, int[] queryDownSig
					, LinkedList<Integer> sigNumsfifo
					, double[] WTKScomb
					, double[] sortedSigs
					, int[] sortedSigsGene) {
				this.isParent = isParent;
				this.totalSortedSigFileNum = totalSortedSigFileNum;
				this.queryUpSig = new int[queryUpSig.length];
				
				int count = 0;
				
				for (int v : queryUpSig)
					this.queryUpSig[count++] = v;
				
				this.queryDownSig = new int[queryDownSig.length];
				
				count = 0;
				
				for (int v : queryDownSig)
					this.queryDownSig[count++] = v;
				
				this.sigNumsfifo = sigNumsfifo;
				this.WTKScomb = WTKScomb;
				this.sortedSigs = sortedSigs;
				this.sortedSigsGene = sortedSigsGene;
				
				dValup = 1.0 / (sortedSigsGene.length - queryUpSig.length);
				dValdown = 1.0 / (sortedSigsGene.length - queryDownSig.length);
			}
			
			protected void computeDirectly() throws InterruptedException {
				
				// Calculate WTKScomb values until there is no remained signature in signature number fifo.
				do {
					synchronized (sigNumsfifo) {						
						if (sigNumsfifo.isEmpty()) 
							break;
					}
					
					// Pick a signature.
					int sigNum = 0;
					
					synchronized(sigNumsfifo) {
						sigNum = sigNumsfifo.pop();
					}

					// Get scores and genes of a picked signature.							
					// Calculate a start index.
					int startIndex = sigNum - totalSortedSigFileNum * NUM_SORTED_SIGNATURES;
					
					synchronized(sortedSigs) {
						synchronized(sortedSigsGene) {												
							for (int i = 0; i < NUM_GENES; i++) {
								scores[i] = sortedSigs[startIndex + i];
								genes[i] = sortedSigsGene[startIndex + i];
							}
						}
					}

					// Extract query genes' indexes in signature.
					queryUpIndexes = extractQueryGeneIndexes(queryUpSig, genes);
					queryDownIndexes = extractQueryGeneIndexes(queryDownSig, genes);
					
					// Calculate each summation value for WTKSup and WTKSdown.
					for (int idx : queryUpIndexes) {
						WTKSupSum += Math.abs(scores[idx]);
					}
					
					for (int idx : queryDownIndexes) {
						WTKSdownSum += Math.abs(scores[idx]);
					}
					
					// Calculate WTKScomb via MCMC.
					rnd = new Random(sigNum);
					
					// Calculate WTKSup.
					double WTKSup;
					
					{
						int[] geneFreq = new int[NUM_GENES]; 
						int gIdx = rnd.nextInt(NUM_GENES);
						geneFreq[gIdx] += 1;
						
						double t;
						
						for (int i = 0; i < numIter; i++) {
							t = 1.0 / (c * Math.log(i + t0));
							double u = rnd.nextDouble();
							int gIdxP = sampleFromProposalFun(gIdx);
							
							double accept = Math.min(1.0
									, (Math.pow(Math.abs(WTKSupFun(gIdxP)), 1.0 / t) * proposalFun(gIdxP, gIdx)) 
										/ (Math.pow(Math.abs(WTKSupFun(gIdx)), 1.0 / t) * proposalFun(gIdx, gIdxP)));
							
							if (u < accept) {
								gIdx = gIdxP;
							} 
						
							geneFreq[gIdx] += 1;	
						}
						
						// Get a gene index with max frequency.
						int gIdxWithMax = 0;
						
						for (int i = 1; i < geneFreq.length; i++) { //?
							if (geneFreq[gIdxWithMax] < geneFreq[i]) {
								gIdxWithMax = i; 
							}
						}
						
						// Get signed WTKSup.
						WTKSup = WTKSupFun(gIdxWithMax); //?
					}
					
					// Calculate WTKSdown.
					double WTKSdown;
					
					{
						int[] geneFreq = new int[NUM_GENES]; 
						int gIdx = rnd.nextInt(NUM_GENES);
						geneFreq[gIdx] += 1;
						
						double t;
						
						for (int i = 0; i < numIter; i++) {
							t = 1.0 / (c * Math.log(i + t0));
							double u = rnd.nextDouble();
							int gIdxP = sampleFromProposalFun(gIdx);
							
							double accept = Math.min(1.0
									, (Math.pow(Math.abs(WTKSdownFun(gIdxP)), 1.0 / t) * proposalFun(gIdxP, gIdx)) 
										/ (Math.pow(Math.abs(WTKSdownFun(gIdx)), 1.0 / t) * proposalFun(gIdx, gIdxP)));
							
							if (u < accept) {
								gIdx = gIdxP;
							} 
						
							geneFreq[gIdx] += 1;	
						}
						
						// Get a gene index with max frequency.
						int gIdxWithMax = 0;
						
						for (int i = 1; i < geneFreq.length; i++) { //?
							if (geneFreq[gIdxWithMax] < geneFreq[i]) {
								gIdxWithMax = i; 
							}
						}
						
						// Get signed WTKSdown.
						WTKSdown = WTKSdownFun(gIdxWithMax); //?
					}
					
					// Calculate WTKScomb.
					double WTKScombV;
					
					if ((WTKSup >= 0.0 && WTKSdown >= 0.0) ||
							(WTKSup < 0.0 && WTKSdown < 0.0)) {
						WTKScombV = 0.0; //?
					} else {
						WTKScombV = (WTKSup - WTKSdown) / 2.0;
					}
					
					// Update WTKScomb.
					synchronized(WTKScomb) {
						WTKScomb[sigNum] = WTKScombV;
					}
				} while(true);
			}
			
			// Get a sample from a proposal function.
			private int sampleFromProposalFun(int geneIndex) {
				int v = (int)(rnd.nextGaussian() * sigma + geneIndex);
				
				if (v < 0) {
					v = 0;
				}
				
				return Math.min(v, NUM_GENES - 1);
			}
			
			// Proposal function.
			private double proposalFun(int gIdx, int gIdxP) {
				double v = Math.pow(Math.E, -1.0 * Math.pow(gIdx - gIdxP, 2.0) / (2.0 * sigma)) 
						/ (sigma * Math.sqrt(2.0 * Math.PI));
				return v; //?
			}
			
			// Extract query gene' indexes in a signature.
			private int[] extractQueryGeneIndexes(int[] queryGenes, int[] sigGenes) {
				int[] queryGeneIndexes = new int[queryGenes.length];
				List<Integer> remainedQueryGenes = new ArrayList<Integer>();
				
				for (int v : queryGenes) {
					remainedQueryGenes.add(v);
				}
				
				// Search for query genes in a signature.
				int count = 0;
				for (int i = 0; i < sigGenes.length; i++) {
					List<Integer> removedQueryGenes = new ArrayList<Integer>();
					
					for (Integer v : remainedQueryGenes) {
						if (genes[i] == v.intValue()) {
							queryGeneIndexes[count++] = i;
							removedQueryGenes.add(v);
						}
					}
					
					remainedQueryGenes.removeAll(removedQueryGenes);
				}
				
				return queryGeneIndexes;
			}
			
			// WTKSup function.
			private double WTKSupFun(int geneIndex) { // Zero index based?
				
				// Check a calculated WTKS value for a gene index close to or same with an input index.
				Set<Integer> calculatedGeneIndexes = WTKSupMap.keySet();
				int searchedIndex = -1;
				int indexMinDiff = Integer.MIN_VALUE;
								
				for (int idx : calculatedGeneIndexes) {
					int indexDiff = idx - geneIndex;
						
					if (indexDiff == 0) {
						return WTKSupMap.get(idx);
					}
					
					if (indexDiff < 0 && indexMinDiff < indexDiff) { //?
						indexMinDiff = indexDiff;
						searchedIndex = idx;
					}
				}
				
				// Check there is a searched index.
				if (searchedIndex == -1) { // None.
					
					// Calculate WTKSss.
					// Sg summation.
					double sgSum = 0.0;
					int numSgSum = 0;
					
					for (int idx : queryUpIndexes) {
						if (idx <= geneIndex) {
							sgSum += scores[idx];
							numSgSum += 1;
						} else {
							break;
						}
					}
					
					// decrement value summation.
					double rSum = -1.0 * ((geneIndex + 1) - numSgSum) * dValup;
					
					// WTKSss.
					WTKSupMap.put(geneIndex, sgSum / WTKSupSum + rSum);
					return WTKSupMap.get(geneIndex);
				}
				
				// Calculate WTKSss.
				// Sg summation.
				double sgSum = 0.0;
				int numSgSum = 0;
				
				for (int idx : queryUpIndexes) {
					if (searchedIndex < idx && idx <= geneIndex) {
						sgSum += scores[idx];
						numSgSum += 1;
					} else {
						break;
					}
				}
				
				// decrement value summation.
				double rSum = -1.0 * ((geneIndex + 1) - (searchedIndex + 1) - numSgSum) * dValup; //?
				
				// WTKSss.
				return WTKSupMap.get(searchedIndex) + sgSum / WTKSupSum + rSum;
			}
			
			// WTKSdown function.
			private double WTKSdownFun(int geneIndex) { // Zero index based?
				
				// Check a calculated WTKS value for a gene index close to or same with an input index.
				Set<Integer> calculatedGeneIndexes = WTKSdownMap.keySet();
				int searchedIndex = -1;
				int indexMinDiff = Integer.MIN_VALUE;
								
				for (int idx : calculatedGeneIndexes) {
					int indexDiff = idx - geneIndex;
						
					if (indexDiff == 0) {
						return WTKSdownMap.get(idx);
					}
					
					if (indexDiff < 0 && indexMinDiff < indexDiff) { //?
						indexMinDiff = indexDiff;
						searchedIndex = idx;
					}
				}
				
				// Check there is a searched index.
				if (searchedIndex == -1) { // None.
					
					// Calculate WTKSss.
					// Sg summation.
					double sgSum = 0.0;
					int numSgSum = 0;
					
					for (int idx : queryDownIndexes) {
						if (idx <= geneIndex) {
							sgSum += scores[idx];
							numSgSum += 1;
						} else {
							break;
						}
					}
					
					// decrement value summation.
					double rSum = -1.0 * ((geneIndex + 1) - numSgSum) * dValdown;
					
					// WTKSss.
					WTKSdownMap.put(geneIndex, sgSum / WTKSdownSum + rSum);
					return WTKSdownMap.get(geneIndex);
				}
				
				// Calculate WTKSss.
				// Sg summation.
				double sgSum = 0.0;
				int numSgSum = 0;
				
				for (int idx : queryDownIndexes) {
					if (searchedIndex < idx && idx <= geneIndex) {
						sgSum += scores[idx];
						numSgSum += 1;
					} else {
						break;
					}
				}
				
				// decrement value summation.
				double rSum = -1.0 * ((geneIndex + 1) - (searchedIndex + 1) - numSgSum) * dValdown; //?
				
				// WTKSss.
				return WTKSdownMap.get(searchedIndex) + sgSum / WTKSdownSum + rSum;
			}
			
			@Override
			protected void compute() {
				if (!isParent) {
					try {
						computeDirectly();
					} catch (InterruptedException e) {
						printMsg(e.getMessage());
					}
					
					return;
				}
				
				// Compute WTKScomb values parallelly.
				// Create tasks according to the number of cores.
				List<ParallelWTKScombCalculator> tasks = new ArrayList<ParallelWTKScombCalculator>();
								
				for (int i = 0; i < NUM_THREADS; i++) {
					tasks.add(new ParallelWTKScombCalculator(false
							, totalSortedSigFileNum
							, queryUpSig
							, queryDownSig
							, sigNumsfifo
							, WTKScomb
							, sortedSigs
							, sortedSigsGene));
				}
				
				invokeAll(tasks);
			}
		}
		
		/** Constructor. */
		public FastWTKSCalModel(int[] geneList) {
			
			// Initialize.
			// Gene list.
			this.geneList = geneList;
			
			// Sort and save the CMAP signature matrix with descending order 
			// according to signature score rank.
			//printMsg("Sort and save the CMAP signature matrix with descending order.");
			//sortSaveSigs();
		}
		
		// Sort and save the CMAP signature matrix with descending order 
		// according to signature score rank.
		private void sortSaveSigs() {
			
			// About not remained sorted signatures.
			printMsg("About not remained sorted signatures");
			for (int i = 0; i < NUM_TOTAL_SORTED_SIGNATURES_FILES - 1; i++) {
				printMsg("Process the " + String.valueOf(i) + " sorted signature file.");
				double[] sortedSigs = new double[NUM_SORTED_SIGNATURES * NUM_GENES];
				int[] sortedSigsGene = new int[NUM_SORTED_SIGNATURES * NUM_GENES];
				
				for (int j = 0; j < NUM_SORTED_SIGNATURES; j++) {
					
					// Calculate a signature index.
					long sigIndex = i * NUM_SORTED_SIGNATURES + j;
					
					// Load scores and ranks by signature.
					double[] scores = CMAPLib.loadFromDoubleFile(SCORES_BY_SIG_FILE_NAME
							, sigIndex * NUM_GENES, NUM_GENES);
					int[] ranks = CMAPLib.loadFromIntFile(RANKS_BY_SIG_FILE_NAME
							, sigIndex * NUM_GENES, NUM_GENES);
					
					// Sort scores with descending order.
					for (int k = 0; k < ranks.length; k++) {
						sortedSigs[j * NUM_GENES + ranks[k] - 1] = scores[k];
						sortedSigsGene[j * NUM_GENES + ranks[k] - 1] = geneList[k];
					}
				}

				// Save sorted signatures into a file.			
				int result = CMAPLib.saveDoubleFile(SORTED_SIGNATURES_FILE_PREFIX + String.valueOf(i), sortedSigs);
				
				if (result == -1) {
					printMsg("Failed to save a sorted signature file.");
				}
				
				result = CMAPLib.saveIntFile(SORTED_SIGNATURES_GENE_FILE_PREFIX + String.valueOf(i), sortedSigsGene);

				if (result == -1) {
					printMsg("Failed to save a sorted signature gene file.");
				}
			}

			// About remained sorted signatures.
			printMsg("About not remained sorted signatures");
			double[] sortedSigs = new double[NUM_REMAINED_SORTED_SIGNATURES * NUM_GENES];
			int[] sortedSigsGene = new int[NUM_REMAINED_SORTED_SIGNATURES * NUM_GENES];
			
			printMsg("Process the " + String.valueOf(NUM_TOTAL_SORTED_SIGNATURES_FILES - 1) + " sorted signature file.");
			for (int j = 0; j < NUM_REMAINED_SORTED_SIGNATURES; j++) {
				
				// Calculate a signature index.
				long sigIndex = (NUM_TOTAL_SORTED_SIGNATURES_FILES - 1) * NUM_SORTED_SIGNATURES + j;
				
				// Load scores and ranks by signature.
				double[] scores = CMAPLib.loadFromDoubleFile(SCORES_BY_SIG_FILE_NAME
						, sigIndex * NUM_GENES, NUM_GENES);
				int[] ranks = CMAPLib.loadFromIntFile(RANKS_BY_SIG_FILE_NAME
						, sigIndex * NUM_GENES, NUM_GENES); // One based?
				
				// Sort scores with descending order.
				for (int k = 0; k < ranks.length; k++) {
					sortedSigs[j * NUM_GENES + ranks[k] - 1] = scores[k];
					sortedSigsGene[j * NUM_GENES + ranks[k] - 1] = geneList[k];
				}
			}

			// Save sorted signatures into a file.	
			int result = CMAPLib.saveDoubleFile(SORTED_SIGNATURES_FILE_PREFIX 
					+ String.valueOf(NUM_TOTAL_SORTED_SIGNATURES_FILES - 1), sortedSigs);
			
			if (result == -1) {
				printMsg("Failed to save a sorted signature file.");
			}
				
			result = CMAPLib.saveIntFile(SORTED_SIGNATURES_GENE_FILE_PREFIX 
					+ String.valueOf(NUM_TOTAL_SORTED_SIGNATURES_FILES - 1), sortedSigsGene);
			
			if (result == -1) {
				printMsg("Failed to save a sorted signature gene file.");
			}
		}
		
		// Calculate WTKS combination values for all signatures.
		public double[] calWTKScomb(int[] queryUpSig, int[] queryDownSig) {
			double[] WTKScomb = new double[NUM_SIGS];
			
			// About not remained sorted signatures.
			for (int i = 0; i < NUM_TOTAL_SORTED_SIGNATURES_FILES - 1; i++) {
				printMsg("Process the " + String.valueOf(i) + " sorted signature file.");
				
				// Load sorted signatures and genes.
				double[] sortedSigs 
					= CMAPLib.loadFromDoubleFile(SORTED_SIGNATURES_FILE_PREFIX + String.valueOf(i), 0, -1);
				int[] sortedSigsGene 
					= CMAPLib.loadFromIntFile(SORTED_SIGNATURES_GENE_FILE_PREFIX + String.valueOf(i), 0, -1);
				
				// Create a signature number list.
				LinkedList<Integer> sigNumsfifo = new LinkedList<Integer>(); // Java SE version?
				
				for (int j = 0; j < NUM_SORTED_SIGNATURES; j++) {
					sigNumsfifo.add(i * NUM_SORTED_SIGNATURES + j);
				}
				
				// Calculate WTKScomb values parallelly.
				ParallelWTKScombCalculator pTask = 
						new ParallelWTKScombCalculator(true
								, i
								, queryUpSig
								, queryDownSig
								, sigNumsfifo
								, WTKScomb
								, sortedSigs
								, sortedSigsGene);
				
				ForkJoinPool pool = new ForkJoinPool();
				pool.invoke(pTask);
			}
			
			// About remained sorted signatures.
			printMsg("Process the " + String.valueOf(NUM_TOTAL_SORTED_SIGNATURES_FILES - 1) + " sorted signature file.");
			// Load sorted signatures and genes.
			double[] sortedSigs 
				= CMAPLib.loadFromDoubleFile(SORTED_SIGNATURES_FILE_PREFIX 
						+ String.valueOf(NUM_TOTAL_SORTED_SIGNATURES_FILES - 1), 0, -1);
			int[] sortedSigsGene 
				= CMAPLib.loadFromIntFile(SORTED_SIGNATURES_GENE_FILE_PREFIX 
						+ String.valueOf(NUM_TOTAL_SORTED_SIGNATURES_FILES - 1), 0, -1);
			
			// Create a signature number list.
			LinkedList<Integer> sigNumsfifo = new LinkedList<Integer>(); //?
			
			for (int j = 0; j < NUM_REMAINED_SORTED_SIGNATURES; j++) {
				sigNumsfifo.add((NUM_TOTAL_SORTED_SIGNATURES_FILES - 1) * NUM_SORTED_SIGNATURES + j);
			}
			
			// Calculate WTKScomb values parallelly.
			ParallelWTKScombCalculator pTask = 
					new ParallelWTKScombCalculator(true
							, NUM_TOTAL_SORTED_SIGNATURES_FILES - 1
							, queryUpSig
							, queryDownSig
							, sigNumsfifo
							, WTKScomb
							, sortedSigs
							, sortedSigsGene);
			
			ForkJoinPool pool = new ForkJoinPool();
			pool.invoke(pTask);
			
			return WTKScomb;
		}
	}
	
	public FastWTKSCalModel fastWTKSCalModel;
	
	/** Init. */
	public int init(int[] geneList) {
		
		// Initialize the fast WTKS calculation model.
		printMsg("Initialize the fast WTKS calculation model.");
		printMsg("#Genes: " + String.valueOf(geneList.length));
		fastWTKSCalModel = new FastWTKSCalModel(geneList);
		
		return 0;
	}
	
	/** Get WTKS combination values. */
	public double[] getWTKScomb(String[] up, String[] down) {
		
		// Get gene lists for up and down query signature.
		int[] queryUpSig = convStrGenes2IntGenes(up);
		int[] queryDownSig = convStrGenes2IntGenes(down);
		
		// Calculate WTKS combination values for all signatures.
		printMsg("Get WTKS combination values");
		return fastWTKSCalModel.calWTKScomb(queryUpSig, queryDownSig);
	}
	
	// Convert a String genes list into an integer genes list.
	private int[] convStrGenes2IntGenes(String[] genes) {
		int[] result = new int[genes.length];
		int count = 0;
		
		for (String v : genes) {
			result[count++] = Integer.parseInt(v);
		}
		
		return result;
	}
	
	// Print a message.
	public static void printMsg(String msg) {
		if (DEBUG)
			System.out.println(msg);
	}
}
