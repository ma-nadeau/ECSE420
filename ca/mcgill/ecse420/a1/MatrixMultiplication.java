package ca.mcgill.ecse420.a1;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import java.util.*;
import javax.swing.*;
import org.jfree.chart.*;
import org.jfree.chart.plot.*;
import org.jfree.data.xy.*;

public class MatrixMultiplication {

	private static final int NUMBER_THREADS = 1;
	private static final int MATRIX_SIZE = 2000;

	public static void main(String[] args) {

		// Generate two random matrices, same size
		double[][] a = generateRandomMatrix(MATRIX_SIZE, MATRIX_SIZE);
		double[][] b = generateRandomMatrix(MATRIX_SIZE, MATRIX_SIZE);
		sequentialMultiplyMatrix(a, b);
		parallelMultiplyMatrix(a, b);
	}

	/**
	 * 1.1.
	 * Returns the result of a sequential matrix multiplication
	 * The two matrices are randomly generated
	 * 
	 * @param a is the first matrix
	 * @param b is the second matrix
	 * @return the result of the multiplication
	 */
	public static double[][] sequentialMultiplyMatrix(double[][] a, double[][] b) {

		int rowsA = a.length; // number of rows in A
		int rowsB = b.length; // number of rows in B
		int colsA = a[0].length; // number of columns in A
		int colsB = b[0].length; // number of columns in B

		// check if multiplication is possible, i.e., number of columns in A must equal
		// number of rows in B
		if (colsA != rowsB) {
			throw new IllegalArgumentException(
					"Incompatible matrix sizes: A is " + rowsA + "x" + colsA + ", B is " + b.length + "x" + colsB);
		}

		double[][] c = new double[rowsA][colsB]; // result matrix of size rowsA x colsB

		for (int i = 0; i < rowsA; i++) { // for each row of A
			for (int j = 0; j < colsB; j++) { // for each column of B
				c[i][j] = 0; // initialize the result cell
				for (int k = 0; k < colsA; k++) { // for each element in the row of A and column of B
					c[i][j] += a[i][k] * b[k][j]; // c[i][j] = Σ (a[i][k] · b[k][j]), for k = 0 to n − 1
				}
			}
		}
		return c;
	}

	/**
	 * Wrapper for parallel matrix multiplication with configurable number of
	 * threads
	 * 
	 * @param a is the first matrix
	 * @param b is the second matrix
	 * @return the result of the multiplication
	 */
	public static double[][] parallelMultiplyMatrix(double[][] a, double[][] b, int numberOfThreads) {
		// delegate to the configurable-thread version using the default NUMBER_THREADS
		int rowsA = a.length; // number of rows in A
		int rowsB = b.length; // number of rows in B
		int colsA = a[0].length; // number of columns in A
		int colsB = b[0].length; // number of columns in B

		// check if multiplication is possible, i.e., number of columns in A must equal
		// number of rows in B
		if (colsA != rowsB) {
			throw new IllegalArgumentException(
					"Incompatible matrix sizes: A is " + rowsA + "x" + colsA + ", B is " + b.length + "x" + colsB);
		}

		double[][] c = new double[rowsA][colsB]; // result matrix of size rowsA x colsB

		// TODO: implement parallel matrix multiplication using NUMBER_THREADS
		ExecutorService executor = Executors.newFixedThreadPool(numberOfThreads);
		executor.shutdown();

		return c;
	}

	/**
	 * 1.2.
	 * Returns the result of a concurrent matrix multiplication
	 * The two matrices are randomly generated
	 * 
	 * @param a is the first matrix
	 * @param b is the second matrix
	 * @return the result of the multiplication
	 */
	public static double[][] parallelMultiplyMatrix(double[][] a, double[][] b) {

		return parallelMultiplyMatrix(a, b, NUMBER_THREADS);

	}

	/**
	 * 1.3.
	 * Measure execution time (ms) for sequential and parallel multiplication.
	 * Prints average times over 'trials' runs and the speedup.
	 * 
	 * @param matrixSize size of (square) matrices
	 * @param trials     number of trials to average over
	 */
	public static void measureExecutionTimes(int matrixSize, int trials, int numThreads) {
		System.out.println(
				"Measuring matrix multiplication with size=" + matrixSize + ", trials=" + trials + ", threads="
						+ numThreads);

		// generate matrices once and reuse for fairness
		double[][] a = generateRandomMatrix(matrixSize, matrixSize);
		double[][] b = generateRandomMatrix(matrixSize, matrixSize);

		long sequentialTotal = 0;
		long parallelTotal = 0;

		for (int t = 0; t < trials; t++) {

			// Sequential
			long start = System.nanoTime(); // start time
			sequentialMultiplyMatrix(a, b);
			long end = System.nanoTime(); // stop time
			sequentialTotal += (end - start); // sum up times over trials

			// Parallel
			start = System.nanoTime(); // start time
			parallelMultiplyMatrix(a, b, numThreads);
			end = System.nanoTime(); // stop time
			parallelTotal += (end - start); // sum up times over trials
		}

		double seqAvgMs = sequentialTotal / (trials * 1_000_000.0); // convert ns to ms
		double parAvgMs = parallelTotal / (trials * 1_000_000.0);
		double speedup = seqAvgMs / parAvgMs;

		System.out.printf("Average sequential: %.2f ms\n", seqAvgMs);
		System.out.printf("Average parallel  : %.2f ms\n", parAvgMs);
		System.out.printf("Speedup (seq/par) : %.2fx\n", speedup);
	}


	/**
	 * Plots execution time as a function of number of threads
	 * 
	 * @param results map of number of threads to execution time
	 */
	private static void plotComputeTimeAsAFunctionOfThreads(Map<Integer, Double> results) {
		XYSeries series = new XYSeries("Execution Time");

		for (Map.Entry<Integer, Double> entry : results.entrySet()) {
			series.add(entry.getKey(), entry.getValue());
		}

		XYDataset dataset = new XYSeriesCollection(series);

		JFreeChart chart = ChartFactory.createXYLineChart(
				"Matrix Multiplication Execution Time vs Number of Threads",
				"Number of Threads",
				"Average Execution Time (ms)",
				dataset,
				PlotOrientation.VERTICAL,
				true,
				true,
				false
		);

		ChartPanel panel = new ChartPanel(chart);
		panel.setPreferredSize(new java.awt.Dimension(800, 600));

		JFrame frame = new JFrame("Performance Graph");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setContentPane(panel);
		frame.pack();
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);

	}

	/**
	 * 1.4.
	 * Vary the number of threads used in parallel multiplication for matrix sizes
	 * of `4000x4000`.
	 * Plot the execution time as a function of number of threads.
	 */
	private static void experimentWithDifferentThreadCounts() {
		int matrixSize = 4000; // size of (square) matrices
		int trials = 5; // number of trials to average over

		// generate matrices once and reuse for fairness
		double[][] a = generateRandomMatrix(matrixSize, matrixSize);
		double[][] b = generateRandomMatrix(matrixSize, matrixSize);

		int[] threadCounts = { 1, 2, 4, 8, 16, 32, 64 };

		// Map to store results of thread count to average execution time
		Map<Integer, Double> results = new LinkedHashMap<>();

		try {
			for (int threads : threadCounts) { // for each thread count
				long totalTime = 0; // total time for all trials
				for (int t = 0; t < trials; t++) { // for each trial
					long start = System.nanoTime(); // start time
					parallelMultiplyMatrix(a, b, threads); 
					long end = System.nanoTime(); // stop time
					totalTime += (end - start); // sum up times over trials
				}
				double avgMs = totalTime / (trials * 1_000_000.0); 
				results.put(threads, avgMs);
				System.out.println("Threads: " + threads + " | Avg Time (ms): " + avgMs);

			}
		} catch (Exception e) {
			e.printStackTrace();
		} 
		plotComputeTimeAsAFunctionOfThreads(results);

	}


	private static void plotExecutionTimeAsFunctionOfMatrixSizeForParallelAndSequential(
			Map<Integer, Double> sequentialResults,
			Map<Integer, Double> parallelResults) {

		XYSeries seqSeries = new XYSeries("Sequential Execution Time");
		XYSeries parSeries = new XYSeries("Parallel Execution Time");

		for (Map.Entry<Integer, Double> entry : sequentialResults.entrySet()) {
			seqSeries.add(entry.getKey(), entry.getValue());
		}

		for (Map.Entry<Integer, Double> entry : parallelResults.entrySet()) {
			parSeries.add(entry.getKey(), entry.getValue());
		}

		XYDataset dataset = new XYSeriesCollection();
		((XYSeriesCollection) dataset).addSeries(seqSeries);
		((XYSeriesCollection) dataset).addSeries(parSeries);

		JFreeChart chart = ChartFactory.createXYLineChart(
				"Matrix Multiplication Execution Time vs Matrix Size",
				"Matrix Size (N x N)",
				"Average Execution Time (ms)",
				dataset,
				PlotOrientation.VERTICAL,
				true,
				true,
				false
		);

		ChartPanel panel = new ChartPanel(chart);
		panel.setPreferredSize(new java.awt.Dimension(800, 600));

		JFrame frame = new JFrame("Performance Graph");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setContentPane(panel);
		frame.pack();
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);

	}

	
	/** 1.5.
	 *  Vary the size of the matrices being multiplied as 100x100, 200x200, 500x500, 1000x1000, 2000x2000, 3000x3000 and 4000x4000,
	 * and plot the execution time as a function of matrix size for both parallel and sequential multiplication in one graph.
	 * 
	 */
	private static void experimentWithDifferentMatrixSizes() {
		int[] matrixSizes = { 100, 200, 500, 1000, 2000, 3000, 4000 };
		int trials = 5; // number of trials to average over

		// Maps to store results of matrix size to average execution time
		Map<Integer, Double> sequentialResults = new LinkedHashMap<>();
		Map<Integer, Double> parallelResults = new LinkedHashMap<>();

		try {
			for (int size : matrixSizes) { // for each matrix size
				// generate matrices once and reuse for fairness
				double[][] a = generateRandomMatrix(size, size);
				double[][] b = generateRandomMatrix(size, size);

				// Sequential multiplication
				long sequentialTotal = 0;
				for (int t = 0; t < trials; t++) {
					long start = System.nanoTime();
					sequentialMultiplyMatrix(a, b);
					long end = System.nanoTime();
					sequentialTotal += (end - start);
				}
				double seqAvgMs = sequentialTotal / (trials * 1_000_000.0);
				sequentialResults.put(size, seqAvgMs);

				// Parallel multiplication
				long parallelTotal = 0;
				for (int t = 0; t < trials; t++) {
					long start = System.nanoTime();
					parallelMultiplyMatrix(a, b);
					long end = System.nanoTime();
					parallelTotal += (end - start);
				}
				double parAvgMs = parallelTotal / (trials * 1_000_000.0);
				parallelResults.put(size, parAvgMs);

				System.out.println("Matrix Size: " + size + " | Seq Avg Time (ms): " + seqAvgMs + " | Par Avg Time (ms): " + parAvgMs);
			}
		} catch (Exception e) {
			e.printStackTrace();
		} 

		// Plotting can be implemented similarly to the previous method
		plotExecutionTimeAsFunctionOfMatrixSizeForParallelAndSequential(sequentialResults, parallelResults);
	}

	/**
	 * Populates a matrix of given size with randomly generated integers between
	 * 0-10.
	 * 
	 * @param numRows number of rows
	 * @param numCols number of cols
	 * @return matrix
	 */
	private static double[][] generateRandomMatrix(int numRows, int numCols) {
		double matrix[][] = new double[numRows][numCols];
		for (int row = 0; row < numRows; row++) {
			for (int col = 0; col < numCols; col++) {
				matrix[row][col] = (double) ((int) (Math.random() * 10.0));
			}
		}
		return matrix;
	}
}
