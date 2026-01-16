package ca.mcgill.ecse420.a1;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

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

	/** Wrapper for parallel matrix multiplication with configurable number of threads
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
	 * Measure execution time (ms) for sequential and parallel multiplication.
	 * Prints average times over 'trials' runs and the speedup.
	 */
	public static void measureExecutionTimes(int matrixSize, int trials, int numThreads) {
		System.out.println(
			"Measuring matrix multiplication with size=" + matrixSize + ", trials=" + trials + ", threads=" + numThreads
		);

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
	 * Populates a matrix of given size with randomly generated integers between
	 * 0-10.
	 * 
	 * @param numRows number of rows
	 * @param numCols number of cols
	 * @return matrix
	 */
	private static double[][] generateRandomMatrix (int numRows, int numCols) {
             double matrix[][] = new double[numRows][numCols];
        for (int row = 0 ; row < numRows ; row++ ) {
            for (int col = 0 ; col < numCols ; col++ ) {
                matrix[row][col] = (double) ((int) (Math.random() * 10.0));
            }
        }
        return matrix;
    }
}
