import java.util.*;
import java.io.*;

//A program to accomplish N-K_M RBF network training.
public class RBFNets
{
	public static void main(String[] args)
	{
		System.out.println("Input Vector N and Output Vector M depends on training data and will hard-coded!");
		System.out.println("Number of training patterns will also be hard coded!");
		
		//number Of Training Patterns, as it is shown in trainRBF.pat file
		int P = 999;
		
		//output vector size M, it is shown as 2 in training.dat file
		int M = 1;
		
		//input vector size N, it is 4 in training.dat file
		int N = 2;
		
		//we will get Number of Neurons K from standard input
		System.out.print("\nEnter the number of neurons in RBF Layer: ");
		Scanner input = new Scanner(System.in);
		int K = input.nextInt();

		//learning rate between RBF Layer and Output Layer
		double eta_K_M = 0.01;
		
		//teacher_M_Matrix for pattern P and output M
		double[][] teacher_M_Matrix = new double[P][M];
		
		//inputVector_N_Matrix for pattern P and input N. There is no bias for input layer
		double[][] inputVector_N_Matrix = new double[P][N];
		
		File file = new File("trainRBF.pat");
		File result = new File("result.txt");
		File learningCurve = new File("learning.curve");
		
		PrintWriter writer = null;
		PrintWriter outWriter = null;
		BufferedReader reader = null;
		
		String line = null;
		int pattern = 0;

		// Setting the seeds for the different random number generators
		int seedWeights  = 2345;
		int seedCenters  = 2345;
		int seedWidths   = 2345;	

		try 
		{	
			reader = new BufferedReader(new FileReader(file));
			result.createNewFile();
			learningCurve.createNewFile();
			writer = new PrintWriter(result);
			outWriter = new PrintWriter(learningCurve);	
			
			while((line = reader.readLine()) != null) {
				String[] lineComponents = line.split(",");
				
				for (int i = 0; i < N; i ++)
				{
					inputVector_N_Matrix[pattern][i] = Double.parseDouble(lineComponents[i]);
				} // for
				
				//Here, 2 will be the respective teacher_M_value
				teacher_M_Matrix[pattern][0] = Double.parseDouble(lineComponents[2]);
				
				pattern ++;
			} // while
			
			//Generate weights, we will use random number generator
			//We are setting seed here to reproduce our results
			Random random = new Random(seedWeights);
			
			//Weight Matrix to hold weights between RBF Layer and Output Layer
			double[][] weightMatrix_K_M = new double[K + 1][M];

			//Now generate weights between RBF Layer and Output Layer
			for (int i = 0; i < K + 1; i++)
			{
	
				for (int j = 0; j < M; j++)
				{	
					//random.nextDouble() - 0.5 will guarantee that random number will be between -0.5 and 0.5
					weightMatrix_K_M[i][j] = (random.nextDouble() - 0.5) ;
				} // for
			} // for	
			
			//Generated weights 
			System.out.println("\nGenerated Weights between RBF Layer and Output Layer:");
			
			for (int j = 0; j < M; j++ )
			{
				for (int i = 0; i < K + 1; i++)
				{
					System.out.println("W_k" + i + "" + "_m" + (j + 1) + " = " + weightMatrix_K_M[i][j]);
				} // for
				
				System.out.println();
			} // for
			
			//matrix to hold the values for each element of c_k vector. we will swap the order n and k during matrix creation
			//double[][] centerVectorC_nk_Matrix = new double[N][K];
			double[][] centerVectorC_nk_Matrix = new double[K][N];
			
			//sometimes random number generator can generate same indices which is not useful
			//we will store each randomly generated index in arraylist and recheck newly generated index each time
			//to change it to the different index
			ArrayList<Integer> patIndices = new ArrayList<Integer>();
			
			if (P < K)
			{
				System.out.println("Number of patterns cannot be smaller than number of neurons in RBF Layer.");
				System.out.println("In the case of P < K, our random uniform selection of centers will not work!");
				System.exit(0);
			} // if
			
			//Print out which method of adjustment is used for centers c_k and widths s_k
			System.out.println("Choosen random uniform selection for centers and widths:");
			
			//setting seed for generating indices
			Random rand = new Random(seedCenters); 

			//we will randomly assign centers from the training data
			for (int k = 0; k < K; k ++)
			{	
				//generates random pattern index
				int randomPat = (int) (rand.nextDouble() * P);
				
				while(patIndices.contains(randomPat))
				{
					randomPat = (int) (rand.nextDouble() * P);
				} // while
				
				patIndices.add(randomPat);
				
				//System.out.println("Random Pattern index: " + randomPat);
				writer.write("Random Pattern index: " + randomPat + "\n");
				
				for (int n = 0; n < N; n ++)
				{
					//now assign the values of elements of x vector to c_k vector
					//centerVectorC_nk_Matrix[n][k] = inputVector_N_Matrix[randomPat][n];
					centerVectorC_nk_Matrix[k][n] = inputVector_N_Matrix[randomPat][n];
				} // for
			} // for
			
			System.out.println("------------------ CENTERS -------------");
			
			//print out chosen centers
			for (int k = 0; k < K; k ++)
			{	
				System.out.print("C" + (k + 1) + " = " + "(");
				
				for (int n = 0; n < N; n ++)
				{
					if (n != N - 1)
						//System.out.print(centerVectorC_nk_Matrix[n][k] + ",");
						System.out.print(centerVectorC_nk_Matrix[k][n] + ",");
					else
						//System.out.print(centerVectorC_nk_Matrix[n][k]);
						System.out.print(centerVectorC_nk_Matrix[k][n]);
				} // for
				
				System.out.println(")");
			} // for
			
			System.out.println("----------------------------------------");
			
			
			//width array s_k to hold the width of gaussian bells
			double[] widthsS_k = new double[K];
			
			//random number with seed set for widths.
			Random r = new Random(seedWidths);
			
			//now it is time to set S_k width for each gaussian bell
			for (int k = 0; k < K; k ++)
			{	
				widthsS_k[k] = 0.3 * (r.nextDouble() * 1.0); // seed will be between 0 and 1 and will be multiplied by 0.3.
			} // for
			
			//print out widths
			System.out.println("\n--------------- WIDTHS ---------------");
			
			for (int k = 0; k < K; k ++)
			{	
				System.out.println("S" + (k + 1) + " = " + widthsS_k[k]);
			} // for
			
			System.out.println("----------------------------------------");
			
			//Temporary weight matrix to hold weights between RBF Layer and Output Layer
			double[][] newWeightMatrix_K_M = new double[K + 1][M];
	
			//For each pattern we will have different r_k values
			double[][] outVectorR_k = new double[P][K + 1];
			
			//initialize the BIAS
			for (int pat = 0; pat < P; pat ++)
			{
				outVectorR_k[pat][0] = 1;
			} // for
			
			
			
			// *** Calculate the global error for generated weights and adjusted centers c_k and widths s_k.
			// *** To do this we have to forward through the net once with initial weights, centers and widths
			
			double globalError = 0;
			
			for (int pat = 0; pat < P; pat ++)
			{
				//we will start k from 1 to assgin R_k values starting from R1, because R0 is 1.
				for (int k = 1; k < K + 1; k ++)
				{
					//determines the value of outVectorR_k for this pattern "pat"
					outVectorR_k[pat][k] = gaussianBellFunction(N, inputVector_N_Matrix[pat], centerVectorC_nk_Matrix[k - 1], widthsS_k[k - 1]);
				} // for
				
				globalError += errorFunctionPerPattern(teacher_M_Matrix[pat], weightMatrix_K_M, outVectorR_k[pat], K + 1, M);
			} // for
			
			globalError /= P; //globalError must be calculated with respect to the number of training examples
			
			//write initial global error to the file.
			outWriter.write(globalError + "\n");
			
			System.out.println("\nGlobal Error for initial weights, centers and widths = " + globalError);
			
			//the variable to hold temporary global errors till we reach the global minima
			double globalErrorAfterCertainIterations = globalError;
			
			globalError = 0;
			
			// *** END of global error calculation
			
			
			
			
			//we will do learning through iterations.
			int iterationCount = 1000;
			
			//random number generator for "every epoch a new random sequence"
			Random randSeq = new Random();
			
			for (int iteration = 0; iteration < iterationCount; iteration ++)
			{
				//Algorithm to achieve every epoch a new random sequence initialization for training patterns
				//Create a sequence of indices
				//Number of Indexes will be equal to number of patterns
				int[] sequence = new int[P];
				
				//initialize each array element with index
				for (int i = 0; i < sequence.length; i++)
				{
					sequence[i] = i;
				} // for
				
				//now shuffle the indexes in sequence array
				for (int i = 0; i < sequence.length; i++)
				{
					//Generates a random sequence between and i and sequence.length
					//But the value of sequence.length is exclusive
					int rIndex = randSeq.nextInt(sequence.length - i) + i;
					
					//swap
					int temp = sequence[rIndex];
					sequence[rIndex] = sequence[i];
					sequence[i] = temp;
					
				} // for
				
				//we will go through each pattern and apply weight changes as in the single step learning
				//In each iteration we will get different sequence of patterns
				for (int pat = 0; pat < P; pat++)
				{
					// choose a random pattern for single step learning
					int randomPatternIndex = sequence[pat];
					
					//initialize the BIAS explicitly
					outVectorR_k[randomPatternIndex][0] = 1;
				
					//we will start k from 1 to assgin R_k values starting from R1, because R0 is 1.
					for (int k = 1; k < K + 1; k ++)
					{
						outVectorR_k[randomPatternIndex][k] = gaussianBellFunction(N, inputVector_N_Matrix[randomPatternIndex], 
															  centerVectorC_nk_Matrix[k - 1], widthsS_k[k - 1]);
					} // for
				
					//Now calculate the net_M based on outVectorR_k values
					double[] netSum_M = net_sum_array(weightMatrix_K_M, outVectorR_k[randomPatternIndex], K + 1, M);
					
					//Calculate deltas
					double[] delta_M = delta_M_array(teacher_M_Matrix[randomPatternIndex], netSum_M, M);
					
					//now update the weights at output layer
					for (int m = 0; m < M; m ++)
					{
						for (int k = 0; k < K + 1; k ++)
						{
							//new weight can be calculated by w_km_new = w_km_old + eta * delta_m * r_k
							newWeightMatrix_K_M[k][m] = weightMatrix_K_M[k][m] + eta_K_M * delta_M[m] * outVectorR_k[randomPatternIndex][k];
						} // for
					} // for
					
				} // for (patterns)
			
				//Now assign new Weight Matrix as weight Matrix
				weightMatrix_K_M = newWeightMatrix_K_M;
				
				
				// ----------------- WRITE NEWLY GENERATED WEIGHTS TO THE FILE --------------------
				
				writer.write("---------------------------------------------------------\n");
			
				//Generated wieghts 
				writer.write("\nNew Weights between RBF Layer and Output Layer:\n\n");
				//System.out.print("\nNew Weights between RBF Layer and Output Layer:\n\n");
			
				for (int m = 0; m < M; m ++)
				{
					for (int k = 0; k < K + 1; k ++)
					{
						writer.write("W_k" + k + "" + "_m" + (m + 1) + " = " + weightMatrix_K_M[k][m] + "\n");
						//System.out.print("W_k" + k + "" + "_m" + (m + 1) + " = " + weightMatrix_K_M[k][m] + "\n");
					} // for
				
					writer.write("\n");
				} // for
				
				// ----------------- WRITE GLOBAL ERROR TO THE FILE --------------------
				
				// Again forward through the net to calculate global error
				for (int pat = 0; pat < P; pat ++)
				{
					
					//initialize the BIAS explicitly
					outVectorR_k[pat][0] = 1;
				
					//we will start k from 1 to assign R_k values starting from R1, because R0 is 1.
					for (int k = 1; k < K + 1; k ++)
					{
						outVectorR_k[pat][k] = gaussianBellFunction(N, inputVector_N_Matrix[pat], centerVectorC_nk_Matrix[k - 1], widthsS_k[k - 1]);
					} // for
				
					globalError += errorFunctionPerPattern(teacher_M_Matrix[pat], weightMatrix_K_M, outVectorR_k[pat], K + 1, M);
				} // for

				globalError /= P; //globalError must be calculated with respect to the number of training examples
			
				writer.write("\n\nGlobal Error after changing weights = " + globalError + "\n");
				
				//print global error
				System.out.print("\nGlobal Error after changing weights = " + globalError + "\n");
				
				outWriter.write(globalError + "\n");
				
				globalErrorAfterCertainIterations = globalError;
				
				globalError = 0;
				
			} // for (iterations)
			
			//New weights to be printed  
			System.out.println("\nNew Weights between RBF Layer and Output Layer After " + iterationCount + " iterations: \n");
			
			for (int j = 0; j < M; j ++)
			{
				for (int i = 0; i < K + 1; i ++)
				{
					System.out.println("W_k" + i + "" + "_m" + (j + 1) + " = " + weightMatrix_K_M[i][j]);
				} // for
				
				System.out.println();
			} // for
			
			System.out.println("\nGlobal Error after " + iterationCount + " iterations = " + globalErrorAfterCertainIterations + "\n");
			
			System.out.println("Check out the result.txt file for further details!");
			
		} // try
		catch(FileNotFoundException e) 
		{ 
			System.out.println("File Not Found");
		} // catch
		catch(ArrayIndexOutOfBoundsException e) { 
			System.out.println("Array index out of bounds: " + e.getMessage());
			e.printStackTrace();
		} // catch
		catch(IOException e) { 
			System.out.println("Error in reading input files: " + e.getMessage());
			e.printStackTrace();
		} // catch
		catch(Exception ex) //general exception catching
		{
			System.err.println(ex.getMessage());
		} // catch
		finally {
			try {
			
				if (reader != null)
					reader.close();
				
			  } // try
			  catch(Exception e) {
				e.printStackTrace();
			  } // catch
			
			if (writer != null)
			  writer.close();
			  
			if(outWriter != null)
				outWriter.close();
				
		} // finally
	} // main
	
	//helper method to calculate value of Gaussian bell function. Which r_k = g(d_k)
	//On the other hand following function will calculate e^( - (|| x - c_k || ^ 2 / (2 * s_k^2)) )
	private static double gaussianBellFunction(int vectorSize_N, double[] vector_X, double[] vector_C_k, double S_k)
	{
		double sumOfSquaredDifferences = 0;
		
		for (int n = 0; n < vectorSize_N; n ++)
		{
			//sumOfSquaredDifferences is actually  || x - c_k || ^ 2
			sumOfSquaredDifferences += (vector_X[n] - vector_C_k[n]) * (vector_X[n] - vector_C_k[n]);
		} // for
		
		//System.out.println("S_k" + S_k);
		
		return Math.exp( (-1) * (sumOfSquaredDifferences / (2.0 * S_k * S_k )) );
	} // gaussianBellFunction
	
	//identity function f(z) = z.
	private static double identityFunction(double z)
	{
		return z;
	} // identityFunction
	
	//Delta functions to calculate delta values.
	//This is delta_M function to calculate delta at output layer
	private static double[] delta_M_array(double[] teacher_M, double[] netSum_M, int outputVectorSize_M)
	{
		double[] delta_M = new double[outputVectorSize_M];
		
		for (int m = 0; m < outputVectorSize_M; m ++)
		{
			//derivative of identity function is 1. That's why, (teacher_M - y_M) * derivativeOfIdentityWithRespectToNetSum_M = (teacher_M - y_M) * 1
			//Here y_m = identityFunction(netSum_M).
			delta_M[m] = (teacher_M[m] - identityFunction(netSum_M[m])) * 1;
		} // for
		
		return delta_M;
	} // delta_M_array
	
	//Error function to calculate Error Per Pattern in 3 layer N - K - M RBF netowrk
	//For calculating error per pattern we will need r_k values. 
	//This will calculated beforehand and passed to the method during a method call.
	private static double errorFunctionPerPattern(double[] teacher_M, 
								 double[][] weights_KM, 
								 double[] outVectorR_K, int outVectorSize_K, int outputVectorSize_M)
	{
		
		//calculates the netSum_M for the given arguments to net_sum_array method
		double[] netSum_M = net_sum_array(weights_KM, outVectorR_K, outVectorSize_K, outputVectorSize_M);
		
		
		double sum = 0;
		
		//this will calculate the error function for one pattern P
		for (int m = 0; m < outputVectorSize_M; m ++)
		{	
			//Here y_m = identityFunction(netSum_M).
			sum = sum + (teacher_M[m] - identityFunction(netSum_M[m])) * (teacher_M[m] - identityFunction(netSum_M[m]));
			
		} // for
		
		return sum / 2.0;
	} // errorFunctionPerPattern
	
	//function to calculate the net_sum values for output layer, which is the sum of the multiplications w_km * r_k
	private static double[] net_sum_array(double[][] weightMatrix, double[] inputVector, int inputVectorSize, int outputVectorSize)
	{
		double sum = 0;
		
		double[] outputVectorNet = new double[outputVectorSize];
		
		for (int m = 0; m < outputVectorSize; m ++)
		{
			for (int k = 0; k < inputVectorSize; k ++)
			{
				sum = sum + weightMatrix[k][m] * inputVector[k];
			} // for
			
			outputVectorNet[m] = sum;
				
			sum = 0;
		} // for
		
		return outputVectorNet;
	} // net_sum_array
	
} // class RBFNets
