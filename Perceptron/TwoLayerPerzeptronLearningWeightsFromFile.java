import java.util.*;
import java.io.*;

public class TwoLayerPerzeptronLearningWeightsFromFile
{
	public static void main(String[] args)
	{
		//number Of Training Patterns
		int P = 32;
		
		//output vector size M
		int M = 2;
		
		//input vector size N
		int N = 5;
		
		//Learning rate must be between 0.0 and 1.0
		double eta = 0.1;
		
		//teacher_M_Matrix for pattern P and output M
		double[][] teacher_M_Matrix = new double[P][M];
		
		//inputVector_N_Matrix for pattern P and input N;
		double[][] inputVector_N_Matrix = new double[P][N + 1];
		
		//initialize the BIAS weights
		for (int pattern = 0; pattern < P; pattern ++)
		{
			inputVector_N_Matrix[pattern][0] = 1;
		} // for
		
		File file = new File("PA-A-train.dat");
		File weightsFile = new File("weights.dat");
		File result = new File("result.txt");
		PrintWriter writer = null;
		BufferedReader reader = null;
		BufferedReader weightsReader = null;
		
		String line = null;
		int pattern = 0;
		
		try {
			reader = new BufferedReader(new FileReader(file));
			result.createNewFile();
			writer = new PrintWriter(result);
			
			while((line = reader.readLine()) != null) {
				String[] lineComponents = line.split(",");
				
				for (int i = 1; i < N + 1; i ++)
				{
					inputVector_N_Matrix[pattern][i] = Double.parseDouble(lineComponents[i - 1]);
				} // for
				
				//5 and 6 will the respective teacher_M_values
				teacher_M_Matrix[pattern][0] = Double.parseDouble(lineComponents[5]);
				teacher_M_Matrix[pattern][1] = Double.parseDouble(lineComponents[6]);
				
				pattern ++;
			} // while
			
			/*
			//Sample output
			for (int pat = 0; pat < P; pat ++)
			{
				for (int i = 0; i < N + 1; i ++)
					System.out.printf("%.0f,",inputVector_N_Matrix[pat][i]);
				
				System.out.println();
			} // for
			*/
			
			
			
			//Generate weights by obtaining them from file
			//first create our matrix for weights, M will start from zero where e.g. y0 will map y1 in real world
			double[][] weightMatrix = new double[N + 1][M];
		
			weightsReader = new BufferedReader(new FileReader(weightsFile));
			
			String weightsLine = null;
			int m = 0;
			
			//there will be n + 1 weights for each m
			while((weightsLine = weightsReader.readLine()) != null) {
				String[] lineComponents = weightsLine.split(",");
				
				for (int i = 0; i < N + 1; i ++)
				{
					weightMatrix[i][m] = Double.parseDouble(lineComponents[i]);
				} // for
				
				m ++;
			} // while
		
			//Print weights from file
			System.out.println("\nWeights from File:\n");
		
			for (int i = 0; i < N + 1; i ++)
				for (int j = 0; j < M; j ++)
				{
					System.out.println("w" + i + "" + (j + 1) + " = " + weightMatrix[i][j]);
				} // for
			
			
			//Temporary wieght matrix to hold intermediate weights
			double[][] newWeightMatrix = new double[N + 1][M];
			
			
			// *** Calculate the global error for generated weights
			
			double globalError = 0;
			
			for (int pat = 0; pat < P; pat ++)
			{
				globalError += errorFunctionPerPattern(teacher_M_Matrix[pat], weightMatrix, inputVector_N_Matrix[pat], N + 1, M);
			} // for
			
			System.out.println("\nGlobal Error for initial weights = " + globalError);
			
			globalError = 0;
			
			// *** END of global error calculation
			
			double globalErrorAfterCertainIterations = 0;
			int iterationCount = 20000;
			
			//WE WILL LEARN THROUGH ITERATIONS
			for (int iteration = 0; iteration < iterationCount; iteration ++)
			{
				//new weight can be calculated by w_hm_new = w_hm_old + eta * (teacher_m - net_m) * derivative_f_net_m * out_h
				for (int pat = 0; pat < P; pat ++)
				{
					double[] net_M = net_m_array(weightMatrix, inputVector_N_Matrix[pat], N + 1, M);
				
					for (int j = 0; j < M; j ++)
					{
						for (int i = 0; i < N + 1; i ++)
						{
							//newWeightMatrix[i][j] = weightMatrix[i][j] + eta * (teacher_M_Matrix[pat][j] - net_M[j]) * 1 * inputVector_N_Matrix[pat][i];
							newWeightMatrix[i][j] = weightMatrix[i][j] + eta * (teacher_M_Matrix[pat][j] - logisticFunction(net_M[j])) 
													* derivativeOfLogisticFunction(net_M[j]) * inputVector_N_Matrix[pat][i];
						} // for
					}
				} // for
			
				//now assign newWeightMatrix as weightMatrix
				weightMatrix = newWeightMatrix;
			
				//New wieghts
				writer.write("\nNewly Generated Weights:\n");
		
				for (int i = 0; i < N + 1; i ++)
					for (int j = 0; j < M; j ++)
					{
						writer.write("\nw" + i + "" + (j + 1) + " = " + weightMatrix[i][j]);
					} // for
			
				for (int pat = 0; pat < P; pat ++)
				{
					globalError += errorFunctionPerPattern(teacher_M_Matrix[pat], weightMatrix, inputVector_N_Matrix[pat], N + 1, M);
				} // for
			
				writer.write("\n\nGlobal Error after changing weights = " + globalError + "\n");
				
				globalErrorAfterCertainIterations = globalError;
				
				globalError = 0;
				
			} // for
			
			//New wieghts after number of iterations
			System.out.println("\nGenerated Weights After " + iterationCount + " iterations: \n");
				for (int i = 0; i < N + 1; i ++)
					for (int j = 0; j < M; j ++)
					{
						System.out.println("w" + i + "" + (j + 1) + " = " + weightMatrix[i][j]);
					} // for
			
			System.out.println("\nGlobal Error after " + iterationCount + " iterations = " + globalErrorAfterCertainIterations + "\n");
			
			System.out.println("\nCheck out the result.txt file for further details!");
		}
		catch(Exception e) {
			System.err.println("Error while reading a file \"PA-A-train.dat\"" + " or \"weights.dat\"" + ": " + e.getMessage());
		}
		finally {
			try {
			
				if (reader != null)
					reader.close();
				
				if(weightsReader != null)
					weightsReader.close();
			  } // try
			  catch(Exception e) {
			  e.printStackTrace();
			  }
			
			if (writer != null)
			  writer.close();
		} // finally
		
	} // main
	
	private static double errorFunctionPerPattern(double[] teacher_M, 
								 double[][] weightMatrix, double[] inputVector_X, int inputVectorSize_N, int outputVectorSize_M)
	{
		
		double[] net_M = net_m_array(weightMatrix, inputVector_X, inputVectorSize_N, outputVectorSize_M);
		
		
		double sum = 0;
		
		//this will calculate the error function for one pattern P
		for (int j = 0; j < outputVectorSize_M; j ++)
		{
			//sum = sum + (teacher_M[j] - net_M[j]) * (teacher_M[j] - net_M[j]);
			sum = sum + (teacher_M[j] - logisticFunction(net_M[j])) * (teacher_M[j] - logisticFunction(net_M[j]));
		} // for
		
		return sum / 2;
	} // errorFunction
	
	//function to calculate the net_m, which is the sum of the multiplications w_nm * x_n
	private static double[] net_m_array(double[][] weightMatrix, double[] inputVector_X, int inputVectorSize_N, int outputVectorSize_M)
	{
		double sum = 0;
		
		double[] outputVectorNet_M = new double[outputVectorSize_M];
		
		for (int j = 0; j < outputVectorSize_M; j ++)
		{
			for (int i = 0; i < inputVectorSize_N; i ++)
			{
				sum = sum + weightMatrix[i][j] * inputVector_X[i];
			} // for
			
			outputVectorNet_M[j] = sum;
				
			sum = 0;
		} // for
		
		return outputVectorNet_M;
	} // net_m_array
	
	//method to calculate the logisticFunction
	private static double logisticFunction(double z)
	{
		return (1 / (1 + Math.exp(-z))); // return value of sigmoid function
	} // logisticFunction
	
	//method to calculate the derivative of logistic function
	private static double derivativeOfLogisticFunction(double z)
	{
		return logisticFunction(z) * (1 - logisticFunction(z));
	} // derivateOfLogisticFunction
	
} // class TwoLayerPerzeptronLearning