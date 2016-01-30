import java.util.*;
import java.io.*;

/*
	Program to simulate backpropagation of error for 3 layer N - H - M network.
	In this case 4 - 3 - 2 neural network.
	
	In Reality this program covers the two parts of the Programming Assigment PA-B:
	1) Implementation of Multi Layer Perzeptron
	2) Implementation of Backpropagation of Error
*/
public class MlpBackPropagationOfError
{
	public static void main(String[] args)
	{
		//number Of Training Patterns, as it is shown in training.dat file
		int P = 11;
		
		//output vector size M, it is shown as 2 in training.dat file
		int M = 2;
		
		//input vector size N, it is 4 in training.dat file
		int N = 4;
		
		//**** Number of Hidden Neurons in our hidden layer
		int H = 3;
		
		//Learning rate must be between 0.0 and 1.0
		//Eta between input layer and hidden layer
		double eta_N_H = 0.1;
		
		//Eta between hidden layer and ouput layer
		double eta_H_M = 0.05;
		
		//teacher_M_Matrix for pattern P and output M
		double[][] teacher_M_Matrix = new double[P][M];
		
		//inputVector_N_Matrix for pattern P and input N;
		double[][] inputVector_N_Matrix = new double[P][N + 1];
		
		//initialize the BIAS
		for (int pattern = 0; pattern < P; pattern ++)
		{
			inputVector_N_Matrix[pattern][0] = 1;
		} // for
		
		File file = new File("training_preprocessed.dat");
		File result = new File("result.txt");
		File learningCurve = new File("learning.curve");
		
		//We do not have test.dat file, so we use training.dat file as test file.
		File testDataFile = new File("training_preprocessed_for_test.dat");
		
		PrintWriter writer = null;
		PrintWriter outWriter = null;
		BufferedReader reader = null;
		
		BufferedReader testDataReader = null;
		
		String line = null;
		int pattern = 0;
		
		try {
			reader = new BufferedReader(new FileReader(file));
			
			testDataReader = new BufferedReader(new FileReader(testDataFile));
			
			result.createNewFile();
			learningCurve.createNewFile();
			writer = new PrintWriter(result);
			outWriter = new PrintWriter(learningCurve);
			
			while((line = reader.readLine()) != null) {
				String[] lineComponents = line.split(",");
				
				for (int i = 1; i < N + 1; i ++)
				{
					inputVector_N_Matrix[pattern][i] = Double.parseDouble(lineComponents[i - 1]);
				} // for
				
				//Here, 4 and 5 will the respective teacher_M_values
				teacher_M_Matrix[pattern][0] = Double.parseDouble(lineComponents[4]);
				teacher_M_Matrix[pattern][1] = Double.parseDouble(lineComponents[5]);
				
				pattern ++;
			} // while
			
			/*
			//Sample output
			for (int pat = 0; pat < P; pat ++)
			{
				for (int i = 0; i < N + 1; i ++)
					System.out.printf("%.1f  ",inputVector_N_Matrix[pat][i]);
				
				System.out.println();
			} // for
			
			//Sample output for teacher_M_Matrix
			for (int pat = 0; pat < P; pat ++)
			{
				for (int j = 0; j < M; j ++)
					System.out.printf("%.1f  ", teacher_M_Matrix[pat][j]);
				
				System.out.println();
			} // for
			*/
			
			
			//Generate weights, we will use random number generator
			//We are setting seed here to reproduce our results
			Random random = new Random((N + 1) * M);
			
			//Random random = new Random((N + 1) * H);
		
			//We are going to have two weight matrices
			//First for holding weights between input layer and hidden Layer
			//Second will hold weight matrices between hidden layer and output layer
			
			//Weight Matrix to hold weights between Input Layer and Hidden Layer
			double[][] weightMatrix_N_H = new double[N + 1][H];
			
			//Weight Matrix to hold weights between Hidden Layer and Output Layer
			double[][] weightMatrix_H_M = new double[H + 1][M];
		
			//Now generate weights between Input and Hidden Layers
			for (int i = 0; i < N + 1; i ++)
				for (int j = 0; j < H; j ++)
				{
					//random.nextDouble() - 2.0 will guarantee that random number will be between -2.0 and 2.0
					//We also reduce decimal points to 2
					weightMatrix_N_H[i][j] = Double.parseDouble(String.format("%.2f", random.nextDouble() - 2.0));
				} // for
			
			//random = new Random((H + 1) * M);
			
			//Let's do it again for weights between Hidden Layer and Output Layer
			for (int i = 0; i < H + 1; i ++)
				for (int j = 0; j < M; j ++)
				{
					//random.nextDouble() - 2.0 will guarantee that random number will be between -2.0 and 2.0
					//We also reduce decimal points to 2
					weightMatrix_H_M[i][j] = Double.parseDouble(String.format("%.2f", random.nextDouble() - 2.0));
				} // for
			
			
		
			//Generated wieghts 
			System.out.println("\nGenerated Weights between Input and Hidden Layer:\n");
			
			for (int j = 0; j < H; j ++)
			{
				for (int i = 0; i < N + 1; i ++)
				{
					System.out.println("W_x" + i + "" + "_h" + (j + 1) + " = " + weightMatrix_N_H[i][j]);
				} // for
				
				System.out.println();
			} // for
			
			System.out.println("---------------------------------------------------------");
			
			//Generated wieghts 
			System.out.println("\nGenerated Weights between Hidden Layer and Output Layer:\n");
			
			for (int j = 0; j < M; j ++)
			{
				for (int i = 0; i < H + 1; i ++)
				{
					System.out.println("W_h" + i + "" + "_m" + (j + 1) + " = " + weightMatrix_H_M[i][j]);
				} // for
				
				System.out.println();
			} // for
			
			
			//Temporary wieght matrix to hold intermediate weights between Input Layer and Hidden Layer
			double[][] newWeightMatrix_N_H = new double[N + 1][H];
			
			//And temporary weight matrix to hold weights between Hidden Layer and Output Layer
			double[][] newWeightMatrix_H_M = new double[H + 1][M];
			
			
			//For each pattern we will have different out_H values
			double[][] out_H = new double[P][H + 1];
			
			//initialize the BIAS
			for (int pat = 0; pat < P; pat ++)
			{
				out_H[pat][0] = 1;
			} // for
			
			
			
			// *** Calculate the global error for generated weights
			// *** To do this we have to forward through the net once with initial weights
			
			double globalError = 0;
			
			for (int pat = 0; pat < P; pat ++)
			{
				double[] net_H = net_array(weightMatrix_N_H, inputVector_N_Matrix[pat], N + 1, H);
				
				//calculate out_H values
				for (int i = 1; i < H + 1; i ++)
				{
					//out_H[pat][i] = logisticFunction(net_H[i - 1]);
					out_H[pat][i] = hyperTanFunction(net_H[i - 1]);
				} // for
				
				globalError += errorFunctionPerPattern(teacher_M_Matrix[pat], weightMatrix_H_M, out_H[pat], H + 1, M);
			} // for
			
			System.out.println("\nGlobal Error for initial weights = " + globalError);
			
			globalError = 0;
			
			// *** END of global error calculation
			
			
			double globalErrorAfterCertainIterations = 0;
			int iterationCount = 40000;
			
			//WE WILL LEARN THROUGH ITERATIONS
			for (int iteration = 0; iteration < iterationCount; iteration ++)
			{
				//new weight can be calculated by w_hm_new = w_hm_old + eta * (teacher_m - net_m) * derivative_f_net_m * out_h
				for (int pat = 0; pat < P; pat ++)
				{
					double[] net_H = net_array(weightMatrix_N_H, inputVector_N_Matrix[pat], N + 1, H);
					
					//calculate out_H values
					for (int i = 1; i < H + 1; i ++)
					{
						//out_H[pat][i] = logisticFunction(net_H[i - 1]);
						out_H[pat][i] = hyperTanFunction(net_H[i - 1]);
					} // for
					
					out_H[pat][0] = 1;
					
					
					//Now calculate the net_M based on out_H values
					double[] net_M = net_array(weightMatrix_H_M, out_H[pat], H + 1, M);
					
					//Calculate deltas
					double[] delta_M = delta_M_array(teacher_M_Matrix[pat], net_M, M);
					
					double[] delta_H = delta_H_array(delta_M, weightMatrix_H_M, net_H, M, H + 1);
					
					//first update the weights at output layer
					for (int j = 0; j < M; j ++)
					{
						for (int i = 0; i < H + 1; i ++)
						{
							newWeightMatrix_H_M[i][j] = weightMatrix_H_M[i][j] + eta_H_M * delta_M[j] * out_H[pat][i];
						} // for
					} // for
					
					//next update the weights at hidden layer
					for (int j = 0; j < H; j ++)
					{
						for (int i = 0; i < N + 1; i ++)
						{
							//delta_H[j + 1] because delta_H[0] does not exist
							newWeightMatrix_N_H[i][j] = weightMatrix_N_H[i][j] + eta_N_H * delta_H[j + 1] * inputVector_N_Matrix[pat][i];
						} // for
					} // for
					
				} // for
			
				//now assign new Weight Matrix as weight Matrix
				weightMatrix_N_H = newWeightMatrix_N_H;
				
				//AND
				weightMatrix_H_M = newWeightMatrix_H_M;
				
				
				
				// ----------------- WRITE NEWLY GENERATED WEIGHTS TO THE FILE --------------------
				
				//Generated wieghts 
				writer.write("\nNew Weights between Input and Hidden Layer:\n\n");
			
				for (int j = 0; j < H; j ++)
				{
					for (int i = 0; i < N + 1; i ++)
					{
						writer.write("W_x" + i + "" + "_h" + (j + 1) + " = " + weightMatrix_N_H[i][j] + "\n");
					} // for
				
					writer.write("\n");
				} // for
			
				writer.write("---------------------------------------------------------\n");
			
				//Generated wieghts 
				writer.write("\nNew Weights between Hidden Layer and Output Layer:\n\n");
			
				for (int j = 0; j < M; j ++)
				{
					for (int i = 0; i < H + 1; i ++)
					{
						writer.write("W_h" + i + "" + "_m" + (j + 1) + " = " + weightMatrix_H_M[i][j] + "\n");
					} // for
				
					writer.write("\n");
				} // for
				
				
				// ----------------- WRITE GLOBAL ERROR TO THE FILE --------------------
				
				// Again forward through the net to calculate global error
				for (int pat = 0; pat < P; pat ++)
				{
					double[] net_H = net_array(weightMatrix_N_H, inputVector_N_Matrix[pat], N + 1, H);
					
					out_H[pat][0] = 1;
					
					//calculate out_H values
					for (int i = 1; i < H + 1; i ++)
					{
						//out_H[pat][i] = logisticFunction(net_H[i - 1]);
						out_H[pat][i] = hyperTanFunction(net_H[i - 1]);
					} // for
				
					globalError += errorFunctionPerPattern(teacher_M_Matrix[pat], weightMatrix_H_M, out_H[pat], H + 1, M);
				} // for
			
				writer.write("\n\nGlobal Error after changing weights = " + globalError + "\n");
				
				outWriter.write(globalError + "\n");
				
				globalErrorAfterCertainIterations = globalError;
				
				globalError = 0;
				
			} // for
			
			//New wieghts to be printed 
			System.out.println("\nNew Weights between Input and Hidden Layer After " + iterationCount + " iterations: \n");
			
			for (int j = 0; j < H; j ++)
			{
				for (int i = 0; i < N + 1; i ++)
				{
					System.out.println("W_x" + i + "" + "_h" + (j + 1) + " = " + weightMatrix_N_H[i][j]);
				} // for
				
				System.out.println();
			} // for
			
			System.out.println("---------------------------------------------------------");
			
			//New wieghts to be printed  
			System.out.println("\nNew Weights between Hidden Layer and Output Layer After " + iterationCount + " iterations: \n");
			
			for (int j = 0; j < M; j ++)
			{
				for (int i = 0; i < H + 1; i ++)
				{
					System.out.println("W_h" + i + "" + "_m" + (j + 1) + " = " + weightMatrix_H_M[i][j]);
				} // for
				
				System.out.println();
			} // for
			
			System.out.println("\nGlobal Error after " + iterationCount + " iterations = " + globalErrorAfterCertainIterations + "\n");
			
			System.out.println("Check out the result.txt file for further details!");
			
			System.out.print("\nDo you want to test the performance of your Neural Network (Y/N): ");
			
			
			//-------------------- TESTING PROCEDURES -------------------------------
			
			//get Scanner object to scan the standard input
			Scanner scanner = new Scanner(System.in);
			
			String feedback = scanner.next();
			
			if (feedback.equalsIgnoreCase("Y"))
			{
				String testDataLine = null;
				int testDataPattern = 0;
			
				while((testDataLine = testDataReader.readLine()) != null) {
					String[] lineComponents = testDataLine.split(",");
				
					for (int i = 1; i < N + 1; i ++)
					{
						inputVector_N_Matrix[testDataPattern][i] = Double.parseDouble(lineComponents[i - 1]);
					} // for
				
					//set the BIAS to 1.
					inputVector_N_Matrix[testDataPattern][0] = 1;
					
					//Here, 4 and 5 will the respective teacher_M_values
					teacher_M_Matrix[testDataPattern][0] = Double.parseDouble(lineComponents[4]);
					teacher_M_Matrix[testDataPattern][1] = Double.parseDouble(lineComponents[5]);
				
					testDataPattern ++;
				} // while
				
				//We will output nice matrix to show pattern number, teacher_m and f(net_m) values.
				//These are the headers
				System.out.println("\nPattern        Teacher_m_1         Y_m_1         Teacher_m_2         Y_m_2 ");
				
				//Now forward through the net to calculate f(net_m) that will be compared by y_m in the test data.
				for (int pat = 0; pat < P; pat ++)
				{
					double[] net_H = net_array(weightMatrix_N_H, inputVector_N_Matrix[pat], N + 1, H);
					
					//calculate out_H values
					for (int i = 1; i < H + 1; i ++)
					{
						//out_H[pat][i] = logisticFunction(net_H[i - 1]);
						out_H[pat][i] = hyperTanFunction(net_H[i - 1]);
					} // for
					
					//intitialize the BIAS 1.
					out_H[pat][0] = 1;
					
					//Now calculate the net_M based on out_H values
					double[] net_M = net_array(weightMatrix_H_M, out_H[pat], H + 1, M);
					
					System.out.printf("%7d        %11.2f        %7.3f        %11.2f        %7.3f", (pat + 1), teacher_M_Matrix[pat][0], 
									 hyperTanFunction(net_M[0]), teacher_M_Matrix[pat][1], hyperTanFunction(net_M[1]));
					
					//go to the new line
					System.out.println();
					
				} // for
				
			} // if
			
			
		} // try
		catch(Exception e) {
			System.err.println("Error while reading a file \"training_preprocessed.dat\": " + e.getMessage());
			System.err.println("OR Error while reading a file \"trainining_preprocessed_for_test.dat\": " + e.getMessage());
		}
		finally {
			try {
			
				if (reader != null)
					reader.close();
					
				if (testDataReader != null)
					testDataReader.close();
				
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
	
	//Delta functions to calculate delta values in each layer for each pattern
	//This is delta_M function to calculate delta at output layer
	private static double[] delta_M_array(double[] teacher_M, double[] net_M, int outputVectorSize_M)
	{
		double[] delta_M = new double[outputVectorSize_M];
		
		for (int i = 0; i < outputVectorSize_M; i ++)
		{
			//delta_M[i] = (teacher_M[i] - logisticFunction(net_M[i])) * derivativeOfLogisticFunction(net_M[i]);
			delta_M[i] = (teacher_M[i] - hyperTanFunction(net_M[i])) * derivativeOfHyperTanFunction(net_M[i]);
		} // for
		
		return delta_M;
	} // delta_M_array
	
	
	//This is delta_H function to calculate delta at the hidden layer through back propagation
	//we will not consider weights w_h0_k<something>. That's why the loop for H will start from 1
	//And actually passed value for H is H + 1.
	private static double[] delta_H_array(double[] delta_K, double[][] weightMatrix_H_K, double[] net_H, int K, int H)
	{
		
		double sum = 0;
		
		double[] delta_H = new double[H];
		
		for (int i = 1; i < H; i ++)
		{	
			for (int j = 0; j < K; j ++)
			{
				sum += delta_K[j] * weightMatrix_H_K[i][j];
			} // for
			
			//calculated net_H indices start from 0 that's why I used i - 1.
			//and we do not need delta_H[0], that's why i has started from 1.
			
			//delta_H[i] = sum * derivativeOfLogisticFunction(net_H[i - 1]);
			delta_H[i] = sum * derivativeOfHyperTanFunction(net_H[i - 1]);
			
			sum = 0;
		} // for
		
		return delta_H;
	} // delta_M_array
	
	
	//Error function to calculate Error Per Pattern in 3 layer N - H - M netowrk
	//For calculating error per pattern we will need out_h values. 
	//This will calculated beforehand and passed to the method during a method call.
	private static double errorFunctionPerPattern(double[] teacher_M, 
								 double[][] weightBetweenHiddenAndOutputLayers, 
								 double[] outVector_H, int outVectorSize_H, int finalOutputVectorSize_M)
	{
		
		//calculates the net_M for the given arguments to net_array method
		double[] net_M = net_array(weightBetweenHiddenAndOutputLayers, outVector_H, outVectorSize_H, finalOutputVectorSize_M);
		
		
		double sum = 0;
		
		//this will calculate the error function for one pattern P
		for (int j = 0; j < finalOutputVectorSize_M; j ++)
		{			
			//sum = sum + (teacher_M[j] - logisticFunction(net_M[j])) * (teacher_M[j] - logisticFunction(net_M[j]));
			sum = sum + (teacher_M[j] - hyperTanFunction(net_M[j])) * (teacher_M[j] - hyperTanFunction(net_M[j]));
			
		} // for
		
		return sum / 2;
	} // errorFunction
	
	
	//function to calculate the net values for any layer, which is the sum of the multiplications w * out
	private static double[] net_array(double[][] weightMatrix, double[] inputVector, int inputVectorSize, int outputVectorSize)
	{
		double sum = 0;
		
		double[] outputVectorNet = new double[outputVectorSize];
		
		for (int j = 0; j < outputVectorSize; j ++)
		{
			for (int i = 0; i < inputVectorSize; i ++)
			{
				sum = sum + weightMatrix[i][j] * inputVector[i];
			} // for
			
			outputVectorNet[j] = sum;
				
			sum = 0;
		} // for
		
		return outputVectorNet;
	} // net_array
	
	/*
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
	*/
	
	// ------- IMPORTANT OUTPUT VALUES ARE IN THE INTERVAL OF -1 AND 1, SO HYPERBOLIC TANGENT MUST BE USED FOR THAT. -----------
	
	private static double hyperTanFunction(double z)
    {
		if (z < -20.0) return -1.0; // approximation is correct to 30 decimals
		else if (z > 20.0) return 1.0;
		else return Math.tanh(z);
    } // hyperTanFunction
	
	private static double derivativeOfHyperTanFunction(double z)
    {
		return 1 - hyperTanFunction(z) * hyperTanFunction(z);
    } // derivativeOfHyperTanFunction
	
} // class BackPropagationOfError