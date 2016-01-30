import java.util.*;
import java.io.*;

public class MSOMTraining
{
	//dimension of the input
	private static final int N = 2;
	//number of training patterns
	private static final int P = 1000;
	
	//initial learning rate
	private static final double eta_init = 0.3;
	//final learning rate
	private static final double eta_final = 0.05;
	
	//number of Max steps to be performed during learning
	private static final int stepTMax = 40000;
	
	//writers to write learning results
	private static PrintWriter learningWriter = null;
	private static PrintWriter resultWriter = null;
	
	//static initialization block to initialize global writers
	static
	{
		try
		{
			File learning = new File("learning.txt");
			File result = new File("PA-E.net");
			
			//create the files
			learning.createNewFile();
			result.createNewFile();
			
			//initialize printwriters
			learningWriter = new PrintWriter(learning);
			resultWriter = new PrintWriter(result);
		}
		catch(IOException ex)
		{
			ex.printStackTrace();
		} // catch
		catch(Exception ex)
		{
			ex.printStackTrace();
		} // catch
		finally
		{
			//do nothing
			//writers will be closed at the end of the main method, there are not throwing exceptions
		} // finally
		
	} // static
	
	public static void main(String[] args)
	{
		Scanner input = new Scanner(System.in);
		
		//input vector
		double[][] inputVector_N_Matrix = new double[P][N];
		
		File file = new File("train_PA-E.dat");
		
		BufferedReader reader = null;
		
		String line = null;
		int pattern = 0;

		// Setting the seeds for the different random number generators
		//int seedWeights  = 2345;
		int seedCenters  = 2345;
		int seedWidths   = 2345;	

		try 
		{	
			reader = new BufferedReader(new FileReader(file));
			
			
			while((line = reader.readLine()) != null) {
				String[] lineComponents = line.split(",");
				
				for (int i = 0; i < N; i ++)
				{
					inputVector_N_Matrix[pattern][i] = Double.parseDouble(lineComponents[i]);
				} // for
				
				pattern ++;
			} // while
		}
		catch(FileNotFoundException ex)
		{
			ex.printStackTrace();
		} // catch
		catch(IOException ex)
		{
			ex.printStackTrace();
		} // catch
		catch(Exception ex)
		{
			ex.printStackTrace();	
		} // catch
		finally
		{
			try {
			
				if (reader != null)
					reader.close();
				
			  } // try
			  catch(Exception e) {
				e.printStackTrace();
			  } // catch
		} // finally
		
		System.out.print("Enter number of Self-Organizing Maps: ");
		
		//number of SOMs
		int M = input.nextInt();
		
		//array to hold dimension for each SOM
		int[] dimensionG = new int[M];
		
		//edge length array to hold edge length for each SOM depending on a dimension
		int[][] edgeLengthsF_MG = new int[M][];
		
		//array to hold number of neurons in each som
		int[] numberOfNeuronsK_M = new int[M];
		
		//now get each dimension from standard input
		for (int m = 0; m < M; m++)
		{
			System.out.print("Enter dimension for SOM " + (m + 1) + ", g = ");
			dimensionG[m] = input.nextInt();
			
			//assign an array with the given dimension
			edgeLengthsF_MG[m] = new int[dimensionG[m]];
			
			//current dimension for this SOM
			int currentDimension = dimensionG[m];
			
			//now get edge lengths from the standard input
			for (int g = 0; g < currentDimension; g ++)
			{
				System.out.print("Enter edge length F_m" + (m + 1) + "g" + (g + 1) + " = ");
				edgeLengthsF_MG[m][g] = input.nextInt();
			} // for
			
		} // for
		
		//calculate number of neurons in each SOM
		for (int m = 0; m < M; m ++)
		{
			//current dimension for this SOM
			int currentDimension = dimensionG[m];
			
			//initialize the number of neurons to 1, before multiplication occurs
			numberOfNeuronsK_M[m] = 1;
			
			for (int g = 0; g < currentDimension; g ++)
			{
				numberOfNeuronsK_M[m] *= edgeLengthsF_MG[m][g];
			} // for
		} // for
		
		/*
		//Sample output
		
		System.out.println("---------------------------------");
		
		//print out dimensions
		for (int m = 0; m < M; m ++)
		{
			System.out.println("Dimension for SOM M" + (m + 1) + ", g = " + dimensionG[m]);	
		} // for
		
		System.out.println("---------------------------------");
		
		//now print out edge length in each SOM
		for (int m = 0; m < M; m++)
		{
			System.out.println("For SOM M" + (m + 1) + ":");
			
			//current dimension for this SOM
			int currentDimension = dimensionG[m];
			
			//now print out edge lengths in the standard output
			for (int g = 0; g < currentDimension; g ++)
			{
				System.out.println("Edge length F_m" + (m + 1) + "g" + (g + 1) + " = " + edgeLengthsF_MG[m][g]);
			} // for
			
			System.out.println("------------------------------------------");
			
		} // for
		
		//now print out number of neurons in each som
		for (int m = 0; m < M; m ++)
		{
			System.out.println("No. Of Neurons in SOM M" + (m + 1) + " = " + numberOfNeuronsK_M[m]);
		} // for
		*/
		
		//3 dimensional array to hold components of center vectors of each SOM M
		double[][][] centerVector_Matrix_MKN = new double[M][][];
		
		//assign second component iteratively. Third component will be N.
		for (int m = 0; m < M; m ++)
		{
			//current number of neurons
			int currentK = numberOfNeuronsK_M[m];
			centerVector_Matrix_MKN[m] = new double[currentK][N];
		} // for
		
		int totalK = 0;
		
		//calculate total number of neurons
		for (int m = 0; m < M; m ++)
		{
			totalK += numberOfNeuronsK_M[m];
		} // for
		
		//print and write total number of neurons
		System.out.println("Total Number of Neurons in " + M + "-SOM = " + totalK);
		learningWriter.println("Total Number of Neurons in " + M + "-SOM = " + totalK);
		resultWriter.println("Total Number of Neurons in " + M + "-SOM = " + totalK);
		
		
		//sometimes random number generator can generate same indices which is not useful
		//we will store each randomly generated index in arraylist and recheck newly generated index each time
		//to change it to the different index
		ArrayList<Integer> patIndices = new ArrayList<Integer>();
			
		if (P < totalK)
		{
			System.out.println("Number of patterns cannot be smaller than number of neurons in M-SOM.");
			System.out.println("In the case of P < K, our random uniform selection of centers will not work!");
			System.exit(0);
		} // if
			
		//Print out which method of adjustment is used for centers c_k and widths s_k
		System.out.println("Choosen random uniform selection for centers and widths:");
			
		//setting seed for generating indices
		Random rand = new Random(seedCenters); 

		//we will randomly assign centers from the training data
		for (int m = 0; m < M; m ++)
		{	
			//centerVector_Matrix_InThisSOM will have dimensions noOfNeurons_K * N
			double[][] centerVector_Matrix_InThisSOM = centerVector_Matrix_MKN[m];
			
			//we will do this operation for each SOM, and each som has different number of neurons.
			for (int k = 0; k < centerVector_Matrix_InThisSOM.length; k ++)
			{
				//generates random pattern index
				int randomPat = (int) (rand.nextDouble() * P);
				
				while(patIndices.contains(randomPat))
				{
					randomPat = (int) (rand.nextDouble() * P);
				} // while
				
				patIndices.add(randomPat);
				
				//System.out.println("Random Pattern index: " + randomPat);
				//writer.write("Random Pattern index: " + randomPat + "\n");
				
				for (int n = 0; n < N; n ++)
				{
					//now assign the values of elements of x vector to c_k vector
					centerVector_Matrix_MKN[m][k][n] = inputVector_N_Matrix[randomPat][n];
				} // for
			} // for
			
		} // for
		
		//boundaries are defined by minimum and maximum. We will use lengths.
		//We will say that if the length is smaller that means the vector is minimum
		double[] minLength = new double[M];
		double[] maxLength = new double[M];
		
		//initialize min and max itertively
		for (int m = 0; m < M; m ++)
		{
			//we will take first neuron of SOM m as a starting point
			minLength[m] = maxLength[m] = vectorLength(N, centerVector_Matrix_MKN[m][0]);
		} // for
		
		//now define max vector and min vector for each SOM m
		double[][] minCenterVector = new double[M][N];
		double[][] maxCenterVector = new double[M][N];
		
		//now go through each neuron of each som to find min and max center vectors		
		//find out boundary of centers of each SOM m
		for (int m = 0; m < M; m ++)
		{
			//centerVector_Matrix_InThisSOM will have dimensions noOfNeurons_K * N
			double[][] centerVector_Matrix_InThisSOM = centerVector_Matrix_MKN[m];
			
			for (int k = 0; k < centerVector_Matrix_InThisSOM.length; k ++)
			{
				//find out minimum
				if (vectorLength(N, centerVector_Matrix_MKN[m][k]) <= minLength[m])
				{
					minLength[m] = vectorLength(N, centerVector_Matrix_MKN[m][k]);
					
					for (int n = 0; n < N; n ++)
					{
						minCenterVector[m][n] = centerVector_Matrix_MKN[m][k][n];
					} // for
					
				} // if
				
				//find out maximum
				if (vectorLength(N, centerVector_Matrix_MKN[m][k]) >= maxLength[m])
				{
					maxLength[m] = vectorLength(N, centerVector_Matrix_MKN[m][k]);
					
					for (int n = 0; n < N; n ++)
					{
						maxCenterVector[m][n] = centerVector_Matrix_MKN[m][k][n];
					} // for
					
				} // if
				
			} // for
		} // for
		
		//also write the boundaries to the result file
		System.out.println("------------------ BOUNDARIES -------------");
		resultWriter.println("------------------ BOUNDARIES -------------");
		
		for (int m = 0; m < M; m ++)
		{	
			
			System.out.print("Boundary for SOM m" + (m + 1) + " = (");
			resultWriter.print("Boundary for SOM m" + (m + 1) + " = (");
			
			//print out boundaries
			for (int n = 0; n < N; n ++)
			{	
				if (n != N - 1)
				{
					System.out.print("n=" + (n + 1) + " => " + minCenterVector[m][n] + " to " + maxCenterVector[m][n] + ", ");
					resultWriter.print("n=" + (n + 1) + " => " + minCenterVector[m][n] + " to " + maxCenterVector[m][n] + ", ");
				} // if
				else
				{
					System.out.print("n=" + (n + 1) + " => " + minCenterVector[m][n] + " to " + maxCenterVector[m][n]);
					resultWriter.print("n=" + (n + 1) + " => " + minCenterVector[m][n] + " to " + maxCenterVector[m][n]);
				} // else
				
			} // for
			
			System.out.println(")");
			resultWriter.println(")");
				
		} // for	
		
		System.out.println("----------------------------------------");
		resultWriter.println("----------------------------------------");
		
		//**** WRITE  INITIAL CENTERS *****
		
		resultWriter.println("------------------ INITIAL CENTERS -------------");
		learningWriter.println("------------------ INITIAL CENTERS -------------");
		
		for (int m = 0; m < M; m ++)
		{	
			//centerVector_Matrix_InThisSOM will have dimensions noOfNeurons_K * N
			double[][] centerVector_Matrix_InThisSOM = centerVector_Matrix_MKN[m];
			
			//print out chosen centers
			for (int k = 0; k < centerVector_Matrix_InThisSOM.length; k ++)
			{	
				resultWriter.print("C_m" + (m + 1) + "k" + (k + 1) + " = " + "(");
				learningWriter.print("C_m" + (m + 1) + "k" + (k + 1) + " = " + "(");
				
				for (int n = 0; n < N; n ++)
				{
					if (n != N - 1)
					{
						resultWriter.print(centerVector_Matrix_MKN[m][k][n] + ",");
						learningWriter.print(centerVector_Matrix_MKN[m][k][n] + ",");
					} // if
					else
					{
						resultWriter.print(centerVector_Matrix_MKN[m][k][n]);
						learningWriter.print(centerVector_Matrix_MKN[m][k][n]);
					} // else
				} // for
				
				resultWriter.println(")");
				learningWriter.println(")");
			} // for
			
			resultWriter.println();
			learningWriter.println();
				
		} // for	
		
		resultWriter.println("----------------------------------------");
		learningWriter.println("----------------------------------------");
		
		//**** END OF INITIAL CENTER WRITING *****
		
		
		
		//print out centers which will also write the centers to the file
		printCenters(M, N, centerVector_Matrix_MKN, false);
		
		
		
		//width array s_mk to hold the width of gaussian bells in each SOM
		double[][] widthsS_MK = new double[M][];
		
		//assign second component iteratively.
		for (int m = 0; m < M; m ++)
		{
			//current number of neurons
			int currentK = numberOfNeuronsK_M[m];
			widthsS_MK[m] = new double[currentK];
		} // for
		
		//random number with seed set for widths.
		Random r = new Random(seedWidths);
			
		//now it is time to set S_mk width for each gaussian bell in each SOM m.
		for (int m = 0; m < M; m ++)
		{
			//current number of neurons
			int currentNoOfNeurons_K = widthsS_MK[m].length;
			
			for (int k = 0; k < currentNoOfNeurons_K; k ++)
			{	
				widthsS_MK[m][k] = 0.3 * (r.nextDouble() * 1.0); // seed will be between 0 and 1 and will be multiplied by 0.3.
			} // for
		} // for
		
		
		
		//print out widths, which can also write to a file
		printWidths(M, widthsS_MK);
		
		
		
		//calculate distances in grid beforehand. Array to hold distances between neuron i and neuron j in each SOM M.
		double[][][] distancesInSOM_M = new double[M][][];
		
		//assign values iteratively through calling helper method
		for (int m = 0; m < M; m ++)
		{
			//pass edge lengths and dimension of this SOM
			distancesInSOM_M[m] = distancesInGrid(edgeLengthsF_MG[m], dimensionG[m]);
			
			//check whether or not we get null. If we get null, one of the dimensions entered is bigger than 5
			if (distancesInSOM_M[m] == null)
			{
				System.out.println("One of the dimensions entered is bigger than 5");
				System.out.println("Please enter the enter dimension 1 <= g <= 5!");
				System.exit(0);
			} // if
		} // for
		
		//we define the response r_k array beforehand to minimize memory usage in iterations. Response r_k will be defined for each SOM M-SOM
		double[][] responseR_MK = new double[M][];
				
		//iteratively assign the second component K
		for (int m = 0; m < M; m ++)
		{
			//current number of neurons in SOM m
			int currentK = numberOfNeuronsK_M[m];
			responseR_MK[m] = new double[currentK];
		} // for
		
		
		//temporary center vector matrix array to hold the intermediate values before update step
		double[][][] newCenterVector_Matrix_MKN = new double[M][][];
		
		//initialize it to the current center vector matrix array
		
		/*
		for (int m = 0; m < M; m ++)
		{
			//current number of neurons
			int currentK = numberOfNeuronsK_M[m];
			newCenterVector_Matrix_MKN[m] = new double[currentK][N];
			
			//copy the elements deeply. We avoid object copying which can cause severe problems
			for (int k = 0; k < currentK; k ++)
			{
				for (int n = 0; n < N; n ++)
				{
					newCenterVector_Matrix_MKN[m][k][n] = centerVector_Matrix_MKN[m][k][n];
				} // for
			} // for
			
		} // for
		*/
		
		//assign second component iteratively. Third component will be N.
		for (int m = 0; m < M; m ++)
		{
			////current number of neurons
			int currentK = numberOfNeuronsK_M[m];
			newCenterVector_Matrix_MKN[m] = new double[currentK][N];
		} // for
		
		
		
		//**** Stopping criteria is number of timeStampss
		
		for (int stepT = 0; stepT < stepTMax; stepT ++)
		{
		
		
			//go through each pattern and take each of them as a stimulus and apply it to M_SOM
			for (int pat = 0; pat < P; pat ++)
			{
				//first of all we have to find the winner neuron i which in turn finds winner SOM m
				int winnerSOMIndex = 0;
				int winnerNeuronIndex = 0;
			
				//go through each SOM to find out winner neuron and winner SOM
				for (int m = 0; m < M; m ++)
				{
					//current number of neurons in SOM m
					int currentK = numberOfNeuronsK_M[m];
					for (int k = 0; k < currentK; k ++)
					{
						responseR_MK[m][k] = distanceBetweenVectors(N, inputVector_N_Matrix[pat], centerVector_Matrix_MKN[m][k]);
					} // for
				} // for
			
				//winner neuron i and winnder som m will be based on minimum distance, so find the minimum
				//firstly initialize minimum to first neuron in first som
				double minimumResponse = responseR_MK[0][0];
			
				//these nested for loops will give us the winner neuron i and winnder SOM m
				//then we use winner SOM index to learn only the centers those resides in itself.
				for (int m = 0; m < M; m ++)
				{
					//current number of neurons in SOM m
					int currentK = numberOfNeuronsK_M[m];
					for (int k = 0; k < currentK; k ++)
					{
						if (responseR_MK[m][k] <= minimumResponse)
						{
							minimumResponse = responseR_MK[m][k];
						
							//then update winnderSOMIndex and winnderNeuronIndex of that SOM
							winnerNeuronIndex = k;
							winnerSOMIndex = m;
						} // if
					} // for
				} // for
			
				//after finding winner SOM and winnder neuron, we can apply the learning rule
				//now assign values to new center vector matrix array, which will be newC_k = oldC_k + changeOfC_k;
				for (int k = 0; k < numberOfNeuronsK_M[winnerSOMIndex]; k ++)
				{
					for (int n = 0; n < N; n ++)
					{
						double[] differenceOfXAndC_k = differenceOfVectors(N, inputVector_N_Matrix[pat], centerVector_Matrix_MKN[winnerSOMIndex][k]);
					
						//we will have scalar multiple here which is learningRate(stepT) * neighborhoodFunction(distance(i,j), width, stepT)
						newCenterVector_Matrix_MKN[winnerSOMIndex][k][n] = centerVector_Matrix_MKN[winnerSOMIndex][k][n] 
																		   + learningRate(stepT) 
																		   * neighborhoodFunction(distancesInSOM_M[winnerSOMIndex][winnerNeuronIndex][k],
																							      widthsS_MK[winnerSOMIndex][k], 
																							      stepT)
																		   * differenceOfXAndC_k[n];
					} // for
				} // for
			
			} // for (patterns)
		
			//update the centers
			centerVector_Matrix_MKN = newCenterVector_Matrix_MKN;
		
			//write out new centers
			printCenters(M, N, centerVector_Matrix_MKN, true);
		
		} // for (iterations)
			
		
		//****** ALSO PRINT out AND WRITE centers after stepTMax iterations
		System.out.println("------------------ FINAL CENTERS after " +  stepTMax + " iterations -------------");
		resultWriter.println("------------------ FINAL CENTERS after " +  stepTMax + " iterations -------------");
		
		for (int m = 0; m < M; m ++)
		{	
			//centerVector_Matrix_InThisSOM will have dimensions noOfNeurons_K * N
			double[][] centerVector_Matrix_InThisSOM = centerVector_Matrix_MKN[m];
			
			//print out chosen centers
			for (int k = 0; k < centerVector_Matrix_InThisSOM.length; k ++)
			{	
				System.out.print("C_m" + (m + 1) + "k" + (k + 1) + " = " + "(");
				resultWriter.print("C_m" + (m + 1) + "k" + (k + 1) + " = " + "(");
				
				for (int n = 0; n < N; n ++)
				{
					if (n != N - 1)
					{
						System.out.print(centerVector_Matrix_MKN[m][k][n] + ",");
						resultWriter.print(centerVector_Matrix_MKN[m][k][n] + ",");
					} // if
					else
					{
						System.out.print(centerVector_Matrix_MKN[m][k][n]);
						resultWriter.print(centerVector_Matrix_MKN[m][k][n]);
					} // else
				} // for
				
				System.out.println(")");
				resultWriter.println(")");
			} // for
			
			System.out.println();
			resultWriter.println();
				
		} // for	
		
		System.out.println("----------------------------------------");
		resultWriter.println("----------------------------------------");
		
		//****** END OF FINAL CENTERS WRIITING AND PRINTING
		
		System.out.println("Please check learning.txt file for further details!");
		
		//close writers
		resultWriter.close();
		learningWriter.close();
		
	} // main
	
	//helper method to calculate distances on grid topology, we will calculate distances depending on different topologies
	//where g <= 5
	private static double[][] distancesInGrid(int[] currentSOMEdgeLengthsG, int g)
	{
		switch(g)
		{
			case 1:
			{
				//if g = 1, there is only one edgge length
				int F_mg1 = currentSOMEdgeLengthsG[g - 1];
				//array to hold positions
				double[] positions = new double[F_mg1];
				//in one dimension, the positions are the indices themselves
				for (int i = 0; i < F_mg1; i ++)
				{
					positions[i] = i;
				} // for
				
				//two dimensional distance array to hold a distance between neuron and i and j which are positioned in a grid
				double[][] distances = new double[F_mg1][F_mg1];
				
				//now calculate the distances between each neuron i and j
				for (int i = 0; i < F_mg1; i ++)
				{
					for (int j = 0; j < F_mg1; j ++)
					{
						distances[i][j] = Math.abs(positions[j] - positions[i]);
						//System.out.println("dist(" + (i + 1) + "," + (j + 1) + ") = " + distances[i][j]);
					} // for
				} // for
				
				return distances;
			} // case 1
			case 2:
			{
				//in two dimensions we have 2 edge lengths
				int F_mg1 = currentSOMEdgeLengthsG[g - 2];
				int F_mg2 = currentSOMEdgeLengthsG[g - 1];
				
				//we will hold the relative position in a string in a from e.g. "a=0, b=0"
				//we will have F_mg1 * F_mg2 total positions
				String[] positions = new String[F_mg1 * F_mg2];
				
				//now assign dimensions for each position and hold them in String array
				int index = 0;
				for (int a = 0; a < F_mg1; a ++)
				{	
					for (int b = 0; b < F_mg2; b ++)
					{
						//we will use comma as a delimiter, the we will split on that.
						positions[index] = a + "," + b;
						index ++;
					} // for
				} // for
				
				//new distance array to hold distances between neuron and i and j in the grid
				double[][] distances = new double[F_mg1 * F_mg2][F_mg1 * F_mg2];
				//now assign the distances accordingly user helper function
				for (int i = 0; i < F_mg1 * F_mg2; i ++)
				{	
					for (int j = 0; j < F_mg1 * F_mg2; j ++)
					{
						distances[i][j] = distance(positions[j], positions[i]);
						//System.out.println("dist(" + (i + 1) + "," + (j + 1) + ") = " + distances[i][j]);
					} // for
				} // for
				
				return distances;
			} // case 2
			case 3:
			{
				//in three dimensions we have 3 edge lengths
				int F_mg1 = currentSOMEdgeLengthsG[g - 3];
				int F_mg2 = currentSOMEdgeLengthsG[g - 2];
				int F_mg3 = currentSOMEdgeLengthsG[g - 1];
				
				//we will hold the relative position in a string in a from e.g. "a=0, b=0, c=0"
				//we will have F_mg1 * F_mg2 * F_mg3 total positions
				String[] positions = new String[F_mg1 * F_mg2 * F_mg3];
				
				//now assign dimensions for each position and hold them in String array
				int index = 0;
				for (int a = 0; a < F_mg1; a ++)
				{	
					for (int b = 0; b < F_mg2; b ++)
					{
						for (int c = 0; c < F_mg3; c ++)
						{
							//we will use comma as a delimiter, the we will split on that.
							positions[index] = a + "," + b + "," + c;
							index ++;
						} // for
					} // for
				} // for
				
				//new distance array to hold distances between neuron and i and j in the grid
				double[][] distances = new double[F_mg1 * F_mg2 * F_mg3][F_mg1 * F_mg2 * F_mg3];
				//now assign the distances accordingly user helper function
				for (int i = 0; i < F_mg1 * F_mg2 * F_mg3; i ++)
				{	
					for (int j = 0; j < F_mg1 * F_mg2 * F_mg3; j ++)
					{
						distances[i][j] = distance(positions[j], positions[i]);
						//System.out.println("dist(" + (i + 1) + "," + (j + 1) + ") = " + distances[i][j]);
					} // for
				} // for
				
				return distances;
			} // case 3
			case 4:
			{
				//in 4 dimensions we have 4 edge lengths
				int F_mg1 = currentSOMEdgeLengthsG[g - 4];
				int F_mg2 = currentSOMEdgeLengthsG[g - 3];
				int F_mg3 = currentSOMEdgeLengthsG[g - 2];
				int F_mg4 = currentSOMEdgeLengthsG[g - 1];
				
				//we will hold the relative position in a string in a from e.g. "a=0, b=0, c=0, d=0"
				//we will have F_mg1 * F_mg2 * F_mg3 * F_mg4 total positions
				String[] positions = new String[F_mg1 * F_mg2 * F_mg3 * F_mg4];
				
				//now assign dimensions for each position and hold them in String array
				int index = 0;
				for (int a = 0; a < F_mg1; a ++)
				{	
					for (int b = 0; b < F_mg2; b ++)
					{
						for (int c = 0; c < F_mg3; c ++)
						{
							for (int d = 0; d < F_mg4; d ++)
							{
								//we will use comma as a delimiter, the we will split on that.
								positions[index] = a + "," + b + "," + c + "," + d;
								index ++;
							} // for
						} // for
					} // for
				} // for
				
				//new distance array to hold distances between neuron and i and j in the grid
				double[][] distances = new double[F_mg1 * F_mg2 * F_mg3 * F_mg4][F_mg1 * F_mg2 * F_mg3 * F_mg4];
				//now assign the distances accordingly user helper function
				for (int i = 0; i < F_mg1 * F_mg2 * F_mg3 * F_mg4; i ++)
				{	
					for (int j = 0; j < F_mg1 * F_mg2 * F_mg3 * F_mg4; j ++)
					{
						distances[i][j] = distance(positions[j], positions[i]);
						//System.out.println("dist(" + (i + 1) + "," + (j + 1) + ") = " + distances[i][j]);
					} // for
				} // for
				
				return distances;
			} // case 4
			case 5:
			{
				//in 5 dimensions we have 5 edge lengths
				int F_mg1 = currentSOMEdgeLengthsG[g - 5];
				int F_mg2 = currentSOMEdgeLengthsG[g - 4];
				int F_mg3 = currentSOMEdgeLengthsG[g - 3];
				int F_mg4 = currentSOMEdgeLengthsG[g - 2];
				int F_mg5 = currentSOMEdgeLengthsG[g - 1];
				
				//we will hold the relative position in a string in a from e.g. "a=0, b=0, c=0, d=0, e=0"
				//we will have F_mg1 * F_mg2 * F_mg3 * F_mg4 * F_mg5 total positions
				String[] positions = new String[F_mg1 * F_mg2 * F_mg3 * F_mg4 * F_mg5];
				
				//now assign dimensions for each position and hold them in String array
				int index = 0;
				for (int a = 0; a < F_mg1; a ++)
				{	
					for (int b = 0; b < F_mg2; b ++)
					{
						for (int c = 0; c < F_mg3; c ++)
						{
							for (int d = 0; d < F_mg4; d ++)
							{
								for (int e = 0; e < F_mg5; e ++)
								{
									//we will use comma as a delimiter, the we will split on that.
									positions[index] = a + "," + b + "," + c + "," + d + "," + e;
									index ++;
								} // for
							} // for
						} // for
					} // for
				} // for
				
				//new distance array to hold distances between neuron and i and j in the grid
				double[][] distances = new double[F_mg1 * F_mg2 * F_mg3 * F_mg4 * F_mg5][F_mg1 * F_mg2 * F_mg3 * F_mg4 * F_mg5];
				//now assign the distances accordingly user helper function
				for (int i = 0; i < F_mg1 * F_mg2 * F_mg3 * F_mg4 * F_mg5; i ++)
				{	
					for (int j = 0; j < F_mg1 * F_mg2 * F_mg3 * F_mg4 * F_mg5; j ++)
					{
						distances[i][j] = distance(positions[j], positions[i]);
						//System.out.println("dist(" + (i + 1) + "," + (j + 1) + ") = " + distances[i][j]);
					} // for
				} // for
				
				return distances;
			} // case 5
			
			default: return null;
			
		} // switch
		
	} // distancesInGrid
	
	//helper method to calculate distance between neurons on dimensions g >= 2
	private static double distance(String positionJ, String positionI)
	{
		String[] JComponents = positionJ.split(",");
		int[] realJComponents = new int[JComponents.length];
		//now converts strings to integers
		for (int index = 0; index < JComponents.length; index ++)
		{
			realJComponents[index] = Integer.parseInt(JComponents[index]);
		} // for
		
		//now do the same thing for positionI
		String[] IComponents = positionI.split(",");
		int[] realIComponents = new int[IComponents.length];
		//now converts strings to integers
		for (int index = 0; index < IComponents.length; index ++)
		{
			realIComponents[index] = Integer.parseInt(IComponents[index]);
		} // for
		
		//calculate euclidean distance
		
		double sumOfSquaredDifferences = 0;
		
		//normally realIComponents.length is equal to realJComponents.length, so:
		for (int i = 0; i < realJComponents.length; i ++)
		{
			sumOfSquaredDifferences += (realJComponents[i] - realIComponents[i]) * (realJComponents[i] - realIComponents[i]);
		}
		
		return Math.sqrt(sumOfSquaredDifferences);
	} // distance
	
	//helper method to calculate euclidean distance between two vectors. In this case, ||P_X - M_C_K||, where P is pattern number of X and M is SOM number of C_K
	private static double distanceBetweenVectors(int vectorSize_N, double[] vector_X, double[] vector_C_k)
	{
		double sumOfSquaredDifferences = 0;
		
		for (int n = 0; n < vectorSize_N; n ++)
		{
			//sumOfSquaredDifferences is actually  || x - c_k || ^ 2
			sumOfSquaredDifferences += (vector_X[n] - vector_C_k[n]) * (vector_X[n] - vector_C_k[n]);
		} // for
		
		
		return Math.sqrt(sumOfSquaredDifferences);
	} // distanceBetweenVectors
	
	//neighborhood function to calculate the neighborhood between neuron i and j. Then the value of neighborhood function will be used in learning process.
	//I chose the neighborhood function as gaussian bell function. It is calculated using exp( (-1/2) * distance^2 / width^2)
	//And neighborhood function is time dependent
	private static double neighborhoodFunction(double distance, double width, int stepT)
	{
		//calculate the width for step t using helper function widthAtStepT(int t)
		width = widthAtStepT(width, stepT);
		
		return Math.exp((-1 / (double) 2) * distance * distance / (width * width));
	} // neighborhoodFunction
	
	//helper method to achieve time decresiing size/width for gaussian.
	private static double widthAtStepT(double initialWidth, int stepT)
	{
		//new width is achieved by the formula width_0 * exp(-stepT / stepTMax);
		return initialWidth * Math.exp( (-1)  * stepT / (double) stepTMax);
	} // widthAtStepT
	
	//helper method to calculate learning rate in each stepT, where we will go through eta_init to eta_final
	private static double learningRate(int stepT)
	{
		//we will define variable alfa for flexibility
		//double alfa = stepTMax * Math.sqrt(eta_final / eta_init);
		
		//return time dependent learning rate, where it will calculated as alfa^t * eta_init
		//return Math.pow(alfa, stepT) * eta_init;
		
		return eta_init * Math.exp( (-1)  * stepT / (double) stepTMax);
	} // learningRate
	
	//helper method to achieve difference of two vectors which in turn is a vector
	private static double[] differenceOfVectors(int vectorSize_N, double[] vector_X, double[] vector_C_k)
	{
		//vector to be returned after certain calculations
		double[] newVector = new double[vectorSize_N];
		
		//we will do u - v which is actually u + (-v)
		//Recall that -v is a scalar multiple of –1 times v.
		//now add u + (-v)
		for (int n = 0; n < N; n ++)
		{
			newVector[n] = vector_X[n] + ( (-1) * vector_C_k[n]);
		} // for
		
		return newVector;
		
	} // differenceOfVectors
	
	//helper method to find out a length of any given vector
	private static double vectorLength(int vectorSize_N, double[] vector)
	{
		double sum = 0;
		
		//go through each component, square them and add that square to the sum
		for (int n = 0; n < vectorSize_N; n ++)
		{
			sum += vector[n] * vector[n];
		} // for
		
		return Math.sqrt(sum);
	} // vectorLength
	
	//helper method to print out centers, which accepts number of SOMs and center vectors in these SOMs
	private static void printCenters(int M, int vectorSize_N, double[][][] centerVector_Matrix_MKN, boolean isInIteration)
	{
		if (!isInIteration)
		{
			System.out.println("------------------ CENTERS -------------");
		
			for (int m = 0; m < M; m ++)
			{	
				//centerVector_Matrix_InThisSOM will have dimensions noOfNeurons_K * N
				double[][] centerVector_Matrix_InThisSOM = centerVector_Matrix_MKN[m];
			
				//print out chosen centers
				for (int k = 0; k < centerVector_Matrix_InThisSOM.length; k ++)
				{	
					System.out.print("C_m" + (m + 1) + "k" + (k + 1) + " = " + "(");
				
					for (int n = 0; n < vectorSize_N; n ++)
					{
						if (n != vectorSize_N - 1)
							System.out.print(centerVector_Matrix_MKN[m][k][n] + ",");
						else
							System.out.print(centerVector_Matrix_MKN[m][k][n]);
					} // for
				
					System.out.println(")");
				} // for
			
				System.out.println();
				
			} // for	
		
			System.out.println("----------------------------------------");
		} // if
		else // write the results to the learning.txt file
		{
			learningWriter.println("------------------ CENTERS -------------");
		
			for (int m = 0; m < M; m ++)
			{	
				//centerVector_Matrix_InThisSOM will have dimensions noOfNeurons_K * N
				double[][] centerVector_Matrix_InThisSOM = centerVector_Matrix_MKN[m];
			
				//print out chosen centers
				for (int k = 0; k < centerVector_Matrix_InThisSOM.length; k ++)
				{	
					learningWriter.print("C_m" + (m + 1) + "k" + (k + 1) + " = " + "(");
				
					for (int n = 0; n < vectorSize_N; n ++)
					{
						if (n != vectorSize_N - 1)
							learningWriter.print(centerVector_Matrix_MKN[m][k][n] + ",");
						else
							learningWriter.print(centerVector_Matrix_MKN[m][k][n]);
					} // for
				
					learningWriter.println(")");
				} // for
			
				learningWriter.println();
				
			} // for	
		
			learningWriter.println("----------------------------------------");
		} // else
			
	} // printCenters
	
	//helper method to print out widths
	private static void printWidths(int M, double[][] widthsS_MK)
	{
		//print out widths
		System.out.println("\n--------------- WIDTHS ---------------");
		//also write it to the files
		resultWriter.println("\n--------------- WIDTHS ---------------");
		learningWriter.println("\n--------------- WIDTHS ---------------");
		
		
		for (int m = 0; m < M; m ++)
		{
			//current number of neurons
			int currentNoOfNeurons_K = widthsS_MK[m].length;
		
			for (int k = 0; k < currentNoOfNeurons_K; k ++)
			{	
				System.out.println("S_m" + (m + 1) + "k" + (k + 1) + " = " + widthsS_MK[m][k]);
				
				//again write it to the files, because this method will be executed only once
				resultWriter.println("S_m" + (m + 1) + "k" + (k + 1) + " = " + widthsS_MK[m][k]);
				learningWriter.println("S_m" + (m + 1) + "k" + (k + 1) + " = " + widthsS_MK[m][k]);
			} // for
			
			System.out.println();
			
			//add new lines to the files
			resultWriter.println();
			learningWriter.println();
		} // for
		
		System.out.println("----------------------------------------");
		
		//end width printing in the files with a lot of dashes
		resultWriter.println("----------------------------------------");
		learningWriter.println("----------------------------------------");
	} // printWidths
	
} // class MSOMTraining