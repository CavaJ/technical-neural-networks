import java.util.*;

public class TwoLayerPerzeptron
{
	public static void main(String[] args)
	{
		Scanner inputScanner = new Scanner(System.in);
		
		System.out.println("\n*** Welcome to the Two Layer Perzeptron Implementation ***");
		System.out.println("*** You will enter input vector size \"N\" as an integer which must be less than 101 ***");
		System.out.println("*** You will also enter output vector size \"M\" as an integer which must be less than 30 ***");
		System.out.println("*** Respective weights will be generated automatically between -0.5 and 0.5 ***");
		
		//input vector of size N
		System.out.print("\nPlease enter input vector size N: ");
		int N = inputScanner.nextInt();
		
		//output vector of size M
		System.out.print("Now enter output vector size M: ");
		int M = inputScanner.nextInt();
		
		//size will be N + 1 because of a BIAS weight
		double[] inputVector = new double[N + 1];
		
		//BIAS weight is always 1
		inputVector[0] = 1;
		
		System.out.println("Now enter input values!");
		System.out.println("\nx0 = 1");
		
		for (int i = 1; i < N + 1; i ++)
		{
			System.out.print("x" + i + " = ");
			inputVector[i] = inputScanner.nextDouble();
		} // for
		
		//Generate weights, we will use random number generator
		Random random = new Random();
		
		//first create our matrix for weights, M will start from zero where e.g. y0 will map y1 in real world
		double[][] weightMatrix = new double[N + 1][M];
		
		for (int i = 0; i < N + 1; i ++)
			for (int j = 0; j < M; j ++)
			{
				//random.nextDouble() - 0.5 will guarantee that random number will be between -0.5 and 0.5
				//We also reduce decimal points to 1
				weightMatrix[i][j] = Double.parseDouble(String.format("%.2f", random.nextDouble() - 0.5));
			} // for
				
		
		//Generate wieghts
		System.out.println("\nGenerated Weights:\n");
		
		for (int i = 0; i < N + 1; i ++)
			for (int j = 0; j < M; j ++)
			{
				System.out.println("w" + i + "" + (j + 1) + " = " + weightMatrix[i][j]);
			} // for
		
		//create net*m vector to calculate result of the weighted sum
		double[] net_m = new double[M];
		
		//ym vector map for non-linear transfer function, will detect if ym is 0 or 1
		//we will say ym = 1 if net_m >= 0.0 else ym = 0
		int[] y_m = new int[M];
			
		double sum = 0;
		
		for (int j = 0; j < M; j ++)
		{
			for (int i = 0; i < N + 1; i ++)
			{
				sum = sum + inputVector[i] * weightMatrix[i][j];
			} // for
			
			net_m[j] = sum;
			
			//now define y_m
			y_m[j] = net_m[j] >= 0 ? 1 : 0;
			
			sum = 0;
		} // for
		
		System.out.println("\nWe will say y_m = 1 if net_m >= 0.0 else y_m = 0 !\nOutput:\n");
		
		//Print the result of ym's
		for (int j = 0; j < M; j ++)
		{ 
			System.out.println("net*" + (j + 1) + " = " + net_m[j]);
			System.out.println("y" + (j + 1) + " = " + y_m[j]);
		} // for
		
	} // main
} // class TwoLayerPerceptron