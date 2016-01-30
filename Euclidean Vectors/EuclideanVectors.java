//Class file for the Programming Assignment PA-C

import java.io.*;

public class EuclideanVectors
{
	public static void main(String[] args)
	{
		//number of dimensions
		int d = 1000;
		//number of random vectors
		int P = 10000;
		
		//three PrintWriters for writing distribution of lengths, angles and distances
		PrintWriter lengthWriter = null;
		PrintWriter angleAlfaWriter = null;
		PrintWriter angleWriter = null;
		PrintWriter distanceWriter = null;
		
		try
		{
			lengthWriter = new PrintWriter(new File("lengths.dat"));
			angleAlfaWriter = new PrintWriter(new File("anglesBetweenVectorsAndSpaceDiagonalVector.dat"));
			angleWriter = new PrintWriter(new File("angles.dat"));
			distanceWriter = new PrintWriter(new File("distances.dat"));
		}
		catch(FileNotFoundException ex)
		{
			ex.printStackTrace();
		} // catch
		
		//the array to hold the vectors and their values in each dimension
		double[][] vectors_Pd = new double[P][d];
		
		//length array of vectors that will hold the length of each vector
		double[] lengths = new double[P];
		
		//we will consider vectors as points because, we assume their starting point is at the origin
		//In other words, the vectors are regarded as position vectors with respect to an origin. 
		//initialize them in the interval [0,1]
		for (int i = 0; i < P; i ++)
			for (int j = 0; j < d; j ++)
			{
				vectors_Pd[i][j] = Math.random();
			} // for
		
		//now calculate the length of each vector
		for (int i = 0; i < P; i ++)
		{
			double squaredSum = 0;
			
			for (int j = 0; j < d; j ++)
			{
				squaredSum += vectors_Pd[i][j] * vectors_Pd[i][j];
			} // for
			
			lengths[i] = Math.sqrt(squaredSum);
			
			//write to the file
			lengthWriter.println(lengths[i]);
		} // for
		
		//for equally distributed vectors the excpected value for the length of vectors is the average of lengths
		
		double sumOfLengths = 0;
		
		for (int i = 0; i < P; i ++)
		{
			sumOfLengths += lengths[i];
		} // for
		
		double expectedValueVectorLengths = sumOfLengths / P;
		
		//calculate the variance and standard deviation
		double sumOfLengthMinusAverageSquared = 0;
		for (int i = 0; i < P; i ++)
		{
			sumOfLengthMinusAverageSquared += (lengths[i] - expectedValueVectorLengths) * (lengths[i] - expectedValueVectorLengths);
		} // for
		
		double sampleVariance = sumOfLengthMinusAverageSquared / (P - 1);
		double std = Math.sqrt(sampleVariance);
		
		double sumInChiDistribution = 0;
		
		//now we are going to calculate chi distribution of vector lengths
		for (int i = 0; i < P; i ++)
		{
			sumInChiDistribution += ((lengths[i] - expectedValueVectorLengths) / std) * ((lengths[i] - expectedValueVectorLengths) / std);
		} // for
		
		double chiDistribution = Math.sqrt(sumInChiDistribution);
		
		System.out.println("Expected value for vector lengths: " + expectedValueVectorLengths);
		System.out.println("Chi Distribution for vector lengths: " + chiDistribution);
		
		
		//------------------------------------ ANGLES ----------------------------------------------
		//We are asked to calculate the angle between each vector and space diagonal Vector
		//First we need calculate spacde diagonal position vector
		double[] spaceDiagonalVector = new double[d];
		
		//For two position vectors A(x1, y1) and B(x2, y2)the diagonal vector calculated as M( (x1 + x2) / 2, (y1 + y2) / 2) ).
		for (int j = 0; j < d; j ++)
		{
			double sum = 0;
			
			for (int i = 0; i < P; i ++)
			{
				sum += vectors_Pd[i][j];
			} // for
			
			//assign the average j'th cooradinate to spaceDiagonalVector j'th coordinate
			spaceDiagonalVector[j] = sum / P;
		} // for
		
		//Let's calculate the length of space diagonal vector
		double sumOfSquares = 0;
		
		for (int j = 0; j < d; j ++)
		{
				sumOfSquares += spaceDiagonalVector[j] * spaceDiagonalVector[j];
		} // for
		
		double lengthOfSpaceDiagonalVector = Math.sqrt(sumOfSquares);
		
		//array to hold angles between each vector and space diagonal vector
		double[] angles = new double[P];
		
		//now calculate the angle in radians between the random vectors and space diagonal vector
		//cos(theta) is dot product of two vectors divided multiplication of their length
		for (int i = 0; i < P; i ++)
		{
			double sumInNominator = 0;
			
			for (int j = 0; j < d; j ++)
			{
				sumInNominator += vectors_Pd[i][j] * spaceDiagonalVector[j];
			} // for
			
			double multiplicationInDenominator = lengths[i] * lengthOfSpaceDiagonalVector;
			
			//we will arccosine to find the angle, it will be in radians
			angles[i] = Math.acos(sumInNominator / multiplicationInDenominator);
			
			//write to the file
			angleAlfaWriter.println(angles[i]);
		} // for
		
		//now it is time to calculate the expected value and chi distribution of the angles
		double sumOfTheAngles = 0;
		
		for (int i = 0; i < P; i ++)
		{
			sumOfTheAngles += angles[i];
		} // for
		
		double expectedValueOfTheAngles = sumOfTheAngles / P;
		
		//calculate the variance and standard deviation
		double sumOfTheAngleMinusAverageSquared = 0;
		for (int i = 0; i < P; i ++)
		{
			sumOfTheAngleMinusAverageSquared += (angles[i] - expectedValueOfTheAngles) * (angles[i] - expectedValueOfTheAngles);
		} // for
		
		double sampleVarianceOfTheAngles = sumOfTheAngleMinusAverageSquared / (P - 1);
		double stdOfTheAngles = Math.sqrt(sampleVarianceOfTheAngles);
		
		double sumInChiDistributionOfTheAngles = 0;
		
		//now we are going to calculate chi distribution of vector lengths
		for (int i = 0; i < P; i ++)
		{
			sumInChiDistributionOfTheAngles += ((angles[i] - expectedValueOfTheAngles) / stdOfTheAngles) 
												* ((angles[i] - expectedValueOfTheAngles) / stdOfTheAngles);
		} // for
		
		double chiDistributionOfTheAngles = Math.sqrt(sumInChiDistributionOfTheAngles);
		
		System.out.println("Ex. value for angles between each vector and diagonal vector: " + expectedValueOfTheAngles);
		System.out.println("Chi Dist. for angles between each vector and diagonal vector: " + chiDistributionOfTheAngles);
		
		
		//--------------------- In the below code we have calculated angle between each vector. ---------------------------------
		//we define P_prime to speed up calculation time. We will take it 2000 rather than 10000.
		
		int P_prime = 2000;
		
		//The array to hold arc cosine values of angles between two vectors. The values will be in radians
		//The number of angles will n * (n - 1) / 2
		//will describe an angle between vector i and j
		double[][] anglesInRadians = new double[P_prime - 1][P_prime];
		
		for (int i = 0; i < P_prime - 1; i ++)
		{
			for (int j = i + 1; j < P_prime; j ++)
			{
				//nominator of the fraction
				double sumInNominator = 0;
				
				for (int k = 0; k < d; k ++)
				{
					sumInNominator += vectors_Pd[i][k] * vectors_Pd[j][k];
				} // for
				
				//denominator will be multiplication of length for 2 vectors
				double denominator = lengths[i] * lengths[j];
				
				//cos(theta) = sumInNominator / denominator, theta will be acos(sumInNominator / denominator)
				anglesInRadians[i][j] = Math.acos(sumInNominator / denominator);
				
				//write to the file
				angleWriter.println(anglesInRadians[i][j]);
			} // for
		} // for
		
		
		int noOfAngles = P_prime * (P_prime - 1) / 2;
		double sumOfAngles = 0;
		
		for (int i = 0; i < P_prime - 1; i ++)
		{
			for (int j = i + 1; j < P_prime; j ++)
			{
				sumOfAngles += anglesInRadians[i][j];
			} // for
		} // for
		
		double expectedValueOfAngles = sumOfAngles / noOfAngles;
		
		//calculate the variance and standard deviation
		double sumOfAngleMinusAverageSquared = 0;
		
		for (int i = 0; i < P_prime - 1; i ++)
		{
			for (int j = i + 1; j < P_prime; j ++)
			{
				sumOfAngleMinusAverageSquared += (anglesInRadians[i][j] - expectedValueOfAngles) * (anglesInRadians[i][j] - expectedValueOfAngles);
			} // for
		} // for
		
		double sampleVarianceOfAngles = sumOfAngleMinusAverageSquared / (noOfAngles - 1);
		double stdOfAngles = Math.sqrt(sampleVarianceOfAngles);
		
		double sumInChiDistributionOfAngles = 0;
		
		//now we are going to calculate chi distribution of angles
		for (int i = 0; i < P_prime - 1; i ++)
		{
			for (int j = i + 1; j < P_prime; j ++)
			{
				sumInChiDistributionOfAngles += ((anglesInRadians[i][j] - expectedValueOfAngles) / stdOfAngles) 
											  * ((anglesInRadians[i][j] - expectedValueOfAngles) / stdOfAngles);
			} // for
		} // for
		
		double chiDistributionOfAngles = Math.sqrt(sumInChiDistributionOfAngles);
		
		System.out.println("Expected value for angles: " + expectedValueOfAngles);
		System.out.println("Chi Distribution for angles: " + chiDistributionOfAngles);
		
		//-------------------------------------- DISTANCES -------------------------------------------------
		
		//array to hold the distances. Number of distances will be P_prime * (P_prime - 1) / 2
		double[][] distances = new double[P_prime - 1][P_prime];
		
		for (int i = 0; i < P_prime - 1; i ++)
		{
			for (int j = i + 1; j < P_prime; j ++)
			{
				//sum of squared differences for each coordinate of this vector
				double sumOfSquaredDifferences = 0;
				
				//if have two vectos A(x1, y1) and B(x2, y2). The sum will be like this: (x1-x2)^2 + (y1 - y2)^2
				for (int k = 0; k < d; k ++)
				{
					sumOfSquaredDifferences += (vectors_Pd[i][k] - vectors_Pd[j][k]) * (vectors_Pd[i][k] - vectors_Pd[j][k]);
				} // for
				
				distances[i][j] = Math.sqrt(sumOfSquaredDifferences);
				
				//write to the file
				distanceWriter.println(distances[i][j]);
			} // for
		} // for
		
		/*
		System.out.println("Number of Distances Assumed: " + P_prime * (P_prime - 1) /2);
		int noOfDistances = 0;
		
		for (int i = 0; i < P_prime - 1; i ++)
		{
			for (int j = i + 1; j < P_prime; j ++)
			{
				System.out.println("Distance " + (noOfDistances + 1) + ": " + distances[i][j]);
				noOfDistances ++;
			} // for
		} // for
		
		System.out.println("Number of Distances After Printing: " + noOfDistances);
		*/
		
		
		int noOfDistances = P_prime * (P_prime - 1) / 2;
		double sumOfDistances = 0;
		
		for (int i = 0; i < P_prime - 1; i ++)
		{
			for (int j = i + 1; j < P_prime; j ++)
			{
				sumOfDistances += distances[i][j];
			} // for
		} // for
		
		double expectedValueOfDistances = sumOfDistances / noOfDistances;
		
		//calculate the variance and standard deviation
		double sumOfDistanceMinusAverageSquared = 0;
		
		for (int i = 0; i < P_prime - 1; i ++)
		{
			for (int j = i + 1; j < P_prime; j ++)
			{
				sumOfDistanceMinusAverageSquared += (distances[i][j] - expectedValueOfDistances) * (distances[i][j] - expectedValueOfDistances);
			} // for
		} // for
		
		double sampleVarianceOfDistances = sumOfDistanceMinusAverageSquared / (noOfDistances - 1);
		double stdOfDistances = Math.sqrt(sampleVarianceOfDistances);
		
		double sumInChiDistributionOfDistances = 0;
		
		//now we are going to calculate chi distribution of distances
		for (int i = 0; i < P_prime - 1; i ++)
		{
			for (int j = i + 1; j < P_prime; j ++)
			{
				sumInChiDistributionOfDistances += ((distances[i][j] - expectedValueOfDistances) / stdOfDistances) 
											  * ((distances[i][j] - expectedValueOfDistances) / stdOfDistances);
			} // for
		} // for
		
		double chiDistributionOfDistances = Math.sqrt(sumInChiDistributionOfDistances);
		
		System.out.println("Expected value for distances: " + expectedValueOfDistances);
		System.out.println("Chi Distribution for distances: " + chiDistributionOfDistances);
		
		//close all writers to finish IO operation
		lengthWriter.close();
		angleAlfaWriter.close();
		angleWriter.close();
		distanceWriter.close();
		
	} // main
} // class EuclideanVectors