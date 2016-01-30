import java.util.Random;

public class RandomSequenceGenerator
{
	public static void main(String[] args)
	{
		//number of training patterns
		int P = 1000;
		
		//Assume all initialization is done here.........
		
		
		//Here is learning process through iterations
		int iterationCount = 40000;
		Random rand = new Random();
		
		for (int iteration = 0; iteration < iterationCount; iteration ++)
		{	
			
			//Create a random sequence of indexes, number of Indexes will be equal to number of patterns.
			int[] sequence = new int[P];
			for (int i = 0; i < sequence.length; i ++)
			{
				sequence[i] = i;
			} // for
			
			
			//Now Shuffle the indexes in sequence array
			for (int i = 0; i < sequence.length; i ++)
			{
				//nextInt(int bound) method returns something between 0 and bound, but bound exclusive.
				//generates a random integer between i and sequence.length, but value of sequence.length is exclusive.
				int r = rand.nextInt(sequence.length - i) + i;
				
				//simple swap operation
				int tmp = sequence[r];
				sequence[r] = sequence[i];
				sequence[i] = tmp;
			} // for
			
			//Now go through the patterns. In each iteration we will get different sequence of patterns
			for (int pat = 0; pat < P; pat ++)
			{
				int randomPatternIndex = sequence[pat];
				
				//Use randomPatternIndex to get values of different patterns and do other stuff........
			} // for
			
		} // for
		
	} // main
	
} // class
