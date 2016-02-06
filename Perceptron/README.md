There are three class files that does the job, code will work on JAVA SE 1.6 or later. 
For the 2nd and 3rd class file I have implemented perzeptron with 6 inputs including BIAS and two outputs. 
These can be easily seen from the original "PA-A-train.dat" file.


1) TwoLayerPerzeptron.class => This just implements two-layer perceptron using Step function that can be 0 or 1. 
Useful comments are provided in the source file TwoLayerPerzeptron.java .
To run this program, go to the location where "TwoLayerPerzeptron.java" resides in cmd and type the following:

java TwoLayerPerzeptron

Program is textual user interface that will help you on working with the program.

2) TwoLayerPerzeptronLearning.class => This just implements two-layer perceptron using Logistic function. Weights are generated randomly.
Also initial weights are changing using Widrow-Hoff learning rule. Useful comments are provided in the source file TwoLayerPerzeptronLearning.java .
To run this program, go to the location where "TwoLayerPerzeptronLearning.java" resides in cmd. And do not forget "PA-A-train.dat" must be in the
same directory as "TwoLayerPerzeptronLearning.java". Type the following:

java TwoLayerPerzeptronLearning

After running you can check the "result.txt" file to see the convergernce of global error.

3) TwoLayerPerzeptronLearningWeightsFromFile.class => This just implements two-layer perceptron using Logistic function. 
Weights are taken from file "weights.dat". Also initial weights are changing using Widrow-Hoff learning rule.
Useful comments are provided in the source file TwoLayerPerzeptronLearningWeightsFromFile.java .
To run this program, go to the location where "TwoLayerPerzeptronLearningWeightsFromFile.java" resides in cmd. 
And do not forget "PA-A-train.dat" and "weights.dat" must be in the
same directory as "TwoLayerPerzeptronLearningWeightsFromFile.java". Type the following:

java TwoLayerPerzeptronLearningWeightsFromFile

After running you can check the "result.txt" file to see the convergernce of global error.

NOTE: I have changed "PA-A-train.dat" file to read it in more convenient way. Original can be found "Original PA-A-train.dat" file.
