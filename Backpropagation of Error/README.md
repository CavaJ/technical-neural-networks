There is only one class file that does the job, code will work on JAVA SE 1.6 or later. 
For the java classs file I have implemented multi layer perzeptron with 4 inputs including BIAS and two outputs. 
These can be easily seen from the original "training.dat" file. In the hidden layer there are 3 neurons. Which makes my network 4 - 3 - 2 MLP.


MlpBackPropagationOfError.class => This just implements three-layer perceptron and back propagation of error using hyperbolic tangent function. 
The y_m values in training file could be -0.7 and 0.8, so hyper tan function is the only choice which is between -1 and 1. When logistic function
is applied, the global error does not converge more than 4.54. The Weights are generated between -2.0 and 2.0 and I have set the seed (N+1)*M
which always generates the same weights. This is useful to reproduce the same thing again and for error handling. I have set different learning
rates between different layers. All weights are changing using delta rule. In every step I have written global error to the file named "learning.curve". 
Then I used gnuplot to depict the convergence of global error. You can check "ConvergenceOfError.png" file. After learning completed testing is done
to evaluate the performance of the network. Useful comments are provided in the source file MlpBackPropagationOfError.java .
To run this program, go to the location where "MlpBackPropagationOfError.java" resides in cmd. And do not forget "training_preprocessed.dat" 
and "training_preprocessed_for_test.dat" file must be in the same directory as "MlpBackPropagationOfError.java". Type the following:

java MlpBackPropagationOfError

After running you can check the "result.txt" file to see the convergernce of global error. 
Program has a textual user interface that will help you on working with the program.


NOTE: I have changed "training.dat" file to read it in more convenient way. Original can be found "traininig.dat" file. Preprocessed data can be found 
"training_preprocessed.dat" file. I have just used "," (comma) as a delimiter and nothing more. This procedure have been also applied to the
"training_preprocessed_for_test.dat" testing file.
