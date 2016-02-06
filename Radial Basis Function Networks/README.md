To run the program just type:

java RBFNets

in the command line. Please make sure that trainRBF.pat file is in the same directory as RBFNets class file. Inside the program you can change seeds to get input distributions that can help you find a global optimum. Number of neurons in RBF Layer is adjustable. I implemented every epoch a new random sequence of training patterns. Also I was calculating
global error with respect to 999 examples.

As you run the program result.txt file will be generated automatically. You can use learning.curve file to depict a graph of gradient descent of Error with respect to epochs.
