# generateMCA
The project is used to generate mixed covering arrays (MCAs) to test configurable software.

The code in branch optimize is the latest. You should specify the strenth of MCA (usually 2~6, set it small if you want to
get MCA quickly), the number of parameters of the tested software, the maximum number of candidate neighbors (usually 10), and
the number of possible values for each parameter in the file testcase.txt and they are separated by spaces. For example:
2 12 10 4 4 4 3 3 2 2 2 2 2 2 2
The result MCAs will be generated in a txt file.
The testing time function can be set in function getRuntime().

The paper
Wang, Y., Zhou, M., Song, X., Gu, M., & Sun, J. (2017). 
Constructing Cost-Aware Functional Test-Suites Using Nested Differential Evolution Algorithm. 
IEEE Transactions on Evolutionary Computation.
describes our work in details.
