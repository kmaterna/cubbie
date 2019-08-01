# What happens in numpy linalg if one of your model parameters never gets touched in 
# your linear system, living squarely in the null space? 

import numpy as np 


data = np.array([2, 0, 1, 1]);
matrix = np.array([[1, 1, 0],[1, -1, 0],[2, -1, 0],[3, -1, 0]])

print(np.shape(data));
print(np.shape(matrix));

m = np.linalg.lstsq(matrix,data, 0.01);
print(m);  
