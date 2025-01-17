from ShallowTree.ShallowTree import ShallowTree
from ShallowTree.RandomThresholdTree import RandomThresholdTree
from ExKMC.Tree import Tree
from sklearn.preprocessing import OneHotEncoder
import os
import pandas as pd
import numpy as np
import time
from scipy.io import savemat

# Set the path to the target folder (can be relative or absolute)
current_directory = os.path.dirname(__file__)
file_path = os.path.join(current_directory, 'Sequence2BinaryData')
data_obj_list = ['activity','aslbu','auslan2','context','epitope','gene',
                  'news','pioneer','question','reuters','robot','skating','unix','webkb']
data_K_list = [2, 7, 10, 5, 2, 2, 5, 3, 2, 4, 2, 7, 4, 3]
ODSs = [item + "_binary" for item in data_obj_list]

SHA_pi = []
IMM_pi = []
Execution_times = np.zeros((14, 2))
MaxDepth = np.zeros((14, 2))
AvgLeafDepth = np.zeros((14, 2))
NumLeaf = np.zeros((14, 2))

num_runs = 50

for i in range(0, 14):
    print(i)
    # Read the file content into a DataFrame
    txt_path = os.path.join(file_path, ODSs[i] + ".txt")
    df = pd.read_csv(txt_path, delimiter='\t', header=None, index_col=False)
    ## Initialize OneHotEncoder with dtype=int to ensure numerical output
    # encoder = OneHotEncoder(dtype=int)
    # df_encoded = encoder.fit_transform(df)
    # data = df_encoded.toarray()  # This is a NumPy array
    data = df.values
    K = data_K_list[i]

    for run in range(num_runs):
        # ExShallow
        start_SHA = time.time()
        SHA = ShallowTree(K)
        # Construct the tree, and return cluster labels
        SHA_pi.append(SHA.fit_predict(data))
        end_SHA = time.time()
        Execution_times[i, 0] = Execution_times[i, 0] + (end_SHA - start_SHA)
        # Max Depth, Average Leaf Depth
        MaxDepth[i, 0] = MaxDepth[i, 0] + SHA._max_depth()
        AvgLeafDepth[i, 0] = AvgLeafDepth[i, 0] + SHA.average_leaf_depth()
        NumLeaf[i, 0] = NumLeaf[i, 0] + SHA.count_leaves()
        # Tree plot saved to filename
        # SHA.plot('eg_SHA_pi')

        # IMM #
        start_IMM = time.time()
        IMM = Tree(k=K)
        # Construct the tree, and return cluster labels
        IMM_pi.append(IMM.fit_predict(data))
        end_IMM = time.time()
        Execution_times[i, 1] = Execution_times[i, 1] + (end_IMM - start_IMM)
        # Max Depth, Average Leaf Depth
        MaxDepth[i, 1] = MaxDepth[i, 1] + IMM._max_depth()
        AvgLeafDepth[i, 1] = AvgLeafDepth[i, 1] + IMM.average_leaf_depth()
        NumLeaf[i, 1] = NumLeaf[i, 1] + IMM.count_leaves()
        # Tree plot saved to filename
        # IMM.plot('eg_IMM_pi')

MaxDepth = MaxDepth/num_runs
Execution_times = Execution_times/num_runs
AvgLeafDepth = AvgLeafDepth/num_runs
NumLeaf = NumLeaf/num_runs

savemat('SHA_pi.mat', {'SHA_pi': SHA_pi})
savemat('IMM_pi.mat', {'IMM_pi': IMM_pi})
savemat('MaxDepth.mat', {'MaxDepth': MaxDepth})
savemat('Execution_times.mat', {'Execution_times': Execution_times})
savemat('AvgLeafDepth.mat', {'AvgLeafDepth': AvgLeafDepth})
savemat('NumLeaf.mat', {'NumLeaf': NumLeaf})