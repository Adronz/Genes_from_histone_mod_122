import numpy as np
import math

#read in data 
def read_fasta(file_path):
    with open(file_path, 'r') as file:
        sequences = []
        sequence = ''
        for line in file:
            if line.startswith('>'):
                continue
            else:
                sequences.append(line.strip())
    return sequences


'''
x : the number of ChIP-seq reads for only mark X in the interval was above a significance
threshold
y : the number of ChIP-seq reads for only mark Y in the interval was above a significance
threshold
z : the number of ChIP-seq reads for both marks X and Y in the interval was above a
significance threshold
n : the number of ChIP-seq reads for neither marks X and Y in the interval was above a
significance threshold
The task is to predict 200bp intervals that are most likely overlapping annotated gene bodies of
protein coding genes.
'''

#create forward-backward algorithm 

def compute_forward(x:list, states:list, emission_prob, transition, alpha_idx):
     fwd = np.zeros((len(states), len(x)))
     #set initial alphas for states 0 and 1
     for i in range(len(states)):
          fwd[i][0] = np.log(1/len(states))
#compute the rest of alpha_t
     for t in range(1,len(x)):
         for j in range(len(states)):
             values = []
             #sum the probabilities of the new alpha in state j coming from either of the previous states * the emission for that state
             for i in range(len(states)):
                 values.append(fwd[i][t-1] + np.log(emission_prob[j][alpha_idx[x[t]]]) + np.log(transition[i][j]))
             fwd[j][t] = np.logaddexp(values[0],values[1])

     return fwd

def compute_backward(x, states, emission_prob, transition, alpha_idx):
    bwd = np.zeros((len(states),len(x)))
    #initialize the last column to -1
    for i in range(len(states)):
        bwd[i][len(x)-1] = 0

    for t in range(len(x)-2,-1,-1):
        for j in range(len(states)):
          values = []
          for k in range(len(states)):
              values.append(bwd[k][t+1] + np.log(transition[j,k]) + np.log(emission_prob[j][alpha_idx[x[t+1]]]))
          
          bwd[j][t] = np.logaddexp(values[0],values[1])
    return bwd
    
        
def forward_backward(x, states, emission_prob,transition, alpha_idx):
    fwd = compute_forward(x, states, emission_prob, transition, alpha_idx)
    bwd = compute_backward(x, states, emission_prob, transition, alpha_idx)
  

    print(fwd[0][400])
    print(bwd[0][(len(x)-400)])


    # print(f_sink)
    soft_decode = np.zeros((len(states),len(x)))
    for i in range(len(x)):
        f_sink = -np.inf  
        for j in range(len(states)):
            soft_decode[j, i] = fwd[j, i] + bwd[j, i]
            f_sink = np.logaddexp(f_sink, soft_decode[j, i])  # Update f_sink with the sum over states
            # print(round(f_sink))
   
        for k in range(len(states)):
            soft_decode[k,i] -= f_sink
            soft_decode[k,i] = np.exp(soft_decode[k,i])

    for i in range(len(x)):
        print(soft_decode[0][i])
    return soft_decode
    


# Example usage
x = read_fasta('/Users/AdrianHanson/Downloads/input.fasta')
states = [0,1]
transition_prob = np.array([[0.9, 0.1], [0.1, 0.8]])
emission_prob = np.array([[0.25, 0.25, 0.1, 0.4], [0.25, 0.25, 0.4, 0.1]])
alphabet = ['x','y','z','n']
alpha_idx = {'x': 0, 'y': 1, 'z': 2,'n': 3}

result = forward_backward(x,states,emission_prob,transition_prob,alpha_idx)


# Assuming result is a 2D NumPy array with shape (2, n)
# Transpose the array to make it (n, 2)
result = result.T

# Add an index to each row
result_with_index = np.hstack((result, np.arange(1, result.shape[0] + 1)[:, None]))

# Sort the rows based on the second column (probability) in descending order
sorted_result = result_with_index[result_with_index[:, 1].argsort()][::-1]

# Extract the top 50000 indices
top_indices = sorted_result[:50000, 2]

# Sort the indices in ascending order
top_indices_sorted = np.sort(top_indices)

# Convert the sorted indices to a list of integers
answers = top_indices_sorted.astype(int).tolist()

# Write the sorted indices to a file
with open("/Users/AdrianHanson/CS 122/predictions.txt", "w") as f:
    for index in top_indices_sorted:
        f.write(f"{index}\n")
