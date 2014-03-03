import sys
import numpy as np
import random

import IPython

MAX_RANGE = 5.
NUM_LANDMARKS = 3
main_landmarks = []

main_landmarks += [30, -10]
main_landmarks += [70, 12.5]
main_landmarks += [20, 10]

def generate_landmarks(num_landmarks, num_examples):
    if num_landmarks < NUM_LANDMARKS:
        return
    
    
    f = open('landmarks.txt','w')
    for example in xrange(num_examples):
        landmarks = []
        i = 0
        for i in xrange(2*num_landmarks):
            curr_index = i % len(main_landmarks)
            new_pos = main_landmarks[curr_index] + random.uniform(-MAX_RANGE/2, MAX_RANGE/2)
            landmarks.append(new_pos)
            
        f.write(' '.join([str(l) for l in landmarks]) + '\n')
        
    f.close()
        
        

if __name__ == '__main__':
    num_landmarks = 3 if len(sys.argv) < 2 else int(sys.argv[1])
    num_examples = 10 if len(sys.argv) < 3 else int(sys.argv[2])
    generate_landmarks(num_landmarks, num_examples)

