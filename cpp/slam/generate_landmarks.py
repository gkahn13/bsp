import sys
import numpy as np
import random
import time

import IPython

MAX_RANGE = 5.
NUM_LANDMARKS = 3 # must be .5*len(main_landmarks)
main_landmarks = []

main_landmarks += [30, -10]
main_landmarks += [70, 12.5]
main_landmarks += [20, 10]

def generate_landmarks(num_landmarks, num_examples):
    if num_landmarks < NUM_LANDMARKS:
        return
    
    
    f = open('landmarks.txt','w')
    f.write(time.asctime()+'\n')
    for example in xrange(num_examples):
        landmarks = [l + random.uniform(-MAX_RANGE/2, MAX_RANGE/2) for l in main_landmarks]
        for i in xrange(num_landmarks - NUM_LANDMARKS):
            #curr_index = i % len(main_landmarks)
            curr_main_landmark = random.randint(0, NUM_LANDMARKS-1)
            new_pos_x = main_landmarks[2*curr_main_landmark] + random.uniform(-MAX_RANGE/2, MAX_RANGE/2)
            new_pos_y = main_landmarks[2*curr_main_landmark+1] + random.uniform(-MAX_RANGE/2, MAX_RANGE/2)
            landmarks += [new_pos_x, new_pos_y]
            
        f.write(' '.join([str(l) for l in landmarks]) + '\n')
        
    f.close()
        
        

if __name__ == '__main__':
    num_landmarks = 3 if len(sys.argv) < 2 else int(sys.argv[1])
    num_examples = 10 if len(sys.argv) < 3 else int(sys.argv[2])
    generate_landmarks(num_landmarks, num_examples)

