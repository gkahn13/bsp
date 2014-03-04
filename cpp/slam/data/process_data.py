import os
import IPython

from collections import defaultdict
import math
import numpy as np

attrs = ['sum_cov_trace','waypoint_distance_error','solve_time','initialization_time']

class File:
    attrs = ['sum_cov_trace','waypoint_distance_error','solve_time','initialization_time']
    slam_types = ['slam-belief', 'slam-state', 'slam-control']
    
    def __init__(self, file_name, slam_type, file_time):
        self.file_name = file_name
        self.slam_type = slam_type
        self.file_time = file_time
        
        # self.example[example_number][attr]
        self.example = defaultdict(lambda: defaultdict(float))
        self.attr_vals = defaultdict(float)
        
        self.num_landmarks = 0
        self.num_examples = 0
        
        self.process()
        
    def process(self):
        f = open(self.file_name,'r')
        with open(self.file_name) as f:
            for line in f:
                if self.num_landmarks == 0:
                    # first line is the number of landmarks
                    self.num_landmarks = int(line)
                    continue
                
                attr, val = line.split(' ')
                if attr == attrs[0]:
                    self.num_examples += 1
                self.example[self.num_examples-1][attr] += float(val)
                
    def get(self, attr):
        return [self.example[i][attr] for i in xrange(self.num_examples)]
    
    def printAverages(self):
        for attr in File.attrs:
            avg = sum([self.example[i][attr] for i in xrange(self.num_examples)])/float(self.num_examples)
            print(attr + ': ' + str(avg))
        print('')
            
    @staticmethod
    def printStatistics(files):
        landmark_numbers = sorted(list(set([f.num_landmarks for f in files])))
        
        for num_landmarks in landmark_numbers:
            files_l = [f for f in files if f.num_landmarks == num_landmarks]
            combined_num_examples = float(sum([f.num_examples for f in files_l]))
            if combined_num_examples > 0:
                print('Number of landmarks: ' + str(num_landmarks))
                for attr in File.attrs:
                    attr_vals = []
                    for f in files_l:
                        attr_vals += f.get(attr)
                    avg = sum(attr_vals) / combined_num_examples
                    print(attr + ' avg: ' + str(avg))
                    # sum(all the data centered)
                    sd = np.sqrt(np.sum(np.power((np.array(attr_vals)-avg),2)) / combined_num_examples)
                    print(attr + ' sd: ' + str(sd))
                print('')
                
        
    
def process_data():
    curr_path = os.path.abspath('.')
    dir_list = os.listdir(curr_path)
    
    data_files = [f for f in dir_list if '.txt' in f]
    
    slam_types = [f.split('.txt')[0].split('_',1)[0] for f in data_files]
    file_times = [f.split('.txt')[0].split('_',1)[1] for f.split('.txt')[0] in data_files]
    
    files = [File(data_file, slam_type, file_time) for data_file, slam_type, file_time in zip(data_files, slam_types, file_times) if slam_type in File.slam_types]
    
    belief_files = [file for file in files if file.slam_type == 'slam-belief']
    state_files = [file for file in files if file.slam_type == 'slam-state']
    control_files = [file for file in files if file.slam_type == 'slam-control']
    
    print('belief_files statistics')
    File.printStatistics(belief_files)
    print('state_files average')
    File.printStatistics(state_files)
    print('control_files average')
    File.printStatistics(control_files)
    
    #IPython.embed()
    

if __name__ == '__main__':
    process_data()
    