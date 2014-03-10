import os
import IPython

from collections import defaultdict
import math
import numpy as np

import matplotlib
import matplotlib.pyplot as plt

attrs = ['sum_cov_trace','waypoint_distance_error','solve_time','initialization_time']

class Data:
	def __init__(self):
		self.data = []
		
	def add(self, d):
		self.data.append(d)
		
	@property
	def mean(self):
		if len(self.data) > 0:
			return sum(self.data)/float(len(self.data))
		return None
		
	@property
	def sd(self):
		m = self.mean
		if m is None:
			return None
			
		return np.std(np.array(self.data))

class File:
    attrs = ['sum_cov_trace','waypoint_distance_error','solve_time','initialization_time','total_time'] # 'failure'
    slam_types = ['slam-traj', 'slam-belief', 'slam-state', 'slam-control']
    
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
                if attr == File.attrs[0]:
                    self.num_examples += 1
                self.example[self.num_examples-1][attr] += float(val)
                
            for i in xrange(self.num_examples):
                self.example[i]['total_time'] = self.example[i]['solve_time'] + self.example[i]['initialization_time']
                
    def get(self, attr):
        return [self.example[i][attr] for i in xrange(self.num_examples)]
        
            
# assume one file per landmark number
class FileGroup:
	def __init__(self, files):
		self.files = files
		self.slam_type = self.files[0].slam_type
		
	def getFileWithLandmark(self, num_landmarks):
		for f in self.files:
			if f.num_landmarks == num_landmarks:
				return f
		return None
		
	def getStats(self, num_landmarks, attr):
		f = self.getFileWithLandmark(num_landmarks)
		
		d = Data()
		for i in xrange(f.num_examples):
			if f.example[i]['failure'] == 0:
				d.add(f.example[i][attr])
				
		return d.mean, d.sd
		
	def getCostStats(self, num_landmarks):
		return self.getStats(num_landmarks, 'sum_cov_trace')
		
	def getTimeStats(self, num_landmarks):
		return self.getStats(num_landmarks, 'solve_time')
		
	def compareAttr(self, otherFileGroup, num_landmarks, attr):
		f_self = self.getFileWithLandmark(num_landmarks)
		f_other = otherFileGroup.getFileWithLandmark(num_landmarks)
		
		d = Data()
		for i in xrange(f_self.num_examples):
			if f_self.example[i]['failure'] == 0 and f_other.example[i]['failure'] == 0:
				d.add(f_self.example[i][attr] / f_other.example[i][attr])
				
		return d.mean, d.sd
		
	def compareCost(self, otherFileGroup, num_landmarks):
		return self.compareAttr(otherFileGroup, num_landmarks, 'sum_cov_trace')
		
	def compareTime(self, otherFileGroup, num_landmarks):
		return self.compareAttr(otherFileGroup, num_landmarks, 'total_time')	
        
    
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
    traj_files = [file for file in files if file.slam_type == 'slam-traj']
    
    #beliefFG = FileGroup(belief_files)
    stateFG = FileGroup(state_files)
    controlFG = FileGroup(control_files)
    trajFG = FileGroup(traj_files)
    
    landmarks = [3,4,5,6]
    
    print('############ Absolute statistics #########')
    for num_landmarks in landmarks:
    	print('Number of landmarks: ' + str(num_landmarks))
    	for fg in [trajFG, stateFG, controlFG]:
    		cost_avg, cost_sd = fg.getCostStats(num_landmarks)
    		time_avg, time_sd = fg.getTimeStats(num_landmarks)
    		
    		print(fg.slam_type)
    		print('Cost: {0} +- {1}'.format(cost_avg, cost_sd))
    		print('Time: {0} +- {1} ms'.format(time_avg, time_sd))
    print('\n')
    
    print('############ Relative statistics #############')
    cost_fig = plt.figure(1)
    cost_ax = cost_fig.add_subplot(111)
    time_fig = plt.figure(2)
    time_ax = time_fig.add_subplot(111)
    
    state_comp_times = []
    control_comp_times = []
    for fg in [stateFG, controlFG]:
    	cost_comp_avgs, cost_comp_sds = [], []
    	time_comp_avgs, time_comp_sds = [], []
        for num_landmarks in landmarks:
            print('Number of landmarks: ' + str(num_landmarks))
            cost_comp_avg, cost_comp_sd = fg.compareCost(trajFG, num_landmarks)
            time_comp_avg, time_comp_sd = fg.compareTime(trajFG, num_landmarks)
            
            cost_comp_avgs.append(cost_comp_avg)
            cost_comp_sds.append(cost_comp_sd)
            time_comp_avgs.append(time_comp_avg)
            time_comp_sds.append(time_comp_sd)
    		
            print(fg.slam_type + ' compared with ' + trajFG.slam_type)
            print('Cost: {0} +- {1}'.format(cost_comp_avg, cost_comp_sd))
            print('Time: {0} +- {1}'.format(time_comp_avg, time_comp_sd))
            
        cost_ax.errorbar(landmarks, cost_comp_avgs, yerr=cost_comp_sds, elinewidth=2, label=fg.slam_type)
        time_ax.errorbar(landmarks, time_comp_avgs, yerr=time_comp_sds, label=fg.slam_type)
            
    for ax in [cost_ax, time_ax]:
        ax.set_xticks(landmarks)
        ax.set_xlabel('Number of landmarks')
        
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels)
    
    cost_ax.set_ylabel('Cost factor versus trajectory')
    time_ax.set_ylabel('Time factor versus trajectory')
    
    cost_ax.set_title('Cost factor of belief, state, and control versus trajectory')
    time_ax.set_title('Time factor of belief, state, and control versus trajectory')
    
    plt.show(block=False)
    raw_input()
    		
    

if __name__ == '__main__':
    process_data()
    