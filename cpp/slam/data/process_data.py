import os
import IPython

from collections import defaultdict
import math
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

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
    slam_types = ['slam-traj', 'slam-belief', 'slam-state', 'slam-control', 'slam-ilqg', 'slam-control-ham']
    
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
    
    def __str__(self):
        if self.slam_type == 'slam-traj':
            return 'Trajectory'
        if self.slam_type == 'slam-ilqg':
            return 'iLQG'
        if self.slam_type == 'slam-belief':
            return 'Full Coll.'
        if self.slam_type == 'slam-state':
            return 'Partial Coll.'
        if self.slam_type == 'slam-control':
            return 'Shooting'
        if self.slam_type == 'slam-control-ham':
            return 'Shooting Ham.'
        
        
            
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
        
        if f is None:
        	return None, None
        
        d = Data()
        for i in xrange(f.num_examples):
            if f.example[i]['failure'] == 0:
                d.add(f.example[i][attr])
                
        return d.mean, d.sd
        
    def getCostStats(self, num_landmarks):
        return self.getStats(num_landmarks, 'sum_cov_trace')
        
    def getTimeStats(self, num_landmarks):
        return self.getStats(num_landmarks, 'total_time')
    
    def getWaypointErrorStats(self, num_landmarks):
        return self.getStats(num_landmarks, 'waypoint_distance_error')
        
    def compareAttr(self, otherFileGroup, num_landmarks, attr):
        f_self = self.getFileWithLandmark(num_landmarks)
        f_other = otherFileGroup.getFileWithLandmark(num_landmarks)
        
        if f_self is None or f_other is None:
        	return None, None
        
        d = Data()
        for i in xrange(f_self.num_examples):
            if f_self.example[i]['failure'] == 0 and f_other.example[i]['failure'] == 0:
                d.add(f_self.example[i][attr] / f_other.example[i][attr])
                
        return d.mean, d.sd
        
    def compareCost(self, otherFileGroup, num_landmarks):
        return self.compareAttr(otherFileGroup, num_landmarks, 'sum_cov_trace')
        
    def compareTime(self, otherFileGroup, num_landmarks):
        return self.compareAttr(otherFileGroup, num_landmarks, 'total_time')    
    
    def compareWaypointError(self, otherFileGroup, num_landmarks):
        return self.compareAttr(otherFileGroup, num_landmarks, 'waypoint_distance_error')
    
    def __str__(self):
        return str(self.files[0])
    
    @staticmethod
    def toLatexTable(fileGroups, landmarks, attr, shiftFactor=1.):
        latex_str = '\hline & ' + ' & '.join([str(fg) for fg in fileGroups]) + ' \\\\ \n'
        for l in landmarks:
            latex_str += '\hline {0} '.format(l)
            for fg in fileGroups:
                mean, sd = fg.getStats(l, attr)
                if mean is None:
                    latex_str += ' & -'
                else:
                    latex_str += ' & {0:.2f} $\pm$ {1:.2f}'.format(shiftFactor*mean, shiftFactor*sd)
            latex_str += ' \\\\ \n'
        latex_str += '\hline'
        return latex_str
        
    
def process_data():
    curr_path = os.path.abspath('.')
    dir_list = os.listdir(curr_path)
    
    data_files = [f for f in dir_list if '.txt' in f]
    
    slam_types = [f.split('.txt')[0].split('_',1)[0] for f in data_files]
    file_times = [f.split('.txt')[0].split('_',1)[1] for f.split('.txt')[0] in data_files]
    
    files = [File(data_file, slam_type, file_time) for data_file, slam_type, file_time in zip(data_files, slam_types, file_times) if slam_type in File.slam_types]
    
    traj_files = [file for file in files if file.slam_type == 'slam-traj']
    ilqg_files = [file for file in files if file.slam_type == 'slam-ilqg' and file.num_landmarks <= 35]
    belief_files = [file for file in files if file.slam_type == 'slam-belief']
    state_files = [file for file in files if file.slam_type == 'slam-state']
    control_files = [file for file in files if file.slam_type == 'slam-control']
    control_ham_files = [file for file in files if file.slam_type == 'slam-control-ham']
    
    trajFG = FileGroup(traj_files)
    ilqgFG = FileGroup(ilqg_files)
    beliefFG = FileGroup(belief_files)
    stateFG = FileGroup(state_files)
    controlFG = FileGroup(control_files)
    controlHamFG = FileGroup(control_ham_files)

    landmarks = [3,4,5,6,10,15,20,25,30,35,40,45,50]
    time_abs_fig = plt.figure()
    time_abs_ax = time_abs_fig.add_subplot(111)
    dist_err_abs_fig = plt.figure()
    dist_err_abs_ax = dist_err_abs_fig.add_subplot(111)
    
    
    
    """
    print([(l, trajFG.getTimeStats(l)) for l in landmarks])
    for l in landmarks:
    	m, s = trajFG.getTimeStats(l)
    	print('{0:.2f} $\pm$ {1:.2f}'.format(m,s))
    return
    """
    print(FileGroup.toLatexTable([trajFG, ilqgFG, beliefFG, stateFG, controlFG, controlHamFG], landmarks, 'total_time', shiftFactor=.001))
    
    print('############ Absolute statistics #########')

    for fg in [controlFG, controlHamFG, stateFG, beliefFG, ilqgFG]:
    	time_abs_avgs, time_abs_sds = [], []
        dist_err_avgs, dist_err_sds = [], []
        for num_landmarks in landmarks:
            cost_avg, cost_sd = fg.getCostStats(num_landmarks)
            time_avg, time_sd = fg.getTimeStats(num_landmarks)
            dist_err_avg, dist_err_sd = fg.getWaypointErrorStats(num_landmarks)
            
            if cost_avg is not None:
            	time_abs_avgs.append(time_avg / 1000.)
            	time_abs_sds.append(time_sd / 1000.)
            	dist_err_avgs.append(dist_err_avg)
             	dist_err_sds.append(dist_err_sd)
            
                print('Number of landmarks: ' + str(num_landmarks))
                print(fg.slam_type)
                print('Cost: {0} +- {1}'.format(cost_avg, cost_sd))
                print('Time: {0} +- {1} ms'.format(time_avg, time_sd))
                print('Waypoint distance error: {0} +- {1} meters'.format(dist_err_avg, dist_err_sd))
                print('')
        
        (_, caps, errlines) = time_abs_ax.errorbar(landmarks[:len(time_abs_avgs)], time_abs_avgs, linewidth=3.0, yerr=time_abs_sds, elinewidth=1.0, capsize=6.0, label=str(fg))
        #for errline in errlines:
        #    errline.set_linestyle('dashed')
        dist_err_abs_ax.errorbar(landmarks[:len(dist_err_avgs)], dist_err_avgs, yerr=dist_err_sds, label=str(fg))
    print('\n')

    
    print('############ Relative statistics #############')
    cost_fig = plt.figure()
    cost_ax = cost_fig.add_subplot(111)
    time_fig = plt.figure()
    time_ax = time_fig.add_subplot(111)
    dist_err_fig = plt.figure()
    dist_err_ax = dist_err_fig.add_subplot(111)
    
    state_comp_times = []
    control_comp_times = []
    for fg in [controlFG, controlHamFG, stateFG, beliefFG, ilqgFG]:
        cost_comp_avgs, cost_comp_sds = [], []
        time_comp_avgs, time_comp_sds = [], []
        dist_err_comp_avgs, dist_err_comp_sds = [], []
        for num_landmarks in landmarks:
            cost_comp_avg, cost_comp_sd = trajFG.compareCost(fg, num_landmarks) # fg.compareCost(trajFG, num_landmarks)
            time_comp_avg, time_comp_sd = trajFG.compareTime(fg, num_landmarks) #fg.compareTime(trajFG, num_landmarks)
            dist_err_comp_avg, dist_err_comp_sd = fg.compareWaypointError(trajFG, num_landmarks)
            
            if cost_comp_avg is not None:
                cost_comp_avgs.append(cost_comp_avg)
                cost_comp_sds.append(cost_comp_sd)
                time_comp_avgs.append(time_comp_avg)
                time_comp_sds.append(time_comp_sd)
                dist_err_comp_avgs.append(dist_err_comp_avg)
                dist_err_comp_sds.append(dist_err_comp_sd)
	            
                print('Number of landmarks: ' + str(num_landmarks))
                print(fg.slam_type + ' compared with ' + trajFG.slam_type)
                print('Cost: {0} +- {1}'.format(cost_comp_avg, cost_comp_sd))
                print('Time: {0} +- {1}'.format(time_comp_avg, time_comp_sd))
                print('Waypoint error: {0} +- {1}'.format(dist_err_comp_avg, dist_err_comp_sd))
                print('')
            
        (_, caps, errlines) = cost_ax.errorbar(landmarks[:len(cost_comp_avgs)], cost_comp_avgs, linewidth=3.0, yerr=cost_comp_sds, elinewidth=1.0, capsize=6.0, label=str(fg))
        #for errline in errlines:
        #    errline.set_linestyle('dashed')
        time_ax.errorbar(landmarks[:len(time_comp_avgs)], time_comp_avgs, yerr=time_comp_sds, label=fg.slam_type)
        dist_err_ax.errorbar(landmarks[:len(dist_err_comp_avgs)], dist_err_comp_avgs, yerr=dist_err_comp_sds, label=str(fg))
            
    # set titles, labels, etc for all the graphs
    for ax in [time_abs_ax, dist_err_abs_ax, cost_ax, time_ax, dist_err_ax]:
        ax.set_xticks(landmarks)
        ax.set_xlabel('Number of landmarks')
        
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc='upper left')
    
    time_abs_ax.set_ylabel('Time (seconds)')
    dist_err_abs_ax.set_ylabel('Waypoint distance error (meters)')
    cost_ax.set_ylabel('Improvement')
    time_ax.set_ylabel('Time factor')
    dist_err_ax.set_ylabel('Waypoint distance error versus trajectory')
    
    time_abs_ax.set_title('Time per number landmarks')
    dist_err_abs_ax.set_title('Waypoint distance error')
    cost_ax.set_title('Cost factor versus trajectory')
    time_ax.set_title('Time factor of belief, state, and control versus trajectory')
    dist_err_ax.set_title('Waypoint distance factor of belief, state, and control versus trajectory')

    def to_percent(y, position):
        # Ignore the passed in position. This has the effect of scaling the default
        # tick locations.
        s = '%2.0f'%(100*y)

        # The percent symbol needs escaping in latex
        if matplotlib.rcParams['text.usetex'] == True:
            return s + r'$\%$'
        else:
            return s + '%' 

    cost_ax.yaxis.set_major_formatter(FuncFormatter(to_percent))
    
    plt.show(block=False)
    raw_input()
            
    

if __name__ == '__main__':
    process_data()
    
