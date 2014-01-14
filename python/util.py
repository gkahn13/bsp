import time
from collections import defaultdict

class Profiler:
    def __init__(self):
        self.startTime = dict()
        self.isRunning = dict()
        self.totalTime = defaultdict(int)

    def start(self, name):
        if (self.isRunning.has_key(name) and not self.isRunning[name]) or (not self.isRunning.has_key(name)):
            self.isRunning[name] = True
            self.startTime[name] = time.time()

    def stop(self, name):
        if (self.isRunning.has_key(name) and self.isRunning[name]):
            self.isRunning[name] = False
            self.totalTime[name] += time.time() - self.startTime[name]

    def allTimes(self):
        return self.totalTime.items()
            
            
        
