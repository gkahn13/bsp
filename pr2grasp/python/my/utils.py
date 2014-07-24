import time

class Timeout():
    def __init__(self, timeout_time):
        """
        timeoutTime is integer of how long until times out
        """
        self.timeout_time = timeout_time
    
    def start(self):
        """
        Restarts timeout every time this method is called
        """
        self.end_time = time.time() + self.timeout_time
    
    def has_timed_out(self):
        """
        returns true if time since start method called is
        greater than the current time
        """
        return time.time() > self.end_time 