import sys, termios, tty
import time
import numpy as np

def press_enter_to_continue(name=None):
    if name is None:
        print('\nPress enter to continue')
    else:
        print('\n{0}: Press enter to continue'.format(name))
    raw_input()
    
class Getch:
    @staticmethod
    def getch(block=True):
        fd = sys.stdin.fileno()
        old_settings = termios.tcgetattr(fd)
        try:
            tty.setraw(sys.stdin.fileno())
            ch = sys.stdin.read(1)
        finally:
            termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
        return ch

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
    
def smaller_angle(angle):
    return (angle + np.pi)%(2*np.pi) - np.pi

def closer_angle(x, a, dir=0):
    """                                                                        
    find angle y (==x mod 2*pi) that is close to a                             
    dir == 0: minimize absolute value of difference                            
    dir == 1: y > x                                                            
    dir == 2: y < x                                                            
    """
    if dir == 0:
        return a + smaller_angle(x-a)
    elif dir == 1:
        return a + (x-a)%(2*np.pi)
    elif dir == -1:
        return a + (x-a)%(2*np.pi) - 2*np.pi