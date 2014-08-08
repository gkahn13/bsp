import numpy as np
import scipy
import time
import heapq

def openraveTransformFromTo(robot, poseMatInRef, refLinkName, targLinkName):
    poseMatInRef = np.array(poseMatInRef)
    
    # ref -> world
    if refLinkName != 'world':
        refFromWorld = robot.GetLink(refLinkName).GetTransform()
    else:
        refFromWorld = np.eye(4)

    # target -> world
    if targLinkName != 'world':
        targFromWorld = robot.GetLink(targLinkName).GetTransform()
    else:
        targFromWorld = np.eye(4)

    # target -> ref
    targFromRef = np.dot(np.linalg.inv(targFromWorld), refFromWorld)

    poseMatInTarg = np.dot(targFromRef, poseMatInRef)
    return np.array(poseMatInTarg)

def plot_point(env, pos_array, size=.01, color=None):
    color = color if color is not None else np.array([0, 1, 0])
    
    handles = env.plot3(points=pos_array,
                        pointsize=size,
                        colors=color,
                        drawstyle=1)
    return handles

def plot_segment(env, start, end, color=(1,0,0)):
    start = np.array(start)
    end = np.array(end)
    
    h = []
    h.append(env.drawlinestrip(points=np.array([start, end]), linewidth=3.0, colors=np.array([color,color])))
    return h

def plot_transform(env, T, s=0.1):
    """
    Plots transform T in openrave environment.
    S is the length of the axis markers.
    """
    T = np.array(T)
    h = []
    x = T[0:3,0]
    y = T[0:3,1]
    z = T[0:3,2]
    o = T[0:3,3]
    h.append(env.drawlinestrip(points=np.array([o, o+s*x]), linewidth=3.0, colors=np.array([(1,0,0),(1,0,0)])))
    h.append(env.drawlinestrip(points=np.array([o, o+s*y]), linewidth=3.0, colors=np.array(((0,1,0),(0,1,0)))))
    h.append(env.drawlinestrip(points=np.array([o, o+s*z]), linewidth=3.0, colors=np.array(((0,0,1),(0,0,1)))))
    return h

def save_view(env, file_name):
    env.GetViewer().SendCommand('SetFiguresInCamera 1') # also shows the figures in the image
    I = env.GetViewer().GetCameraImage(640,480,  env.GetViewer().GetCameraTransform(),[640,640,320,240])
    scipy.misc.imsave(file_name ,I)
    env.GetViewer().SendCommand('SetFiguresInCamera 0')

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

class Getch:
    @staticmethod
    def getch():
        import sys, tty, termios
        fd = sys.stdin.fileno()
        old_settings = termios.tcgetattr(fd)
        try:
            tty.setraw(sys.stdin.fileno())
            ch = sys.stdin.read(1)
        finally:
            termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
        return ch
    
import heapq
 
 
class PriorityQueue(object):
    """Priority queue based on heap, capable of inserting a new node with
    desired priority, updating the priority of an existing node and deleting
    an abitrary node while keeping invariant"""
 
    def __init__(self, heap=[]):
        """if 'heap' is not empty, make sure it's heapified"""
 
        heapq.heapify(heap)
        self.heap = heap
        self.entry_finder = dict({i[-1]: i for i in heap})
        self.REMOVED = '<remove_marker>'
 
    def insert(self, node, priority=0):
        """'entry_finder' bookkeeps all valid entries, which are bonded in
        'heap'. Changing an entry in either leads to changes in both."""
 
        if node in self.entry_finder:
            self.delete(node)
        entry = [priority, node]
        self.entry_finder[node] = entry
        heapq.heappush(self.heap, entry)
 
    def delete(self, node):
        """Instead of breaking invariant by direct removal of an entry, mark
        the entry as "REMOVED" in 'heap' and remove it from 'entry_finder'.
        Logic in 'pop()' properly takes care of the deleted nodes."""
 
        entry = self.entry_finder.pop(node)
        entry[-1] = self.REMOVED
        return entry[0]
 
    def pop(self):
        """Any popped node marked by "REMOVED" does not return, the deleted
        nodes might be popped or still in heap, either case is fine."""
 
        while self.heap:
            priority, node = heapq.heappop(self.heap)
            if node is not self.REMOVED:
                del self.entry_finder[node]
                return priority, node
        return None, None
 