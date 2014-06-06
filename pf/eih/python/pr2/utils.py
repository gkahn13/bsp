import numpy as np
import scipy
import time

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