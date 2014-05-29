import numpy as np

def openraveTransformFromTo(robot, poseMatInRef, refLinkName, targLinkName):
    poseMatInRef = np.array(poseMatInRef)
    
    # ref -> world
    if refLinkName != 'world':
        refFromWorld = robot.GetLink(refLinkName).GetTransform()
    else:
        refFromWorld = np.eye(4)

    # target -> world
    targFromWorld = robot.GetLink(targLinkName).GetTransform()

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