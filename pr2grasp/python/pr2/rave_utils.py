import openravepy as rave

import numpy as np

def transform_from_to(robot, pose_mat_in_ref, ref_link_name, targ_link_name):
    pose_mat_in_ref = np.array(pose_mat_in_ref)
    
    # ref -> world
    if ref_link_name != 'world':
        ref_from_world = robot.GetLink(ref_link_name).GetTransform()
    else:
        ref_from_world = np.eye(4)

    # target -> world
    if targ_link_name != 'world':
        targ_from_world = robot.GetLink(targ_link_name).GetTransform()
    else:
        targ_from_world = np.eye(4)

    # target -> ref
    targ_from_ref = np.dot(np.linalg.inv(targ_from_world), ref_from_world)

    pose_mat_in_targ = np.dot(targ_from_ref, pose_mat_in_ref)
    return np.array(pose_mat_in_targ)

def transform_relative_pose_for_ik(manip, matrix4, ref_frame, targ_frame):
    robot = manip.GetRobot()

    if ref_frame == "world":
        worldFromRef = np.eye(4)
    else:
        ref = robot.GetLink(ref_frame)
        worldFromRef = ref.GetTransform()

    if targ_frame == "end_effector":        
        targFromEE = np.eye(4)
    else:
        targ = robot.GetLink(targ_frame)
        worldFromTarg = targ.GetTransform()
        worldFromEE = manip.GetEndEffectorTransform()    
        targFromEE = np.dot(np.linalg.inv(worldFromTarg), worldFromEE)       

    refFromTarg_new = matrix4
    worldFromEE_new = np.dot(np.dot(worldFromRef, refFromTarg_new), targFromEE)    

    return np.array(worldFromEE_new)

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
    