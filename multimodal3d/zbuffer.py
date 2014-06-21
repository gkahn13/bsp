import time
import random
import colorsys
import numpy as np

import openravepy as rave
import roslib
roslib.load_manifest('tfx')
import tfx

import geometry3d
import utils

import IPython

handles = list()
    
class FOV:
    def __init__(self, robot, KK, height, width, F, max_range):
        """
        @param KK        - camera intrinsics matrix
        @param height    - camera pixel height
        @param width     - camera pixel width
        @param F         - focal distance (meters) 
        @param max_range - max range of sensor
        """
        self.robot = robot
        self.KK = KK
        self.height = height
        self.width = width
        self.F = F
        self.max_range = max_range
        
        f = KK[0,0]
        self.H = F*(height/f)
        self.W = F*(width/f)
        
    def directions(self, origin, subsample=.1):
        """
        Returns rays that emanate from the origin through the image plane
        (in 'world' frame)
        
        @param origin     - tfx.pose of camera origin
        @param subsample  - resolution subsampling in width/height
        """
        h_sub = int(subsample*self.height)
        w_sub = int(subsample*self.width)
        N = w_sub*h_sub
        
        height_offsets = np.linspace(-self.H/2.0, self.H/2.0, h_sub)
        width_offsets = np.linspace(-self.W/2.0, self.W/2.0, w_sub)
        
        height_grid, width_grid = np.meshgrid(height_offsets, width_offsets)
        
        height_grid_vec = height_grid.reshape((N,1))
        width_grid_vec = width_grid.reshape((N,1))
        z_vec = np.zeros((N,1))
        
        offsets = np.hstack((width_grid_vec, height_grid_vec, z_vec))
        
        points_cam = (self.max_range/self.F)*(np.tile(np.array([0,0,self.F]), (N,1)) + offsets)
        
        ref_from_world = self.robot.GetLink(origin.frame).GetTransform()
        
        directions = np.zeros((N,3))
        
        origin_world_pos = utils.openraveTransformFromTo(self.robot, origin.matrix, origin.frame, 'world')[0:3,3]
        global handles
        for i in xrange(N):
            p_cam = origin + points_cam[i,:]
            p_world_pos = np.dot(ref_from_world, np.array(p_cam.matrix))[0:3,3]
            #p_world_pos = utils.openraveTransformFromTo(self.robot, p_cam.matrix, p_cam.frame, 'world')[0:3,3]
 
            direction = np.array(p_world_pos) - np.array(origin_world_pos)
            directions[i,:] = direction
 
        
        return directions
    
    def get_beams(self, origin, subsample=.1):
        """
        Given the origin, returns a list of beams defining the FOV
        
        @return beams - ndarray of Beam objects
        """
        dirs = self.directions(origin, subsample)
        
        origin_world = tfx.pose(utils.openraveTransformFromTo(self.robot, origin.matrix, origin.frame, 'world'), frame='world')
        origin_world_pos = origin_world.position.array
        rays = np.hstack((np.tile(origin_world_pos, (dirs.shape[0],1)), dirs))
        
        env = self.robot.GetEnv()
        with env:
            is_hits, hits = env.CheckCollisionRays(rays, None)
            
        h_sub = int(subsample*self.height)
        w_sub = int(subsample*self.width)
        N = h_sub*w_sub
        
        is_hits = is_hits.reshape((w_sub,h_sub))
        hits = hits.reshape((w_sub,h_sub,6))
        dirs = dirs.reshape((w_sub,h_sub,3))
        zpoints = np.zeros((w_sub,h_sub,3))
        
        global handles
        for i in xrange(h_sub):
            for j in xrange(w_sub):
                if is_hits[j,i]:
                    zpoints[j,i,:] = hits[j,i,:3]
                else:
                    zpoints[j,i,:] = origin_world_pos + dirs[j,i,:]
                #handles.append(utils.plot_point(env, zpoints[j,i,:]))
                #print('(i,j): ({0},{1})'.format(i,j))
                #raw_input()
        
        beams = np.empty((h_sub-1,w_sub-1), dtype=object)
        for i in xrange(h_sub-1):
            for j in xrange(w_sub-1):
                beams[i,j] = geometry3d.Beam(origin_world_pos,
                                             zpoints[j,i+1,:],
                                             zpoints[j,i,:],
                                             zpoints[j+1,i,:],
                                             zpoints[j+1,i+1,:])
                #beams[i,j].plot(env)
                #raw_input()
                #geometry3d.handles = list()
        
        return beams
    
    def get_border(self, beams):
        # border == list of triangles
        border = list()
        
        rows, cols = beams.shape
        # deal with left and right border columns
        for i in xrange(rows):
            border += beams[i,0].get_side('left')
            border += beams[i,cols-1].get_side('right')
            
            """
            if i < rows-1:
                border += [Triangle(beams[i,0].c, beams[i,0].d, beams[i+1,0].b),
                           Triangle(beams[i,0].d, beams[i+1,0].b, beams[i+1,0].a)]
                border += [Triangle(beams[i,cols-1].c, beams[i,cols-1].d, beams[i+1,cols-1].b),
                           Triangle(beams[i,cols-1].d, beams[i+1,cols-1].b, beams[i+1,cols-1].a)]
            """
                
        # deal with top and bottom border rows
        for j in xrange(cols):
            border += beams[0,j].get_side('top')
            border += beams[rows-1,j].get_side('bottom')
            
            """
            if j < cols-1:
                border += [Triangle(beams[0,j].a, beams[0,j].d, beams[0,j+1].b),
                           Triangle(beams[0,j].d, beams[0,j+1].b, beams[0,j+1].c)]
                border += [Triangle(beams[rows-1,j].a, beams[rows-1,j].d, beams[rows-1,j+1].b),
                           Triangle(beams[rows-1,j].d, beams[rows-1,j+1].b, beams[rows-1,j+1].c)]
            """
                
        # connect with left and top
        for i in xrange(rows):
            for j in xrange(cols):
                # left
                if i > 0:
                    border += [geometry3d.Triangle(beams[i-1,j].a, beams[i-1,j].d, beams[i,j].b),
                               geometry3d.Triangle(beams[i-1,j].b, beams[i,j].b, beams[i,j].c)]
                    
                # top
                if j > 0:
                    border += [geometry3d.Triangle(beams[i,j-1].c, beams[i,j-1].d, beams[i,j].b),
                               geometry3d.Triangle(beams[i,j-1].b, beams[i,j].b, beams[i,j].a)]
            
        pruned_border = list()
        for tri in border:
            if tri.area() > geometry3d.epsilon:
                pruned_border.append(tri)
                
        return pruned_border
        
    def signed_distance(self, p, beams, border):
        is_inside = False
        beams_vec = beams.reshape((np.prod(beams.shape),))
        for beam in beams_vec:
            if beam.is_inside(p):
                is_inside = True
                break
        sd_sign = -1 if is_inside else 1
            
        sd = min([tri.distance_to(p) for tri in border])
        
        return sd_sign*sd
        
    def plot(self, beams):
        env = self.robot.GetEnv()
        
        rows, cols = beams.shape
        for i in xrange(rows):
            for j in xrange(cols):
                beams[i,j].plot(env, with_sides=False)
            
        
def random_within(lower, upper):
    return random.random()*(upper - lower) + lower     
        
def test_FOV(M=10):
    env = rave.Environment()
    env.Load('envs/pr2-test.env.xml')
    env.SetViewer('qtcoin')
    env.GetViewer().SendCommand('SetFiguresInCamera 1') # also shows the figures in the image
    time.sleep(1)
    robot = env.GetRobots()[0]
    
    cam = robot.GetAttachedSensor('head_cam').GetSensor()
    type = rave.Sensor.Type.Camera
    cam_geom = cam.GetSensorGeometry(type)
    
    height, width, _ = cam_geom.imagedata.shape
    F = .01 # real focal length in meters
    max_range = 5.
    
    fov = FOV(robot, cam_geom.KK, height, width, F, max_range)
        
    cam_pose = tfx.pose([0,0,0.05], frame='wide_stereo_optical_frame')
    
    start = time.time()
    beams = fov.get_beams(cam_pose)
    print('beams time: {0}'.format(time.time()-start))
    
    start = time.time()
    border = fov.get_border(beams)
    print('border time: {0}'.format(time.time()-start))
    
    table = env.GetKinBody('table')
    base = table.GetLink('base')
    extents = base.Geometry.GetBoxExtents(base.GetGeometries()[0])
    
    table_pos = tfx.pose(table.GetTransform()).position
    # assume table has orientation np.eye(3)
    x_min, x_max = table_pos.x - extents[0], table_pos.x + extents[0]
    y_min, y_max = table_pos.y - extents[1], table_pos.y + extents[1]
    z_min, z_max = table_pos.z + extents[2], table_pos.z + extents[2] + .2
    
    particles = list()
    for i in xrange(M):
        x = random_within(x_min, x_max)
        y = random_within(y_min, y_max)
        z = random_within(z_min, z_max)
        particles.append(np.array([x,y,z]))
        
    signed_distances = list()
    for i in xrange(M):
        print(i)
        signed_distances.append(fov.signed_distance(particles[i], beams, border))
        
    min_sd = min(signed_distances)
    sd_hue = np.array(signed_distances)+abs(min_sd)
    sd_hue = (1/3.0)*(sd_hue/max(sd_hue))
    for p, h in zip(particles, sd_hue):
        rgb = colorsys.hsv_to_rgb(h, 1, 1)
        handles.append(utils.plot_point(env, p, color=rgb))
    
    #fov.plot(beams)
    for tri in border:
        tri.plot(env)
    
    IPython.embed()
    """
    start = time.time()
    directions = fov.directions(cam_pose)
    print('Time: {0}'.format(time.time()-start))
    
    cam_pose_world = tfx.pose(utils.openraveTransformFromTo(robot, cam_pose.matrix, cam_pose.frame, 'world'))
    closest_collisions(env, cam_pose_world.position.array, directions, plot=True)
    IPython.embed()
    """

def test_image_rays():
    env = rave.Environment()
    env.Load('envs/pr2-test.env.xml')
    env.SetViewer('qtcoin')
    env.GetViewer().SendCommand('SetFiguresInCamera 1') # also shows the figures in the image
    time.sleep(1)
    robot = env.GetRobots()[0]
    
    cam = robot.GetAttachedSensor('head_cam').GetSensor()
    type = rave.Sensor.Type.Camera
    cam_geom = cam.GetSensorGeometry(type)
    
    depth = robot.GetAttachedSensor('head_depth').GetSensor()
    type = rave.Sensor.Type.Laser
    depth_geom = depth.GetSensorGeometry(type)
    
    #cam.Configure(rave.Sensor.ConfigureCommand.PowerOn)
    #cam.Configure(rave.Sensor.ConfigureCommand.RenderDataOn)
    
    #cam_pose = tfx.pose(cam.GetTransform())
    #cam_pose.position.z += .32
    
    cam_pose = tfx.pose([0,0,0.05], frame='wide_stereo_optical_frame')
    cam_pose_world = tfx.pose(utils.openraveTransformFromTo(robot, cam_pose.matrix, cam_pose.frame, 'world'))
    img_plane_center = cam_pose + [0, 0, .01]
    
    global handles
    img_plane_world = tfx.pose(utils.openraveTransformFromTo(robot, img_plane_center.matrix, cam_pose.frame, 'world'))
    #handles.append(utils.plot_point(env, img_plane_world.position.array, size=.0005))
    
    height, width, _ = cam_geom.imagedata.shape
    f = cam_geom.KK[0,0]
    F = .01 # real focal length in meters
    
    W = F*(width/f)
    H = F*(height/f)
    
    width_offsets = np.linspace(-W/2.0, W/2.0, 64)
    height_offsets = np.linspace(-H/2.0, H/2.0, 48)
    
    directions = np.zeros((len(width_offsets)*len(height_offsets), 3))
    
    index = 0
    for w_offset in width_offsets:
        for h_offset in height_offsets:
            p = img_plane_center + [w_offset, h_offset, 0]
            p_world = tfx.pose(utils.openraveTransformFromTo(robot, p.matrix, p.frame, 'world'))
            direction = (p_world.position.array - cam_pose_world.position.array)
            direction = 5 * direction/np.linalg.norm(direction)
            directions[index,:] = direction
            index += 1
            
            #closest_collision(env, cam_pose_world.position.array, direction, plot=False)
            #handles.append(utils.plot_point(env, p_world.position.array, size=.0001))
    start_time = time.time()
    closest_collisions(env, cam_pose_world.position.array, directions, plot=False)
    print('Total time: {0}'.format(time.time() - start_time))
    
    IPython.embed()
    rave.RaveDestroy()

if __name__ == '__main__':
    #test_collision()
    test_FOV()
    #test_image_rays()
    
    
    
    
    
    
    
    
    
    
    
    
    
"""
OLD
"""

def closest_collision(env, origin, direction, plot=False):
    rays = np.hstack((origin, direction))
    with env:
        is_hits, hits = env.CheckCollisionRays(np.array([rays]), None)
        
    closest_hit, closest_hit_dist = None, np.inf
    for i in xrange(hits.shape[0]):
        if is_hits[i] and np.linalg.norm(hits[i,:3] - origin) < closest_hit_dist:
            closest_hit = hits[i,:3]
            closest_hit_dist = np.linalg.norm(hits[i,:3] - origin)
            
    if plot:
        global handles
        handles += utils.plot_segment(env, origin, origin+direction)
        if closest_hit is not None:
            handles.append(utils.plot_point(env, closest_hit))
            
    return closest_hit

def closest_collisions(env, origin, directions, plot=False):
    rays = np.hstack((np.tile(origin, (directions.shape[0],1)), directions))
    
    with env:
        is_hits, hits = env.CheckCollisionRays(rays, None)
        
    if plot:
        global handles
        for i in xrange(len(rays)):
            handles += utils.plot_segment(env, origin, origin + directions[i,:])
        
        hits_list = [hits[i,:] for i in xrange(len(is_hits)) if is_hits[i]]
        for hit in hits_list:
            handles.append(utils.plot_point(env, hit))

def test_collision():
    env = rave.Environment()
    env.Load('envs/pr2-test.env.xml')
    env.SetViewer('qtcoin')
    time.sleep(1)
    robot = env.GetRobots()[0]
    
    cam_pose = tfx.pose(robot.GetSensor('head_cam').GetTransform())
    
    mug = env.GetKinBody('mug')
    mug_pose = tfx.pose(mug.GetTransform())
    
    start, end = cam_pose.position.array, mug_pose.position.array
    origin = start
    direction = (end - start) / np.linalg.norm(end-start)
    
    closest_hit = closest_collision(env, origin, 5*direction, plot=True)
    print(closest_hit)
    
    
    IPython.embed()

