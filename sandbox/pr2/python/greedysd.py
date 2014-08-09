import sys
import time
import random
import colorsys
import Queue
import numpy as np

import openravepy as rave
import roslib
roslib.load_manifest('tfx')
import tfx

import geometry3d
import utils

import IPython

handles = list()
    
class Camera:
    FOCAL_LENGTH = .01
    MAX_RANGE = 5.0
    
    WIDTH   =    256
    HEIGHT  =      192
    fx  = WIDTH # WIDTH*2.0;
    fy  = HEIGHT # HEIGHT*2.0;
    cx  = WIDTH/2.0 + 0.5 # WIDTH/2.0 + 0.5;
    cy  = HEIGHT/2.0 + 0.5 # HEIGHT/2.0 + 0.5;
    HEIGHT_M = FOCAL_LENGTH*(HEIGHT/fy);
    WIDTH_M = FOCAL_LENGTH*(WIDTH/fx);
    N = WIDTH*HEIGHT
    KK = np.array([[fx, 0, cx],
                   [0, fy, cy],
                   [0, 0, 1]])
    
    W_SUB = 64
    H_SUB = 48
    fx_sub  = W_SUB*1.0 # W_SUB*2.0;
    fy_sub  = H_SUB*1.0 # H_SUB*2.0;
    cx_sub  = W_SUB/2.0 + 0.5 # W_SUB/2.0 + 0.5;
    cy_sub  = H_SUB/2.0 + 0.5 # H_SUB/2.0 + 0.5;
    H_SUB_M = FOCAL_LENGTH*(H_SUB/fy_sub);
    W_SUB_M = FOCAL_LENGTH*(W_SUB/fx_sub);
    N_SUB = W_SUB*H_SUB
    KK_sub = np.array([[fx_sub, 0, cx_sub],
                       [0, fy_sub, cy_sub],
                       [0, 0, 1]])
    
    def __init__(self, robot, sensor):
        self.robot = robot
        self.sensor = sensor
        
    def directions(self, subsampled=True):
        """
        Returns rays that emanate from the sensor pose through the image plane
        (in 'world' frame)
        """
        if subsampled:
            h, w, n = Camera.H_SUB, Camera.W_SUB, Camera.N_SUB
            H, W = Camera.HEIGHT_M, Camera.WIDTH_M
        else:
            h, w, n = Camera.HEIGHT, Camera.WIDTH, Camera.N
            H, W = Camera.H_SUB_M, Camera.W_SUB_M
        
        height_offsets = np.linspace(-H/2.0, H/2.0, h)
        width_offsets = np.linspace(-W/2.0, W/2.0, w)
        
        height_grid, width_grid = np.meshgrid(height_offsets, width_offsets)
        
        height_grid_vec = height_grid.reshape((n,1))
        width_grid_vec = width_grid.reshape((n,1))
        z_vec = np.zeros((n,1))
        
        offsets = np.hstack((width_grid_vec, height_grid_vec, z_vec))
        
        points_cam = (Camera.MAX_RANGE/Camera.FOCAL_LENGTH)*(np.tile(np.array([0,0,Camera.FOCAL_LENGTH]), (n,1)) + offsets)
        
        sensor_world_transform = self.sensor.GetTransform()
        
        directions = np.zeros((n,3))
        
        sensor_world_pos = sensor_world_transform[:3,3]
        global handles
        for i in xrange(n):
            p_cam = tfx.pose(points_cam[i,:])
            p_world_pos = np.dot(sensor_world_transform, np.array(p_cam.matrix))[0:3,3]
 
            direction = np.array(p_world_pos) - np.array(sensor_world_pos)
            directions[i,:] = direction
            
            #handles += utils.plot_segment(self.robot.GetEnv(), sensor_world_pos, sensor_world_pos + direction)

        
        return directions
    
    def get_hits(self, subsampled=True):
        """ From the current camera pose, returns the ray collision points """
        if subsampled:
            h, w, n = Camera.H_SUB, Camera.W_SUB, Camera.N_SUB
        else:
            h, w, n = Camera.HEIGHT, Camera.WIDTH, Camera.N
        dirs = self.directions(subsampled)
        
        origin_world_pos = self.sensor.GetTransform()[:3,3]
        rays = np.hstack((np.tile(origin_world_pos, (dirs.shape[0],1)), dirs))
        
        env = self.robot.GetEnv()
        with env:
            is_hits, hits = env.CheckCollisionRays(rays, None)
        
        is_hits = is_hits.reshape((w,h)).T
        hits = np.swapaxes(hits.reshape((w,h,6)), 0, 1)
        dirs = np.swapaxes(dirs.reshape((w,h,3)), 0, 1)
        zpoints = np.zeros((h,w,3))
        
        global handles
        for i in xrange(h):
            for j in xrange(w):
                if is_hits[i,j]:
                    zpoints[i,j,:] = hits[i,j,:3]
                else:
                    zpoints[i,j,:] = origin_world_pos + dirs[i,j,:]
                #handles.append(utils.plot_point(env, zpoints[j,i,:]))
                #print('(i,j): ({0},{1})'.format(i,j))
                #raw_input()
                
        return is_hits, zpoints
    
    def get_zbuffer(self, subsampled=True):
        is_hits, hits = self.get_hits(subsampled)
        
        sensor_world_pos = self.sensor.GetTransform()[:3,3]
        
        rows, cols = hits.shape[:2]
        zbuffer = np.zeros((rows,cols))
        for i in xrange(rows):
            for j in xrange(cols):
                zbuffer[i,j] = np.linalg.norm(hits[i,j,:] - sensor_world_pos)
                
        return zbuffer
    
    def get_pixel_from_point(self, x, subsampled=True, is_round=True):
        #x_mat = np.array(tfx.point(x).as_pose().matrix)
        #x_mat[:3,:3] = np.zeros((3,3))
        
        x_mat = np.eye(4)
        x_mat[0:3,3] = tfx.point(x).array
        
        x_mat_tilde = np.dot(np.linalg.inv(self.sensor.GetTransform()), x_mat)
        
        #print x_mat
        #print ''
        #print x_mat_tilde
        
        if subsampled:
            y = np.dot(Camera.KK_sub, x_mat_tilde[0:3,3])
        else:
            y = np.dot(Camera.KK, x_mat_tilde[0:3,3])
        
        y_pixel = y[1]/y[2]
        x_pixel = y[0]/y[2]
        
        if is_round:
            y_pixel, x_pixel = np.round(y_pixel), np.round(x_pixel)
            
        return (y_pixel, x_pixel)
    
    def is_in_fov(self, point, zbuffer, subsampled=True):
        point = tfx.point(point)
        if subsampled:
            h, w = Camera.H_SUB, Camera.W_SUB
        else:
            h, w = Camera.HEIGHT, Camera.WIDTH
            
        y, x = self.get_pixel_from_point(point, subsampled)
        
        if (y < 0) or (y >= h) or (x < 0) or (x >= w):
            return False
        
        if zbuffer[y,x] < (tfx.point(self.sensor.GetTransform()) - point).norm:
            return False
        
        if np.dot(np.linalg.inv(self.sensor.GetTransform()),point.as_pose().matrix)[2,3] < 0:
            return False
        
        return True
        
        
    
class VoxelGrid:
    def __init__(self, camera, pos_center, x_height, y_height, z_height, resolution=512):
        self.camera = camera
        self.resolution = resolution
        self.corner = np.array(pos_center) - np.array([x_height/2., y_height/2., z_height/2.])
        self.dx = x_height / float(resolution)
        self.dy = y_height / float(resolution)
        self.dz = z_height / float(resolution)
        self.radius = min([self.dx, self.dy, self.dz])/10.
        
        self.TSDF = np.ones((resolution, resolution, resolution))
        self.object = None
        self.ODF = -np.inf*np.ones((resolution, resolution, resolution))
        
        self.handles = list()
    
    """
    User-called methods
    """
    
    def update_TSDF(self):
        """ Update TSDF by getting collisions from current camera pose """
        is_hits, hits = self.camera.get_hits()
        rows, cols = is_hits.shape
        for r in xrange(rows):
            for c in xrange(cols):
                if is_hits[r,c]:
                    i, j, k = self.voxel_from_point(hits[r,c,:])
                    self.TSDF[i,j,k] = 0
                    for neighbor, _ in self.get_voxel_neighbors_and_dists([i,j,k]):
                        l, m, n = neighbor
                        self.TSDF[l,m,n] = 0
                    
    def update_ODF(self, obj):
        """ 
        Update object distance-field.
        Call this after update_TSDF so Dijkstra's goes around the environment.
        """
        # TEMP
        #self.handles = list()
        #env = self.camera.sensor.GetEnv()
        self.object = obj
        self.ODF = np.inf*np.ones((self.resolution, self.resolution, self.resolution))
        
        pq = utils.PriorityQueue()
        
        for i in xrange(self.resolution):
            for j in xrange(self.resolution):
                for k in xrange(self.resolution):
                    if self.TSDF[i,j,k] != 0:
                        pq.insert((i,j,k), np.inf)
        
     
        obj_voxel = tuple(self.voxel_from_point(obj))
        
        pq.delete(obj_voxel)
        pq.insert(obj_voxel, 0)
        
        iter = 0
        while True:
            print('iter: {0}'.format(iter))
            iter += 1
            
            curr_dist, curr = pq.pop()
            if curr is None:
                break
            
            #self.handles += [utils.plot_point(env, self.get_voxel_center(curr), size=self.radius)]
            #time.sleep(.005)
            #raw_input()
            
            i, j, k = curr
            self.ODF[i,j,k] = curr_dist
            
            for n, n_dist in self.get_voxel_neighbors_and_dists(curr):
                dist = curr_dist + n_dist
                i, j, k = n
                if dist < self.ODF[i,j,k]:
                    self.ODF[i,j,k] = dist
                    pq.delete(n)
                    pq.insert(n, dist)
        
    def signed_distance_complete(self):
        """
        Finds signed-distance from the current field-of-view to the object
        (where object is implicitly given by the ODF).
        Call this after update_TSDF and update_ODF.
        """
        zbuffer = self.camera.get_zbuffer()
        obj_in_fov = self.camera.is_in_fov(self.object, zbuffer)
        
        min_dist, min_voxel = np.inf, None
        for i in xrange(self.resolution):
            for j in xrange(self.resolution):
                for k in xrange(self.resolution):
                    print('({0},{1},{2})'.format(i,j,k))
                    if self.TSDF[i,j,k] != 0:
                        voxel_center = self.get_voxel_center([i,j,k])
                        if obj_in_fov:
                            if not self.camera.is_in_fov(voxel_center, zbuffer) and self.ODF[i,j,k] < min_dist:
                                min_dist = self.ODF[i,j,k]
                                min_voxel = [i,j,k]
                        else:
                            if self.camera.is_in_fov(voxel_center, zbuffer) and self.ODF[i,j,k] < min_dist:
                                min_dist = self.ODF[i,j,k]
                                min_voxel = [i,j,k]
                                
                        if self.camera.is_in_fov(voxel_center, zbuffer):
                            self.handles += [utils.plot_point(self.camera.sensor.GetEnv(), voxel_center, size=self.radius)]
                            
        sd = min_dist if obj_in_fov else -min_dist
        
        # TEMP
        self.handles += [utils.plot_point(self.camera.sensor.GetEnv(), self.get_voxel_center(min_voxel), size=5*self.radius)]
        
        return sd
    
    
    """
    Support methods
    """
        
    def get_voxel_neighbors_and_dists(self, voxel):
        """ Gets neighbors of voxel that are not part of the environment (i.e. TSDF != 0) """
        voxel = np.array(voxel)
        offsets_and_dists = [([0, 0, -1], self.dz),
                           ([1, 0, -1], np.linalg.norm([self.dx,0,self.dz])),
                           ([-1, 0, -1], np.linalg.norm([self.dx,0,self.dz])),
                           ([0, 1, -1], np.linalg.norm([0,self.dy,self.dz])),
                           ([0, -1, -1], np.linalg.norm([0,self.dy,self.dz])),
                           ([1, 1, -1], np.linalg.norm(np.array([self.dx,self.dy,self.dz]))),
                           ([1, -1, -1], np.linalg.norm(np.array([self.dx,self.dy,self.dz]))),
                           ([-1, 1, -1], np.linalg.norm(np.array([self.dx,self.dy,self.dz]))),
                           ([-1, -1, -1], np.linalg.norm(np.array([self.dx,self.dy,self.dz]))),
                           
                           ([1, 0, 0], self.dx),
                           ([-1, 0, 0], self.dx),
                           ([0, 1, 0], self.dy),
                           ([0, -1, 0], self.dy),
                           ([1, 1, 0], np.linalg.norm(np.array([self.dx,self.dy,0]))),
                           ([1, -1, 0], np.linalg.norm(np.array([self.dx,self.dy,0]))),
                           ([-1, 1, 0], np.linalg.norm(np.array([self.dx,self.dy,0]))),
                           ([-1, -1, 0], np.linalg.norm(np.array([self.dx,self.dy,0]))),
                           
                           ([0, 0, 1], self.dz),
                           ([1, 0, 1], np.linalg.norm([self.dx,0,self.dz])),
                           ([-1, 0, 1], np.linalg.norm([self.dx,0,self.dz])),
                           ([0, 1, 1], np.linalg.norm([0,self.dy,self.dz])),
                           ([0, -1, 1], np.linalg.norm([0,self.dy,self.dz])),
                           ([1, 1, 1], np.linalg.norm(np.array([self.dx,self.dy,self.dz]))),
                           ([1, -1, 1], np.linalg.norm(np.array([self.dx,self.dy,self.dz]))),
                           ([-1, 1, 1], np.linalg.norm(np.array([self.dx,self.dy,self.dz]))),
                           ([-1, -1, 1], np.linalg.norm(np.array([self.dx,self.dy,self.dz])))]
                           
        
        neighbors = list()
        for offset, dist in offsets_and_dists:
            neighbor = voxel + np.array(offset)
            if np.min(neighbor) >= 0 and np.max(neighbor) < self.resolution:
                i, j, k = neighbor
                if self.TSDF[i,j,k] != 0:
                    neighbors.append((tuple(neighbor), dist))
                
        return neighbors
                
    def voxel_from_point(self, point):
        """ Get voxel (x,y,z) closest to point """
        rel_x, rel_y, rel_z = list(point - self.corner)
        voxel = [int(rel_x/self.dx), int(rel_y/self.dy), int(rel_z/self.dz)]
        
        for voxel_coord in voxel:
            if voxel_coord < 0 or voxel_coord >= self.resolution:
                print('point {0} not in the VoxelGrid!'.format(point))
                return None
            
        return voxel
    
    def get_voxel_center(self, voxel):
        i, j, k = list(voxel)
        center = self.corner + np.array([i*self.dx, j*self.dy, k*self.dz])
        return center
    
    """
    Plotting methods
    """
    
    def plot_TSDF(self):
        """ Only plot objects hit (i.e. TSDF[i,j,k] = 0) """
        self.handles = list()
        env = self.camera.sensor.GetEnv()
        for i in xrange(self.resolution):
            for j in xrange(self.resolution):
                for k in xrange(self.resolution):
                    if self.TSDF[i,j,k] == 0:
                        self.handles += [utils.plot_point(env, self.get_voxel_center([i,j,k]), size=self.radius)]
    
    def plot_ODF(self):
        """ Plots ODF """
        self.handles = list()
        env = self.camera.sensor.GetEnv()
        
        max_dist = -np.inf
        for i in xrange(self.resolution):
            for j in xrange(self.resolution):
                for k in xrange(self.resolution):
                    if self.ODF[i,j,k] > max_dist and self.ODF[i,j,k] != np.inf:
                        max_dist = self.ODF[i,j,k]
        
        for i in xrange(self.resolution):
            for j in xrange(self.resolution):
                for k in xrange(self.resolution):
                    dist = self.ODF[i,j,k]
                    if dist >= 0 and dist != np.inf:
                        dist_pct = dist/max_dist
                        rgb = colorsys.hsv_to_rgb((2/3.)*dist_pct,1,1)
                        self.handles += [utils.plot_point(env, self.get_voxel_center([i,j,k]), size=self.radius, color=rgb)]
                    else:
                        self.handles += [utils.plot_point(env, self.get_voxel_center([i,j,k]), size=self.radius, color=[0,0,0])]
    
    def plot_FOV(self):
        env = self.camera.sensor.GetEnv()
        color = [1,0,0]
        
        is_hits, hits = self.camera.get_hits()
        rows, cols = is_hits.shape
        for i in xrange(rows):
            for j in xrange(cols):
                if is_hits[i,j]:
                    voxel = self.voxel_from_point(hits[i,j,:])
                    voxel_center = self.get_voxel_center(voxel)
                    self.handles += [utils.plot_point(env, voxel_center, size=self.radius, color=color)]
                    
        cam_pos = self.camera.sensor.GetTransform()[:3,3]
        self.handles += utils.plot_segment(env, cam_pos, hits[0,0,:], color=color)
        self.handles += utils.plot_segment(env, cam_pos, hits[0,cols-1,:], color=color)
        self.handles += utils.plot_segment(env, cam_pos, hits[rows-1,0,:], color=color)
        self.handles += utils.plot_segment(env, cam_pos, hits[rows-1,cols-1,:], color=color)
    
    def plot_centers(self):
        self.handles = list()
        env = self.camera.sensor.GetEnv()
        for i in xrange(self.resolution):
            for j in xrange(self.resolution):
                for k in xrange(self.resolution):
                    self.handles += [utils.plot_point(env, self.get_voxel_center([i,j,k]), size=self.radius)]
        
        
def test_camera():
    env = rave.Environment()
    env.Load('../envs/pr2-empty.env.xml')
    env.SetViewer('qtcoin')
    #env.GetViewer().SendCommand('SetFiguresInCamera 1') # also shows the figures in the image
    time.sleep(1)
    robot = env.GetRobots()[0]
    
    sensor = robot.GetAttachedSensor('head_cam').GetSensor()

    cam = Camera(robot, sensor)
    
    is_hits, hits = cam.get_hits(subsampled=True)
    zbuffer = cam.get_zbuffer()
    
    height, width, _ = hits.shape
    global handles
    for i in xrange(height):
        for j in xrange(width):
            #print zpoints[i,j,:]
            if is_hits[i,j]:
                if cam.is_in_fov(hits[i,j,:], zbuffer):
                    handles += [utils.plot_point(env, hits[i,j,:], size=.01)]
    
    IPython.embed()
    
def test_voxel_grid():
    env = rave.Environment()
    env.Load('../envs/pr2-empty.env.xml')
    env.SetViewer('qtcoin')
    #env.GetViewer().SendCommand('SetFiguresInCamera 1') # also shows the figures in the image
    time.sleep(1)
    robot = env.GetRobots()[0]
    
    sensor = robot.GetAttachedSensor('head_cam').GetSensor()

    cam = Camera(robot, sensor)

    table = env.GetKinBody('table')
    table_pos = tfx.point(table.GetTransform()).array
    
    vgrid = VoxelGrid(cam, table_pos + np.array([0,0,0]), 1.5, 2, 1, resolution = 50)
    #vgrid.plot_centers()
    
    object = table_pos + np.array([.15,0,-.2])
    vgrid.object = object # TEMP
    h = utils.plot_point(env, object, size=.05, color=[1,0,0])
    
    vgrid.update_TSDF()
    #vgrid.update_ODF(object)
    #vgrid.signed_distance_complete()
    
    #vgrid.plot_TSDF()
    #vgrid.plot_ODF()
    #vgrid.plot_FOV()
    
    IPython.embed()
    
if __name__ == '__main__':
    #test_camera()
    test_voxel_grid()