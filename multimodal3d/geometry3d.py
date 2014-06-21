import time
import numpy as np

import roslib
roslib.load_manifest('tfx')
import tfx
import tf.transformations as tft

import geometry2d
import utils

import IPython

epsilon = 1e-5
handles = list()

class Beam:
    def __init__(self, base, a, b, c, d):
        """
        A pyramid with orign base and points a,b,c,d arranged as
        
        b --- a
        |     |
        |     |
        c --- d
        """
        self.base = np.array(base)
        self.a = np.array(a)
        self.b = np.array(b)
        self.c = np.array(c)
        self.d = np.array(d)
        
    def is_inside(self, p):
        p = np.array(p)
        
        halfplanes = self.get_halfplanes()
        return np.min([h.contains(p) for h in halfplanes])    
    
    def get_halfplanes(self):
        base, a, b, c, d = self.base, self.a, self.b, self.c, self.d

        origins = [(base+a+d)/3.0,
                   (base+b+a)/3.0,
                   (base+c+b)/3.0,
                   (base+d+c)/3.0,
                   (a+b+c+d)/4.0]
    
        normals = [-np.cross(a-base, d-base),
                   -np.cross(b-base, a-base),
                   -np.cross(c-base, b-base),
                   -np.cross(d-base, c-base),
                   -np.cross(b-a, d-a)]
        normals = [n/np.linalg.norm(n) for n in normals]
        
        return [Halfplane(origin, normal) for origin, normal in zip(origins, normals)]
        
    def get_side(self, side):
        """
        @param side - 'right', 'top', 'left', 'bottom', 'front'
        
        @return list of triangles
        """
        if side == 'right':
            return [Triangle(self.base, self.a, self.d)]
        elif side == 'top':
            return [Triangle(self.base, self.a, self.b)]
        elif side == 'left':
            return [Triangle(self.base, self.b, self.c)]
        elif side == 'bottom':
            return [Triangle(self.base, self.c, self.d)]
        elif side == 'front':
            return [Triangle(self.a, self.b, self.c), Triangle(self.a, self.d, self.c)]
        else:
            return None
    
    def plot(self, env, with_sides=True, color=(1,0,0)):
        global handles
        
        base, a, b, c, d = self.base, self.a, self.b, self.c, self.d
        
        if with_sides:
            handles += utils.plot_segment(env, base, a)
            handles += utils.plot_segment(env, base, b)
            handles += utils.plot_segment(env, base, c)
            handles += utils.plot_segment(env, base, d)
        
        handles += utils.plot_segment(env, a, b)
        handles += utils.plot_segment(env, b, c)
        handles += utils.plot_segment(env, c, d)
        handles += utils.plot_segment(env, d, a)
        
class Halfplane:
    def __init__(self, origin, normal):
        self.origin = origin
        self.normal = normal
        
    def contains(self, x):
        return np.dot(self.normal, x - self.origin) >= epsilon
    
    def plot(self, env, color=(0,0,1)):
        global handles
        
        o, n = self.origin, self.normal
        handles += utils.plot_segment(env, o, o + .05*n, color=color)

class Triangle:
    def __init__(self, a, b, c):
        self.a, self.b, self.c = np.array(a), np.array(b), np.array(c)
        
    def align_with(self, target):
        """
        Aligns the normal of this triangle to target
        
        @param target - 3 dim np.array
        @return (rotated triangle, rotation matrix)
        """
        target = np.array(target)
        source = np.cross(self.b - self.a, self.c - self.a)
        source /= np.linalg.norm(source)
    
        rotation = np.eye(3)
        
        dot = np.dot(source, target)
        if not np.isnan(dot):
            angle = np.arccos(dot)
            if not np.isnan(angle):
                cross = np.cross(source, target)
                cross_norm = np.linalg.norm(cross)
                if not np.isnan(cross_norm) and not cross_norm < epsilon:
                    cross = cross / cross_norm
                    rotation = tft.rotation_matrix(angle, cross)[:3,:3]

        return (Triangle(np.dot(rotation, self.a),
                        np.dot(rotation, self.b),
                        np.dot(rotation, self.c)),
                rotation)
        
    def closest_point_to(self, p):
        """
        Find distance to point p
        by rotating and projecting
        then return that closest point unrotated
        """
        p = np.array(p)
        # align with z-axis so all triangle have same z-coord
        tri_rot, rot = self.align_with([0,0,1])
        tri_rot_z = tri_rot.a[-1]
        p_rot = np.dot(rot, p)
        
        p_2d = p_rot[:2]
        tri_2d = geometry2d.Triangle(tri_rot.a[:2], tri_rot.b[:2], tri_rot.c[:2])
        
        if tri_2d.is_inside(p_2d):
            # projects onto triangle, so return difference in z
            return np.dot(np.linalg.inv(rot), np.array(list(p_2d) + [tri_rot_z]))
        else:
            closest_pt_2d = tri_2d.closest_point_to(p_2d)
            
            closest_pt_3d = np.array(list(closest_pt_2d) + [tri_rot_z])
            
            return np.dot(np.linalg.inv(rot), closest_pt_3d)
        
    def distance_to(self, p):
        """
        Find distance to point p
        by rotating and projecting
        """
        closest_pt = self.closest_point_to(p)
        return np.linalg.norm(p - closest_pt)
    
    def area(self):
        tri_rot, rot = self.align_with([0,0,1])
        tri_2d = geometry2d.Triangle(tri_rot.a[:2], tri_rot.b[:2], tri_rot.c[:2])
        return tri_2d.area()
        
    def plot(self, env, color=(1,0,0)):
        global handles
        
        handles += utils.plot_segment(env, self.a, self.b, color)
        handles += utils.plot_segment(env, self.b, self.c, color)
        handles += utils.plot_segment(env, self.c, self.a, color)
        
        
        
def test_align_with():
    t = Triangle([0,0,1.2], [0,1,1.2], [1,0,1.2])
    
    t_rot, rot = t.align_with([0,0,1])
    
    print('t_rot:\n{0}\n{1}\n{2}'.format(t_rot.a, t_rot.b, t_rot.c))
        
def test_distance_to():
    t = Triangle([0,0,0], [0,1,0], [1,0,0])
    
    p = [0, 0, 1]
    dist = t.distance_to(p)
    print('Distance should be 1')
    print('Computed distance is: {0}'.format(dist))
    
    p = [-1, 0, 1]
    dist = t.distance_to(p)
    print('Distance should be sqrt(2)')
    print('Computed distance is: {0}'.format(dist))
        
def test_distance_to_plot():
    import openravepy as rave
    
    env = rave.Environment()
    env.SetViewer('qtcoin')
    
    t = Triangle(np.random.rand(3), np.random.rand(3), np.random.rand(3))
    t.plot(env)
    
    p = 2*np.random.rand(3)
    closest_pt = t.closest_point_to(p)
    
    handles.append(utils.plot_point(env, p, color=(0,0,1)))
    handles.append(utils.plot_point(env, closest_pt))
    
    IPython.embed()
    
def test_beam_inside():
    import openravepy as rave
    
    env = rave.Environment()
    env.SetViewer('qtcoin')
    time.sleep(1)
    
    base = [0,0,0]
    a = [.1, .1, .5]
    b = [-.1, .1, .5]
    c = [-.1, -.1, .5]
    d = [.1, -.1, .5]
    
    beam = Beam(base, a, b, c, d)
    beam.plot(env)
    
    halfplanes = beam.get_halfplanes()
    for h in halfplanes:
        h.plot(env)
    
    """
    p0 = [0, 0, .3]
    p1 = [0, 0, .55]
    p2 = [.2, 0, .3]
    
    for p in [p0, p1, p2]:
        print('is_inside: {0}'.format(beam.is_inside(p)))
        handles.append(utils.plot_point(env, p))
    """
        
    for i in xrange(1000):
        p = np.random.rand(3)
        is_inside = beam.is_inside(p)
        print('is_inside: {0}'.format(is_inside))
        h = utils.plot_point(env, p)
        if is_inside:
            raw_input()
        
    
    IPython.embed()
        
if __name__ == '__main__':
    #test_align_with()
    #test_distance_to()
    #test_distance_to_plot()
    test_beam_inside()