import numpy as np

epsilon = 1e-5

class Triangle:
    def __init__(self, a, b, c):
        self.a, self.b, self.c = np.array(a), np.array(b), np.array(c)
        
    def closest_point_to(self, x):
        if self.is_inside(x):
            return None
        
        min_pt, min_dist = None, np.inf
        for s in self.segments():
            s_min_pt = s.closest_point_to(x)
            if np.linalg.norm(x - s_min_pt) < min_dist:
                min_dist = np.linalg.norm(x - s_min_pt)
                min_pt = s_min_pt
                
        return min_pt
    
    def is_inside(self, x):
        total_area = self.area()
        area0 = Triangle(self.a, self.b, x).area()
        area1 = Triangle(self.b, self.c, x).area()
        area2 = Triangle(self.c, self.a, x).area()
        
        is_correct_area = np.abs(total_area - (area0 + area1 + area2)) < epsilon
        
        return is_correct_area
    
    def area(self):
        a, b, c = self.a, self.b, self.c
        return np.abs((c[0]*(a[1] - b[1]) + a[0]*(b[1] - c[1]) + b[0]*(c[1] - a[1])) / 2.0)
        
    def segments(self):
        return (Segment(self.a, self.b), Segment(self.b, self.c), Segment(self.c, self.a))
        
class Segment:
    def __init__(self, p0, p1):
        self.p0, self.p1 = np.array(p0), np.array(p1)
        
    def closest_point_to(self, x):
        # min_{0<=t<=1} ||t*(p1-p0) + p0 - x||_{2}^{2}
        v = self.p1 - self.p0
        b = self.p0 - x
        
        t = -np.dot(v, b) / np.dot(v, v)
        if (0 <= t <= 1):
            intersection = t*(self.p1 - self.p0) + self.p0
            return intersection
        else:
            if np.linalg.norm(x - self.p0) < np.linalg.norm(x - self.p1):
                return self.p0
            else:
                return self.p1
        