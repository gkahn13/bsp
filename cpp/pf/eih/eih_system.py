import IPython

from pr2 import pr2_sim

class EihSystem:
    def __init__(self):
        self.brett = pr2_sim.PR2('envs/pr2-test.env.xml')
        
        self.DT = 1.0
        
    def dynfunc(self, x, u):
        pass

def test_eih_system():
    brett = pr2_sim.PR2('envs/pr2-test.env.xml')
    IPython.embed()

if __name__ == '__main__':
    test_eih_system()