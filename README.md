Belief Space Planning
=====================

Various examples and implementations of belief space planning
====

point: Ideas first tested in canonical light-dark example. **NOT MAINTAINED**

===

Experiments for [WAFR 2014 submission](http://goldberg.berkeley.edu/pubs/Patil-WAFR2014-CFGBSP.pdf)
arm: For a 6-DOF arm and a fixed camera, finds path from start to end that minimizes uncertainty
parameter: For a dynamical system with unknown parameters, finds control inputs that determines the unknown parameters efficiently 
slam: For a car and landmarks, finds a path that minimizes uncertainty about both the car and landmark positions while reaching specified waypoints

===

**ONGOING:** exploring non-parametric methods, namely particle filters

point-pf: Uncertainty about robot represented by particles, optimize directly over them
pf: Uncertainty about something in the environment
  explore: agents seek a target object located randomly in an environment
  boxes: agent goal is to localize the position of box(es)
  eih: eye-in-hand, localize an object with a kinect like sensor on an end effector

===

The code is not suitable for public consumption yet in its current form. If you are interested in using it, please send Greg (gkahn [at] berkeley.edu) an email saying that you would like to use the code and what you plan to do with it, and we will try to help you out. 
