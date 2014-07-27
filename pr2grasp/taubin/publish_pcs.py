#!/usr/bin/env python

import rospy, roslib
roslib.load_manifest('sensor_msgs')
import sensor_msgs.msg as sm

import argparse

class PublishPCs:
    def __init__(self, subscribe_topic, publish_topic):
        rospy.loginfo('Listening for point-cloud topic: {0}'.format(subscribe_topic))
        self.sub = rospy.Subscriber(subscribe_topic, sm.PointCloud2, self._pc_callback)
        
        rospy.loginfo('Publishing to point-cloud topic: {0}'.format(publish_topic))
        self.pub = rospy.Publisher(publish_topic, sm.PointCloud2)
        
        self.last_pc = None
        while self.last_pc is None and not rospy.is_shutdown():
            rospy.sleep(.05)
    
    def run(self):
        num_published = 0
        while not rospy.is_shutdown():
            rospy.loginfo('Press enter to publish')
            raw_input()
            
            self.pub.publish(self.last_pc)
            
            rospy.loginfo('Number published: {0}'.format(num_published))
            num_published += 1
    
    def _pc_callback(self, msg):
        self.last_pc = msg
    
if __name__ == '__main__':
    rospy.init_node('publish_pcs', anonymous=True)
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--subscribe',type=str,default='/camera/depth_registered/points')
    parser.add_argument('--publish',type=str,default='/save_pc')
    args = parser.parse_args()
    
    publish_pcs = PublishPCs(args.subscribe, args.publish)
    publish_pcs.run()