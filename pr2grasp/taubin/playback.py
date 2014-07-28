import rospy, roslib, rosbag
import tf
roslib.load_manifest('sensor_msgs')
import sensor_msgs.msg as sm

import argparse
import threading

import IPython

# don't forget
# roslaunch localization_sensor.launch

class Playback:
    def __init__(self, bag_name, pc_publish_topic='/camera/depth_registered/points'):
        self.pcs = list()
        self.tfs = list()
        
        self.bag = rosbag.Bag(bag_name + '.bag')
        for topic, msg, t in self.bag.read_messages():
            if topic == 'tf':
                self.tfs.append(msg)
            if topic == 'save_pc':
                self.pcs.append(msg)
                
        self.playback_bag = rosbag.Bag(bag_name + '_playback.bag', 'w')
        
        rospy.loginfo('Number of point-clouds read: {0}'.format(len(self.pcs)))
        rospy.loginfo('Number of tfs saved: {0}'.format(len(self.tfs)))
        
        self.pub = rospy.Publisher(pc_publish_topic, sm.PointCloud2)
        self.br = tf.TransformBroadcaster()
        
        self.tfs = dict()
        self.tfs_thread = threading.Thread(target=self._tfs_loop)
        self.tfs_thread.start()
        
    def _tfs_loop(self):
        while not rospy.is_shutdown():
            tfs = dict(self.tfs)
            for t_list in tfs.values():
                t = t_list[-1]
                trans = [t.transform.translation.x, t.transform.translation.y, t.transform.translation.z]
                rot = [t.transform.rotation.x, t.transform.rotation.y, t.transform.rotation.z, t.transform.rotation.w]
                self.br.sendTransform(trans,
                                      rot,
                                      rospy.Time.now(),
                                      t.child_frame_id,
                                      t.header.frame_id)
                
            rospy.sleep(.001)
            
    def play(self):
        pass
        
    def create_playback_bag(self):
        save_pc_indices = [i for i, (topic, msg, t) in enumerate(self.bag.read_messages()) if topic == 'save_pc']
        print save_pc_indices
        
        bag_msgs = list()
        for (topic, msg, t) in self.bag.read_messages():
            bag_msgs.append((topic, msg, t))
        new_bag_msgs = [None for _ in xrange(save_pc_indices[-1]+1)]
        
        prev_pc_index = 0
        new_bag_msgs[0] = bag_msgs[0]
        for curr_pc_index in save_pc_indices:
            new_bag_msgs[curr_pc_index] = bag_msgs[curr_pc_index]
            #print bag_msgs[prev_pc_index][-1]
            
            tfs = dict()
            for index in xrange(curr_pc_index-1, prev_pc_index, -1):
                topic, msg, t = bag_msgs[index]
                if topic == 'tf':
                    new_bag_msgs[index] = (topic, msg, t)
                    
                    """
                    for transform in msg.transforms:
                        key = (transform.child_frame_id, transform.header.frame_id)
                        if key in tfs.keys():
                            #future_msg = tf.msg.tfMessage(tfs[key])
                            #future_msg.transforms.transforms[0].header.stamp = t
                            new_bag_msgs[index] = (topic, tfs[key], t)
                        else:
                            new_bag_msgs[index] = (topic, msg, t)
                            tfs[key] = msg
                    """
            
            prev_pc_index = curr_pc_index
        
        for i, item in enumerate(new_bag_msgs):
            if item is None:
                print i        
        new_bag_msgs = [item for item in new_bag_msgs if item is not None]
        
        for (topic, msg, t) in new_bag_msgs:
            self.playback_bag.write(topic, msg, t=t)
            
        self.playback_bag.close()
    
    def run_all(self):
        num_pcs_published = 0
        tfs = dict()
        
        for topic, msg, t in self.bag.read_messages():
            
            if topic == 'tf':
                for t in msg.transforms:
                    key = (t.child_frame_id, t.header.frame_id)
                    if key in tfs.keys():
                        self.tfs[key].append(t)
                    else:
                        self.tfs[key] = [t]
                    """
                    # translation, rotation, time, child, parent
                    trans = [t.transform.translation.x, t.transform.translation.y, t.transform.translation.z]
                    rot = [t.transform.rotation.x, t.transform.rotation.y, t.transform.rotation.z, t.transform.rotation.w]
                    self.br.sendTransform(trans,
                                          rot,
                                          rospy.Time.now(),
                                          t.child_frame_id,
                                          t.header.frame_id)
                    """
            if topic == 'save_pc':
                msg.header.stamp = rospy.Time.now()
                self.pub.publish(msg)
                
                """
                rospy.loginfo('publishing loop')
                while not rospy.is_shutdown(): # put this loop on a separate thread
                    
                    #self.br.sendTransform([-0.106925711716, -0.00652973239027, -0.0566985277547],
                    #              [0.501389436235, 0.466364575859, 0.509987956093, 0.520600614921],
                    #              rospy.Time.now(),
                    #              'camera_rgb_optical_frame',
                    #              'r_gripper_tool_frame')
                    
                    for t_list in tfs.values():
                        t = t_list[-1]
                        trans = [t.transform.translation.x, t.transform.translation.y, t.transform.translation.z]
                        rot = [t.transform.rotation.x, t.transform.rotation.y, t.transform.rotation.z, t.transform.rotation.w]
                        self.br.sendTransform(trans,
                                              rot,
                                              rospy.Time.now(),
                                              t.child_frame_id,
                                              t.header.frame_id)
                        
                    rospy.sleep(.001)
                """
                   
                    
                rospy.loginfo('Published pc {0}'.format(num_pcs_published))
                num_pcs_published += 1
                
                print('Press enter')
                raw_input()
            
            if rospy.is_shutdown():
                break
            
            rospy.sleep(.001)
        #for i, pc in enumerate(self.pcs):
        #    rospy.loginfo('Press enter to publish pc {0}'.format(i))
        #    raw_input()
        #    
        #    self.pub.publish(pc)
    
if __name__ == '__main__':
    rospy.init_node('playback', anonymous=True)
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--bag',type=str,default='bags/spray_can')
    args = parser.parse_args()
    
    playback = Playback(args.bag)
    #playback.create_playback_bag()
    
    IPython.embed()