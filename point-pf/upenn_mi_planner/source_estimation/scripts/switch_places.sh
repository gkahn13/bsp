#!/bin/bash

S1=scarab26
S2=scarab32

S1_X=$(rostopic echo -n 1 /$S1/amcl_pose/pose/pose/position/x | head -1)
S1_Y=$(rostopic echo -n 1 /$S1/amcl_pose/pose/pose/position/y | head -1)
S2_X=$(rostopic echo -n 1 /$S2/amcl_pose/pose/pose/position/x | head -1)
S2_Y=$(rostopic echo -n 1 /$S2/amcl_pose/pose/pose/position/y | head -1)

rosrun rvo_move move_client _x:=$S1_X _y:=$S1_Y move_server:=/$S2/move_server __name:=foo &
rosrun rvo_move move_client _x:=$S2_X _y:=$S2_Y move_server:=/$S1/move_server __name:=bar &
wait
