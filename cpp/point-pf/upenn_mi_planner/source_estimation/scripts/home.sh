#!/bin/bash

# Home, James!
rosrun rvo_move move_client __name:=bar _x:=-9.4 _y:=1.0 move_server:=/scarab26/move_server &
rosrun rvo_move move_client __name:=quux _x:=-9.4 _y:=0.5 move_server:=/scarab29/move_server &
rosrun rvo_move move_client __name:=foo _x:=-9.0 _y:=1.5 move_server:=/scarab32/move_server &

