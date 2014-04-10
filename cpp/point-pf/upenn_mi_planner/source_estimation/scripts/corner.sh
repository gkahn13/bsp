#!/bin/bash

# Home, James!
rosrun rvo_move move_client __name:=bar _x:=-7.2 _y:=4.4 move_server:=/scarab26/move_server &
rosrun rvo_move move_client __name:=quux _x:=-7.5 _y:=5.0 move_server:=/scarab29/move_server &
rosrun rvo_move move_client __name:=foo _x:=-8.0 _y:=4.8 move_server:=/scarab32/move_server &

