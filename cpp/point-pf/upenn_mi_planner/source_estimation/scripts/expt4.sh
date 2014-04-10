#!/bin/bash

# Target
rosrun rvo_move move_client __name:=bar  _x:=51   _y:=51 move_server:=/scarab29/move_server &
# Seekers
rosrun rvo_move move_client __name:=quux _x:=70.8 _y:=55.8   move_server:=/scarab32/move_server &
rosrun rvo_move move_client __name:=foo  _x:=72.0  _y:=55.8   move_server:=/scarab26/move_server &

wait