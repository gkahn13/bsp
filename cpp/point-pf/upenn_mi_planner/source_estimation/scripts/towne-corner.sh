#!/bin/bash

# Target
rosrun rvo_move move_client __name:=bar  _x:=54.0   _y:=51.0 move_server:=/scarab29/move_server &
# Seekers
rosrun rvo_move move_client __name:=quux _x:=54.63 _y:=51.8   move_server:=/scarab32/move_server &
rosrun rvo_move move_client __name:=foo  _x:=55.3  _y:=51.0   move_server:=/scarab29/move_server &

