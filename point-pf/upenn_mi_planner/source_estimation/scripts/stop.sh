#!/bin/bash
rosrun rvo_move move_client __name:=bar _x:=-50.5 _y:=4.8 move_server:=/scarab26/move_server &
rosrun rvo_move move_client __name:=quux _x:=-50.5 _y:=4.0 move_server:=/scarab29/move_server &
rosrun rvo_move move_client __name:=foo _x:=-50.0 _y:=4.5 move_server:=/scarab32/move_server &
