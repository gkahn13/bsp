#!/bin/bash

# Home, James!
MS=move_server:=/scarab32/move_server
rosrun rvo_move move_client _x:=-7.2 _y:=5.0 $MS
rosrun rvo_move move_client _x:=1.2 _y:=4.7 $MS
rosrun rvo_move move_client _x:=1.2 _y:=28.2 $MS
rosrun rvo_move move_client _x:=-7.4 _y:=28.2 $MS
rosrun rvo_move move_client _x:=-7.2 _y:=5.0 $MS
