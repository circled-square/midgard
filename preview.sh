#!/bin/bash

while true 
do
    target/release/examples/preview & sleep 3 && xdotool key Escape&
    sleep 2.8
done
