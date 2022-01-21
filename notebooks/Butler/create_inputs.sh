#!/bin/sh
for var in 13141; do
#for var in 13038 13039 13206 13207 13208 13209 13210 13211 13212 13213 13214 13215 13216 13217 13218 13219 13220 13221 13222 13223 13224 13148 13149 13150 13151; do
    echo "Run $var"
    find -L /sps/lsst/groups/FocalPlane/SLAC/run5/$var -name "*.fits" > $var"_list.txt"
done 
