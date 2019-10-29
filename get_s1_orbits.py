# Python script to get orbits from a given date
# Input: something of the format get_s1_orbits(20171202,s1a)
# Copies the correct orbit into the directory that called it. 

from subprocess import call
import sys
import sentinel_utilities

date=sys.argv[1]
sat =sys.argv[2]
eof_dir="/Volumes/Ironwolf/S1_orbits"
# Should find a more automatic way to get this (like from batch.config)

print("Copying orbit file into current directory...")
print("date is: %s " % date);
print("satellite is: %s" % sat);

eof_name = sentinel_utilities.get_eof_from_date_sat(date, sat, eof_dir);
orbit_name='ORBIT_'+sat.upper()+date+'.EOF'
call(['cp '+ eof_name+' .'],shell=True);
print("Successfully copied "+eof_name);
