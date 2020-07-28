

import sys

if len(sys.argv) != 3:
    print("Please provide the name of the search file and bounds file");
    sys.exit(1);
else:
    infile=open(sys.argv[1]);
outfile=open(sys.argv[2],'w');
configline=infile.readline();
temp = configline.split(); 

point_index=temp.index('-p') if '-p' in temp else -1;
region_index=temp.index('-r') if '-r' in temp else -1;

if point_index != -1:
    print("We are drawing a point.")
    point_string=temp[point_index+1];
    temp2=point_string.split('/');
    outfile.write('%s %s\n' % ( temp2[0], temp2[1] ) );
    
elif region_index != -1:
    print("We are drawing a region.")
    region_string=temp[region_index+1];
    temp2=region_string.split('/');
    outfile.write('%s %s\n' % ( temp2[0], temp2[2] ) );
    outfile.write('%s %s\n' % ( temp2[0], temp2[3] ) );
    outfile.write('%s %s\n' % ( temp2[1], temp2[3] ) );
    outfile.write('%s %s\n' % ( temp2[1], temp2[2] ) );
    outfile.write('%s %s\n' % ( temp2[0], temp2[2] ) );
else:
    print("could not find region or point.")

infile.close();
outfile.close();
