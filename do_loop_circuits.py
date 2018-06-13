

# Step 1: Glob the intf_all
# Step 2: Form as many loops as possible
# Step 3: Make images of all the possible loops. 

import numpy as np 
import matplotlib.pyplot as plt 
import glob
import subprocess
import netcdf_read_write

def identify_all_loops():
	loops=[]; 
	edges=[];
	nodes=[];
	directories = glob.glob("intf_all/*_*");
	for item in directories:
		edgepair=item.split('/')[-1];
		edges.append(edgepair);
		nodes.append(edgepair.split('_')[0]);
		nodes.append(edgepair.split('_')[1]);
	nodes=list(set(nodes));

	# Making a graph of nodes and edges
	startnode=[];
	endnode=[];
	for i in range(len(nodes)):
		startnode.append(nodes[i]);
		thisend=[];
		for j in range(len(nodes)):
			if nodes[i]+'_'+nodes[j] in edges or nodes[j]+'_'+nodes[i] in edges:
				thisend.append(nodes[j])
		endnode.append(thisend);

	graph=dict([ (startnode[i], endnode[i]) for i in range(len(startnode)) ]);

	# Now we find loops. 
	# Start in chronological order. Step through the main loop.  
	# Within, have a loop of connectors, and another loop of connectors. If the innermost 
	# layer comes back to the starter, then you have a loop of 3. 
	for node0 in startnode:
		for node1 in graph[node0]:
			possible_list=graph[node1];
			for node2 in possible_list:
				if node2==node0:
					continue;
				elif node0 in graph[node2]:
					loop_try = [node0, node1, node2];
					loop_try = sorted(loop_try);
					if loop_try not in loops:
						loops.append(loop_try);

	print("Number of Images (nodes): % s" % (len(startnode)));
	print("Number of Loops (circuits): % s" % (len(loops)));
	return loops;



def towards_zero_2n_pi(phase_value):
	phase_value = phase_value + np.pi;
	residual = phase_value / (2*np.pi);
	if residual>0:
		phase_jump = np.floor(residual);
	else:
		phase_jump = np.floor(residual);
	return phase_jump;
	# Right now loop 24 is giving trouble... returning NAN for some reason. 


def compute_loops(all_loops, loops_dir, loops_guide, reference_pixel):
	subprocess.call(['mkdir','-p',loops_dir],shell=False);
	ofile=open(loops_dir+loops_guide,'w');
	for i in range(len(all_loops)):
		ofile.write("Loop %d: %s %s %s\n" % (i,all_loops[i][0], all_loops[i][1], all_loops[i][2]) );
	ofile.close();

	unwrapped='unwrap.grd'
	wrapped='phasefilt.grd'

	for i in range(25, len(all_loops)):
		edge1=all_loops[i][0]+'_'+all_loops[i][1];
		edge2=all_loops[i][1]+'_'+all_loops[i][2];
		edge3=all_loops[i][0]+'_'+all_loops[i][2];
		[xdata, ydata, z1 ] = netcdf_read_write.read_grd_xyz('intf_all/'+edge1+'/'+unwrapped);
		z2 = netcdf_read_write.read_grd('intf_all/'+edge2+'/'+unwrapped);
		z3 = netcdf_read_write.read_grd('intf_all/'+edge3+'/'+unwrapped);

		[xdata, ydata, wr_z1 ] = netcdf_read_write.read_grd_xyz('intf_all/'+edge1+'/'+wrapped);
		wr_z2 = netcdf_read_write.read_grd('intf_all/'+edge2+'/'+wrapped);
		wr_z3 = netcdf_read_write.read_grd('intf_all/'+edge3+'/'+wrapped);

		print("Loop "+str(i)+":");

		histdata_raw=[];
		histdata_fix=[];
		znew_raw = np.zeros(np.shape(z1));
		znew_fix = np.zeros(np.shape(z1));
		errorcount=0;
		for j in range(np.shape(z1)[0]):
			for k in range(np.shape(z1)[1]):
				
				if len(reference_pixel)>0:
					wr1=wr_z1[j][k]-wr_z1[reference_pixel[0], reference_pixel[1]];
					z1_adj=z1[j][k]-z1[reference_pixel[0], reference_pixel[1]];
					wr2=wr_z2[j][k]-wr_z2[reference_pixel[0], reference_pixel[1]];
					z2_adj=z2[j][k]-z2[reference_pixel[0], reference_pixel[1]];
					wr3=wr_z3[j][k]-wr_z3[reference_pixel[0], reference_pixel[1]];
					z3_adj=z3[j][k]-z3[reference_pixel[0], reference_pixel[1]];
				else:
					wr1=wr_z1[j][k];
					wr2=wr_z2[j][k];
					wr3=wr_z3[j][k];
					z1_adj=z1[j][k];
					z2_adj=z2[j][k];
					z3_adj=z3[j][k];

				# Using equation from Heresh Fattahi's PhD thesis to isolate unwrapping errors. 
				wrapped_closure_raw=wr_z1[j][k]+wr_z2[j][k]-wr_z3[j][k];
				wrapped_closure_fix=wr1+wr2-wr3;
				offset_before_unwrapping=np.mod(wrapped_closure_fix,2*np.pi);
				if offset_before_unwrapping>np.pi:
					offset_before_unwrapping=offset_before_unwrapping-2*np.pi;  # send it to the -pi to pi realm. 
				
				unwrapped_closure_raw=z1[j][k]+z2[j][k]-z3[j][k];
				unwrapped_closure_fix=z1_adj+z2_adj-z3_adj;

				znew_raw[j][k]=unwrapped_closure_raw - offset_before_unwrapping;
				znew_fix[j][k]=unwrapped_closure_fix - offset_before_unwrapping;

				if ~np.isnan(znew_raw[j][k]):
					histdata_raw.append(znew_raw[j][k]/np.pi);
				if ~np.isnan(znew_fix[j][k]):
					histdata_fix.append(znew_fix[j][k]/np.pi);
				if abs(znew_fix[j][k])>0.5:
					errorcount=errorcount+1;

		errorpixels = round(100*float(errorcount)/len(histdata_fix),2);
		print("Most common raw loop sum: ")
		print(np.median(histdata_raw));
		print("Most common fix loop sum: ")
		print(np.median(histdata_fix));
		print("\n");

		make_plot(xdata, ydata, znew_fix, loops_dir+'phase_closure_'+str(i)+'.eps', errorpixels);
		make_histogram(histdata_fix, loops_dir+'histogram_'+str(i)+'.eps');

	return;



def make_plot(x, y, z, plotname, errorpixels):
	fig = plt.figure();
	plt.imshow(z, cmap='jet',aspect=0.7);
	plt.gca().invert_yaxis()
	plt.gca().invert_xaxis()
	plt.gca().get_xaxis().set_ticks([]);
	plt.gca().get_yaxis().set_ticks([]);
	plt.title('Phase Closure: %.2f %% error pixels' % (errorpixels) );
	plt.gca().set_xlabel("Range",fontsize=16);
	plt.gca().set_ylabel("Azimuth",fontsize=16);
	cb = plt.colorbar();
	cb.set_label('phase residual', size=16);
	plt.savefig(plotname);
	plt.close();
	return;

def make_histogram(histdata, plotname):
	plt.figure();
	plt.hist(histdata,bins=700);
	plt.yscale('log');
	plt.xlabel('Phase Difference (#*pi radians)');
	plt.ylabel('number of pixels');
	plt.savefig(plotname);
	plt.close();
	return;


if __name__=="__main__":
	loops_dir="Phase_Circuits/"
	loops_guide="loops.txt";
	reference_pixel = [177, 253]; # doing correction for phase ambiguity
	# reference_pixel = [];  # not doing any correction for phase ambiguity


	all_loops=identify_all_loops();
	compute_loops(all_loops,loops_dir,loops_guide, reference_pixel);