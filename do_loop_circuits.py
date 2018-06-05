

# Step 1: Glob the intf_all
# Step 2: Form as many loops as possible
# Step 3: Make images of all the possible loops. 

import numpy as np 
import matplotlib.pyplot as plt 
import glob



def form_all_loops():
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
	# print(graph);

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




def show_images(all_loops, loops_dir):
	return;




if __name__=="__main__":
	loops_dir="Phase_Circuits/"
	all_loops=form_all_loops();
	show_images(all_loops,loops_dir);