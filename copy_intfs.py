import glob
import subprocess


intf_all="intf_all_save_for_posterity"
source="phase.grd"
destination="intf_all_remote/"



files=glob.glob(intf_all+"/???????_???????/"+source);
print(files);

for item in files:
	namestring=item.split('/')[-2];
	print(namestring);
	subprocess.call(['mkdir','-p',destination+namestring],shell=False);
	subprocess.call(['cp',item,destination+namestring+'/phase_raw.grd'],shell=False);