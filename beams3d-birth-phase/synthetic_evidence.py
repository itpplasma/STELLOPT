import sys,json,tempfile
from pathlib import Path
import numpy as np,unyt
from a5py import Ascot
mode,fn,out=sys.argv[1:]
mrk=Ascot(fn).data.marker.active.read()
tmp=tempfile.TemporaryDirectory();a=Ascot(str(Path(tmp.name)/'uniform.h5'),create=True)
a.data.create_input('B_TC',bxyz=np.array([0.,2.,0.]),jacobian=np.zeros((3,3)),rhoval=.5)
a.input_init(bfield=True);a.simulation_initmarkers(**mrk)
s=a._inistate;n=len(s) if hasattr(s,'__len__') else mrk['n'];n=mrk['n']
get=lambda k:np.array([getattr(s[i],k) for i in range(n)])
r,p,z=get('rprt'),get('phiprt'),get('zprt');xyz=np.column_stack([r*np.cos(p),r*np.sin(p),z])
pp=np.column_stack([get(k) for k in ['p_r','p_phi','p_z']]);mass=get('mass');vel=pp/np.sqrt(mass**2+np.sum(pp**2,axis=1)/299792458.**2)[:,None]
truthphi=np.array([.1,.4,.8]);truthr=np.array([10.5,10.4,10.6]);truthz=np.array([0,.1,-.1]);truthx=np.column_stack([truthr*np.cos(truthphi),truthr*np.sin(truthphi),truthz])
truthv=np.array([[1e6,2e6,.4e6],[-1e6,1.5e6,-.3e6],[.5e6,-2e6,.6e6]])
# Convert cylindrical physical velocity to Cartesian before comparing bases.
vxyz=np.column_stack([vel[:,0]*np.cos(p)-vel[:,1]*np.sin(p),vel[:,0]*np.sin(p)+vel[:,1]*np.cos(p),vel[:,2]])
tvxyz=np.column_stack([truthv[:,0]*np.cos(truthphi)-truthv[:,1]*np.sin(truthphi),truthv[:,0]*np.sin(truthphi)+truthv[:,1]*np.cos(truthphi),truthv[:,2]])
massref=(mrk['mass']*unyt.amu).to('kg').v
gamma=1/np.sqrt(1-np.sum(truthv**2,axis=1)/299792458.**2);pref=massref*gamma*np.linalg.norm(truthv,axis=1)
report=dict(mode=mode,marker_ids=[1,2,3],position_error_mm=(np.linalg.norm(xyz-truthx,axis=1)*1000).tolist(),velocity_error_km_s=(np.linalg.norm(vxyz-tvxyz,axis=1)/1000).tolist(),momentum_ratio=(np.linalg.norm(pp,axis=1)/pref).tolist(),errors=get('err').astype(int).tolist())
Path(out).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report));a.input_free()
