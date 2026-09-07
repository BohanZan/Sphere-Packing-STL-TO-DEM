"""Independent NumPy/SciPy finite-triangle and sphere checks (no C++ grid reuse)."""
from pathlib import Path
import argparse
import json
import numpy as np
from scipy.spatial import cKDTree


def distances(point, triangles):
    a,b,c=triangles[:,0],triangles[:,1],triangles[:,2]
    ab,ac=b-a,c-a
    n=np.cross(ab,ac)
    nsq=np.einsum('ij,ij->i',n,n)
    h=np.einsum('ij,ij->i',point-a,n)
    projected=point-n*np.divide(h,nsq,out=np.zeros_like(h),where=nsq>0)[:,None]
    signs=[]
    minimum=np.full(len(a),np.inf)
    for v,w in [(a,b),(b,c),(c,a)]:
        e=w-v
        denominator=np.einsum('ij,ij->i',e,e)
        t=np.clip(np.divide(np.einsum('ij,ij->i',point-v,e),denominator,
                            out=np.zeros_like(denominator),where=denominator>0),0,1)
        minimum=np.minimum(minimum,np.linalg.norm(point-v-t[:,None]*e,axis=1))
        signs.append(np.einsum('ij,ij->i',np.cross(e,projected-v),n))
    inside=np.all(np.array(signs)>=0,axis=0)&(nsq>0)
    minimum[inside]=np.minimum(minimum[inside],np.abs(h[inside])/np.sqrt(nsq[inside]))
    return minimum


def validate(triangle_file,sphere_file,tolerance=1e-7):
    triangles=np.loadtxt(triangle_file,delimiter=',').reshape(-1,3,3)
    spheres=np.loadtxt(sphere_file,delimiter=',',skiprows=1,ndmin=2)
    p,r=spheres[:,1:4],spheres[:,4]
    pairs=cKDTree(p).query_pairs(2*max(r)+tolerance,output_type='ndarray')
    pair_gaps=np.linalg.norm(p[pairs[:,0]]-p[pairs[:,1]],axis=1)-r[pairs[:,0]]-r[pairs[:,1]]
    min_pair=float(np.min(pair_gaps)) if len(pair_gaps) else 0.
    # Bounding spheres around each triangle form a separate conservative broad phase.
    centres=triangles.mean(axis=1)
    bounds=np.max(np.linalg.norm(triangles-centres[:,None],axis=2),axis=1)
    positive=np.maximum(bounds,np.finfo(float).tiny)
    bins=np.floor(np.log2(positive)).astype(int)
    groups=[]
    for value in np.unique(bins):
        ids=np.flatnonzero(bins==value)
        groups.append((ids,float(max(bounds[ids])),cKDTree(centres[ids]),cKDTree(centres[ids,:2])))
    min_wall=np.inf
    outside=[]
    for i,(point,radius) in enumerate(zip(p,r)):
        candidates=[];ray_candidates=[]
        for ids,bound,tree,ray_tree in groups:
            candidates.extend(ids[tree.query_ball_point(point,radius+bound+tolerance)])
            ray_candidates.extend(ids[ray_tree.query_ball_point(point[:2],bound+tolerance)])
        if candidates:
            gap=float(min(distances(point,triangles[candidates]))-radius)
            min_wall=min(min_wall,gap)
        tri=triangles[ray_candidates]
        a=tri[:,0];ab=tri[:,1]-a;ac=tri[:,2]-a;ap=point-a
        det=ab[:,0]*ac[:,1]-ab[:,1]*ac[:,0]
        keep=np.abs(det)>1e-20
        a,ab,ac,ap,det=a[keep],ab[keep],ac[keep],ap[keep],det[keep]
        u=(ap[:,0]*ac[:,1]-ap[:,1]*ac[:,0])/det
        v=(ab[:,0]*ap[:,1]-ab[:,1]*ap[:,0])/det
        keep=(u>=-1e-12)&(v>=-1e-12)&(u+v<=1+1e-12)
        hits=a[:,2]+u*ab[:,2]+v*ac[:,2]
        z=np.sort(hits[keep&(hits<point[2]-1e-10)])
        count=(1+np.count_nonzero(np.diff(z)>1e-9)) if len(z) else 0
        if count%2==0:outside.append(i+1)
    result=dict(spheres=len(r),minPairGap=min_pair,minWallGap=float(min_wall),outsideIds=outside,
                tolerance=tolerance,method='Independent cKDTree triangle bounding spheres + finite plane/segment distances; independent ray parity')
    print(json.dumps(result,indent=2),flush=True)
    assert min_pair>=-tolerance and min_wall>=-tolerance and not outside,'Geometry validation failed'
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--triangles',type=Path,required=True)
    parser.add_argument('--spheres',type=Path,required=True);parser.add_argument('--report',type=Path)
    args=parser.parse_args();result=validate(args.triangles,args.spheres)
    if args.report:args.report.write_text(json.dumps(result,indent=2)+'\n')
