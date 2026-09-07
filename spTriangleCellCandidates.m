function [indices,distances,stats]=spTriangleCellCandidates(triangle,grid,radius,halo)
%SPTRIANGLECELLCANDIDATES Conservative finite-triangle distance-band cells.
% Enumerate projected columns, not the triangle's 3-D bounding volume.
% RADIUS is explicit: contact and occupancy grids use different distances.
% DISTANCES are conservative finite-triangle distance lower bounds; numerical
% feature ambiguity can retain extra cells but must never remove a contact.
h=grid.cellSize; counts=grid.cellCount;
tri=double(triangle)-grid.lower;
roundoff=128*eps(max([abs(triangle(:));abs(grid.lower(:));h]))+halo;
R=radius+roundoff;
lo=max(1,ceil((min(tri,[],1)-R)/h+.5)-1);
hi=min(counts,floor((max(tri,[],1)+R)/h+.5)+1);
indices=zeros(0,3); distances=zeros(0,1); stats=struct('enumerated',0);
if any(lo>hi), return; end
edges=tri([2 3 1],:)-tri;
[~,longest]=max(hypot(hypot(edges(:,1),edges(:,2)),edges(:,3)));
P=tri(longest,:); E=edges(longest,:); edgeLength=norm(E);
% A bounded tiny box is cheaper to evaluate in one MATLAB array operation.
% The 4096-cell budget is the same working-block size used below; large
% oblique/long-thin boxes still use surface/capsule enumeration, not volume.
if edgeLength==0 || prod(hi-lo+1)<=4096
    parts={boxIndices(lo,hi)};
else
    third=tri(mod(longest+1,3)+1,:);
    t=max(0,min(1,dot(third-P,E/edgeLength)/edgeLength));
    thickness=norm(third-P-t*E);
    if thickness<=h
        % A thin triangle is contained in its longest-edge capsule. Each
        % dominant-axis slice has a small transverse box, even if diagonal.
        S=R+thickness; [~,axis]=max(abs(E)); other=setdiff(1:3,axis);
        slices=lo(axis):hi(axis); parts=cell(numel(slices),1);
        for k=1:numel(slices)
            centre=(slices(k)-.5)*h;
            ts=sort(([centre-S centre+S]-P(axis))/E(axis));
            t0=max(0,ts(1)); t1=min(1,ts(2));
            if t0>t1, continue; end
            ends=P+[t0;t1]*E;
            lower=max(lo,ceil((min(ends,[],1)-S)/h+.5)-1);
            upper=min(hi,floor((max(ends,[],1)+S)/h+.5)+1);
            lower(axis)=slices(k); upper(axis)=slices(k);
            if any(lower(other)>upper(other)), continue; end
            parts{k}=boxIndices(lower,upper);
        end
    else
        ab=tri(2,:)-tri(1,:); ac=tri(3,:)-tri(1,:);
        scale=max(abs([ab ac])); m=cross(ab/scale,ac/scale);
        [~,axis]=max(abs(m)); m=m/max(abs(m));
        other=setdiff(1:3,axis);
        % Use the actual vertex support interval, including normal roundoff.
        support=(tri-tri(1,:))*m.';
        slab=[min(support)-R*norm(m),max(support)+R*norm(m)];
        nx=hi(other(1))-lo(other(1))+1;
        ny=hi(other(2))-lo(other(2))+1;
        parts=cell(ceil(nx*ny/4096),1);
        for block=1:numel(parts)
            number=((block-1)*4096:min(block*4096,nx*ny)-1).';
            xy=[lo(other(1))+rem(number,nx),lo(other(2))+floor(number/nx)];
            centres=(xy-.5)*h;
            base=(centres-tri(1,other))*m(other).';
            limits=(slab-base)/m(axis)+tri(1,axis);
            lower=max(lo(axis),ceil(min(limits,[],2)/h+.5)-1);
            upper=min(hi(axis),floor(max(limits,[],2)/h+.5)+1);
            lengths=max(0,upper-lower+1);
            rows=repelem((1:numel(lengths)).',lengths);
            starts=cumsum(lengths)-lengths;
            positions=(1:sum(lengths)).'-repelem(starts,lengths)-1;
            ijk=zeros(numel(rows),3);
            ijk(:,other)=xy(rows,:); ijk(:,axis)=lower(rows)+positions;
            parts{block}=ijk;
        end
    end
end
distanceParts=cell(size(parts));
for block=1:numel(parts)
    ijk=parts{block}; if isempty(ijk), continue; end
    stats.enumerated=stats.enumerated+size(ijk,1);
    points=(ijk-.5)*h;
    d=finiteDistance(points,tri);
    keep=d<=R;
    parts{block}=ijk(keep,:); distanceParts{block}=d(keep);
end
indices=vertcat(parts{:}); distances=vertcat(distanceParts{:});
if isempty(indices), indices=zeros(0,3); distances=zeros(0,1); end
end

function indices=boxIndices(lo,hi)
[x,y,z]=ndgrid(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3));
indices=[x(:) y(:) z(:)];
end

function distance=finiteDistance(points,tri)
% T lies in projection_n(T) times its vertex support interval for ANY unit n.
% Distance to that superset is a lower bound even when cross products lose
% orthogonality on very skinny faces. Keep 3-D projected edge distances.
relative=tri-tri(1,:); q=points-tri(1,:);
scale=max(abs(relative(:)));
n=zeros(1,3);
if scale>0, n=cross(relative(2,:)/scale,relative(3,:)/scale); end
nLength=norm(n); slab=zeros(size(points,1),1);
if nLength>0
    n=n/nLength; support=relative*n.'; height=q*n.';
    slab=max(0,max(min(support)-height,height-max(support)));
    relative=relative-support*n; q=q-height*n;
end
projected=inf(size(points,1),1);
for edge=[1 2;2 3;3 1].'
    e=relative(edge(2),:)-relative(edge(1),:); v=q-relative(edge(1),:);
    length=norm(e);
    if length>0
        unit=e/length; v=v-max(0,min(length,v*unit.'))*unit;
    end
    projected=min(projected,hypot(hypot(v(:,1),v(:,2)),v(:,3)));
end
if nLength>0
    [~,axis]=max(abs(n)); other=setdiff(1:3,axis);
    e=relative(2,other)/scale; f=relative(3,other)/scale; v=q(:,other)/scale;
    det=e(1)*f(2)-e(2)*f(1);
    if det~=0
        alpha=(v(:,1)*f(2)-v(:,2)*f(1))/det;
        beta=(e(1)*v(:,2)-e(2)*v(:,1))/det;
        error=128*eps*max(1,(sum(abs(e))+sum(abs(f)))^2/abs(det));
        inside=alpha>=-error & beta>=-error & alpha+beta<=1+error;
        projected(inside)=0;
    end
end
distance=hypot(projected,slab);
end
