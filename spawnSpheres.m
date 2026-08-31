function [ assembly, masses, totVolSph, mI ] = populateSpheres( fileName, nNodes, scaleFactor, smoothFact, writeTec, ...
                                                                writeLammpsTemp, lammpsType, writeStlScaled, rho, findSTLmI )

                                                                [vertices,faces,~,~] = stlRead(fileName);

vertices = vertices*scaleFactor;

%find a bounding box
xlow=min(vertices(:,1));
ylow=min(vertices(:,2));
zlow=min(vertices(:,3));

xhigh=max(vertices(:,1));
yhigh=max(vertices(:,2));
zhigh=max(vertices(:,3));

Lx = xhigh-xlow;
Ly = yhigh-ylow;
Lz = zhigh-zlow;

disp(['Bounding Box Dimensions Lx=', num2str(Lx), '; Ly=', num2str(Ly), ...
                                                  '; Lz=', num2str(Lz)]);

%Calculate some parameters
rSph = Lx/2/nNodes; %Radius of potential spherical cells

nRayY = floor(Ly/Lx*nNodes)+1; %Nomber of place holder lines in y-direction
nRayZ = floor(Lz/Lx*nNodes)+1; %Nomber of place holder lines in z-direction

%Determine the direction and initial origin of eah ray, direction is 
%arbitrarily set to +x, origins are distributed unifromely on x=xlow plane
%This should not be changed
initCoords = zeros(3,nRayY,nRayZ);
initDirs = zeros(3,nRayY,nRayZ);
for k = 1:nRayZ
    for j = 1:nRayY
        initCoords(:,j,k)=[xlow-1e13*eps(xlow); ylow+(j-1)*2*rSph; zlow+(k-1)*2*rSph]; 
        initDirs(:,j,k)=[Lx;0;0]; 
    end
end

%generate spheres along each ray (+x)
initSphCoords = zeros(nNodes+1,nRayY,nRayZ);
for k = 1:nRayZ
    for j = 1:nRayY
        for n = 1:nNodes+1
            initSphCoords(n,j,k)=xlow+(n-1)*2*rSph; 
        end
    end
end

%Just for testing a very simple crude Monte Carlo algorithm to calculte mI
%of the original STL, set findSTLmI to zero since this may take very long
%specially if the stl is not convex
[stlVol,~] = stlVolume(vertices',faces');
%stlMass = stlVol*rho;
if(findSTLmI) 
    [STLmI, STLCoM] = calculateSTLmI(vertices,faces,stlVol,rho,10000);
end

%For each ray determine the intersection
sz = size(faces);
nface= sz(1);

%intersections of each ray witht the collection of trinagles in +x
maxIntersectPerRay = 10; %Increase of crashed
allIntersects = zeros(maxIntersectPerRay, nRayY, nRayZ);
nIntersectRay = zeros(nRayY,nRayZ);

for k = 1:nRayZ
    for j = 1:nRayY
        direction = initDirs(:,j,k);
        origin    = initCoords(:,j,k);
        ncIntersect = 0;
        for facet=1:nface
            vert = vertices(faces(facet,:),:);
            [flag, ~, ~, t] = rayTriangleIntersection(origin, direction, vert(1,:)', vert(2,:)', vert(3,:)', eps(Lx));    
            if(flag)
                intersection = origin + t*direction;
                ncIntersect = ncIntersect+1;
                allIntersects(ncIntersect,j,k) = intersection(1);
            end
        end
        nIntersectRay(j,k) = ncIntersect;
        allIntersects(1:ncIntersect,j,k) = sort(allIntersects(1:ncIntersect,j,k));
    end
end

%bFlagSphere = 0: surface node, <7: inernal node, =7 
bFlagSphere = zeros(nNodes+1,nRayY,nRayZ);

for k = 1:nRayZ
    for j = 1:nRayY
        if(nIntersectRay(j,k)) %otherwise bFlagSphere is already set to 0
            for nn = 1:2:nIntersectRay(j,k)
                for n = 1:nNodes+1
                    if(allIntersects(nn,j,k) < initSphCoords(n,j,k) && ...
                            allIntersects(nn+1,j,k) > initSphCoords(n,j,k) )
                        bFlagSphere(n,j,k) = bFlagSphere(n,j,k) + 1;
                        %check to see if boundary
                        if(n ~=1)
                            if(bFlagSphere(n-1,j,k))
                                bFlagSphere(n-1,j,k) = bFlagSphere(n-1,j,k) + 1;
                                bFlagSphere(n,j,k) = bFlagSphere(n,j,k) + 1;
                            end
                        end
                        if(n ~= nNodes+1)
                            if(bFlagSphere(n+1,j,k))
                                bFlagSphere(n+1,j,k) = bFlagSphere(n+1,j,k) + 1;
                                bFlagSphere(n,j,k) = bFlagSphere(n,j,k) + 1;
                            end
                        end
                        
                        if(j~=1)
                            if(bFlagSphere(n,j-1,k))
                                bFlagSphere(n,j-1,k) = bFlagSphere(n,j-1,k) + 1;
                                bFlagSphere(n,j,k) = bFlagSphere(n,j,k) + 1;
                            end
                        end
                        
                        if(j ~= nRayY)
                            if(bFlagSphere(n,j+1,k))
                                bFlagSphere(n,j+1,k) = bFlagSphere(n,j+1,k) + 1;
                                bFlagSphere(n,j,k) = bFlagSphere(n,j,k) + 1;
                            end
                        end
                        
                        if(k ~= 1)
                            if(bFlagSphere(n,j,k-1))
                                bFlagSphere(n,j,k-1) = bFlagSphere(n,j,k-1) + 1;
                                bFlagSphere(n,j,k) = bFlagSphere(n,j,k) + 1;
                            end
                        end
                        
                        if(k ~= nRayZ)
                            if(bFlagSphere(n,j,k+1))
                                bFlagSphere(n,j,k+1) = bFlagSphere(n,j,k+1) + 1;
                                bFlagSphere(n,j,k) = bFlagSphere(n,j,k) + 1;
                            end
                        end
                        
                    end
                    
                end
            end
        end
        
    end
end

