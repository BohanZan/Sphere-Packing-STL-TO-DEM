function tests = testProfileHotspots
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.oldPath = path;
addpath(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(fileparts(mfilename('fullpath')),'helpers'));
end

function teardownOnce(testCase)
path(testCase.TestData.oldPath);
end

function testParityOnConvexConcaveAndCavity(testCase)
rng(741,'twister');
for kind = 1:3
    model = fixture(kind);
    context = spBuildContext(model,0.2,0,1e-9);
    points = context.lower + rand(600,3).*(context.upper-context.lower);
    points = [points; context.vertices; context.faceCentres; ...
        context.faceCentres+context.tolerance*[1 -1 1]; ...
        context.faceCentres-context.tolerance*[1 -1 1]]; %#ok<AGROW>
    for id = 1:size(points,1)
        verifyEqual(testCase,spExactPointInside(context,points(id,:)), ...
            referenceExactPointInside(context,points(id,:)));
    end
    if kind==3
        verifyFalse(testCase,spExactPointInside(context,[2 2 2]));
        verifyTrue(testCase,spExactPointInside(context,[0.5 0.5 0.5]));
    end
end
end

function testParityNearEdgesAcrossScales(testCase)
rng(912,'twister');
for scale = [1e-4 1 1e4]
    model = fixture(1);
    % Oblique, very thin projection and a large origin offset stress bounds.
    transform = [1 0.2 0.3;0.1 1e-5 0.1;0.2 0.4 1];
    model.vertices = scale*(model.vertices*transform+[1e3 -2e3 3e3]);
    % Radius only selects the hash here; retain all thin/edge points while
    % keeping this parity-only fixture inexpensive with the finer grid.
    context = spBuildContext(model,0.2*scale,0,1e-9);
    for id=1:size(model.faces,1)
        a=model.vertices(model.faces(id,1),:);
        b=model.vertices(model.faces(id,2),:);
        c=model.vertices(model.faces(id,3),:);
        t=context.tolerance;
        weights=[-t/2 0.3;1+t -t/2;0.5 0.5+t/2;0 0;1 0;0 1];
        points=a+weights(:,1).*(b-a)+weights(:,2).*(c-a);
        points(:,3)=points(:,3)+10*t;
        for k=1:size(points,1)
            verifyEqual(testCase,spExactPointInside(context,points(k,:)), ...
                referenceExactPointInside(context,points(k,:)));
        end
    end
end
end

function testRelaxGeometryAndSeededReproducibility(testCase)
for kind=1:3
    context=spBuildContext(fixture(kind),0.2,.2,1e-9);
    options=struct('maxAttempts',60,'maxCompressionSweeps',3, ...
        'shakeSweeps',2,'compressionTolerance',1e-7);
    for globalPass=[true false]
        % Build independent Map states: containers.Map is a handle object.
        rng(21,'twister');
        initial=spEmptyState(context,12);
        [initial,~]=spInitialPlacement(context,initial,linspace(.1,.2,12).',1,options);
        rng(91,'twister');
        % The former stepped solver can tunnel and shared one shake vector.
        % Its trajectory is no longer a correctness oracle. Check seeded
        % repeatability here and independent finite geometry below instead.
        reference=spRelax(context,initial,[0 0 -1],options,globalPass);
        expectedRng=rng;
        rng(21,'twister');
        actual=spEmptyState(context,12);
        [actual,~]=spInitialPlacement(context,actual,linspace(.1,.2,12).',1,options);
        rng(91,'twister');
        actual=spRelax(context,actual,[0 0 -1],options,globalPass);
        verifyEqual(testCase,actual.centres,reference.centres,'AbsTol',0);
        verifyEqual(testCase,actual.radii,reference.radii,'AbsTol',0);
        verifyEqual(testCase,rng,expectedRng);
        for id=1:actual.count
            p=actual.centres(id,:); r=actual.radii(id);
            verifyTrue(testCase,referenceExactPointInside(context,p));
            verifyFalse(testCase,referenceSphereHitsTriangles(context,p,r,1:size(context.faces,1)));
            others=id+1:actual.count;
            d=sqrt(sum((actual.centres(others,:)-p).^2,2));
            verifyTrue(testCase,all(d>=r+actual.radii(others)-2*context.tolerance));
        end
    end
end
end

function testLegacyContextAndObliqueContacts(testCase)
context=spBuildContext(fixture(2),.2,.2,1e-9);
fullContext=context;
% Older saved contexts have no new bounds: keep both query paths usable.
if isfield(context.ray,'xyLower')
    context.ray=rmfield(context.ray,{'xyLower','xyUpper'});
end
if isfield(context.ray,'xyBoundsCells')
    context.ray=rmfield(context.ray,'xyBoundsCells');
end
options=struct('maxAttempts',60,'maxCompressionSweeps',2, ...
    'shakeSweeps',2,'compressionTolerance',1e-7);
direction=[1 2 -3]; direction=direction/norm(direction);
rng(52,'twister');
expected=spEmptyState(context,8);
[expected,~]=spInitialPlacement(context,expected,.2*ones(8,1),1,options);
expected=spRelax(fullContext,expected,direction,options,true);
rng(52,'twister');
actual=spEmptyState(context,8);
[actual,~]=spInitialPlacement(context,actual,.2*ones(8,1),1,options);
actual=spRelax(context,actual,direction,options,true);
verifyEqual(testCase,actual.centres,expected.centres,'AbsTol',0);
end

function model=fixture(kind)
if kind==2
    polygon=[0 0;4 0;4 1;1 1;1 4;0 4];
    model.vertices=[polygon zeros(6,1);polygon 4*ones(6,1)];
    cap=[1 2 3;1 3 4;1 4 6;4 5 6];
    model.faces=[cap(:,[1 3 2]);cap+6];
    for id=1:6
        next=mod(id,6)+1;
        model.faces=[model.faces;id next next+6;id next+6 id+6]; %#ok<AGROW>
    end
else
    model.vertices=[0 0 0;4 0 0;4 4 0;0 4 0;0 0 4;4 0 4;4 4 4;0 4 4];
    model.faces=[1 3 2;1 4 3;5 6 7;5 7 8;1 2 6;1 6 5; ...
        2 3 7;2 7 6;3 4 8;3 8 7;4 1 5;4 5 8];
    if kind==3
        model.vertices=[model.vertices;1+model.vertices/2];
        model.faces=[model.faces;model.faces(:,[1 3 2])+8];
    end
end
end
