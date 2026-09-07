function tests=testStaticGrid
tests=functiontests(localfunctions);
end
function testStableSortedFaceLists(testCase)
grid=spBuildStaticGrid(uint64([9;2;9;2;9]),[8;7;3;1;3],'uint64');
verifyEqual(testCase,spGridLookup(grid,{uint64(9)}),[3 8]);
verifyEqual(testCase,spGridLookup(grid,{uint64(2)}),[1 7]);
verifyEqual(testCase,spGridLookup(grid,{uint64(9),uint64(99),uint64(2)}),[3 8 1 7]);
verifyEqual(testCase,numel(grid.offsets),numel(grid.occupiedKeys)+1);
end
function testEmptyAndWideKeys(testCase)
empty=spBuildStaticGrid(zeros(0,1,'uint64'),zeros(0,1),'uint64');
verifyEmpty(testCase,spGridLookup(empty,{uint64(1)}));
base=uint64(2^54);
grid=spBuildStaticGrid([base+uint64(1);base+uint64(2)],[2;1],'uint64');
verifyEqual(testCase,spGridLookup(grid,{base+uint64(1)}),2);
verifyEqual(testCase,spGridLookup(grid,{base+uint64(2)}),1);
verifyEqual(testCase,numel(grid.offsets),3);
end
function testOverflowCoordinateKeys(testCase)
grid=spBuildStaticGrid({'1,1,1';'2,1,1';'1,1,1'},[5;7;2],'char');
verifyEqual(testCase,spGridLookup(grid,{'1,1,1','2,1,1'}),[2 5 7]);
end
function testLegacyMapLookup(testCase)
map=containers.Map('KeyType','uint64','ValueType','any');
map(uint64(2))=[1 7];
verifyEqual(testCase,spGridLookup(map,{uint64(9),uint64(2)}),[1 7]);
end
function testCharacterKeysWithRowFaceIds(testCase)
grid=spBuildStaticGrid({'a';'a';'b';'b'},[1 2 3 4],'char');
verifyEqual(testCase,spGridLookup(grid,{'a'}),[1 2]);
verifyEqual(testCase,spGridLookup(grid,{'b'}),[3 4]);
end
function testScalarCharacterKey(testCase)
grid=spBuildStaticGrid({'1,1,1';'2,1,1'},[2;7],'char');
verifyEqual(testCase,spGridLookup(grid,'1,1,1'),2);
end
