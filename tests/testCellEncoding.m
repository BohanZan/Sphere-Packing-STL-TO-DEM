function tests=testCellEncoding
tests=functiontests(localfunctions);
end
function testLiteralThreeDimensionalKeys(testCase)
verifyEqual(testCase,spCellKeys([1 1 1;3 1 1;1 2 1;1 1 2],[3 4 5]),uint64([1;3;4;13]));
end
function testLiteralTwoDimensionalKeys(testCase)
verifyEqual(testCase,spCellKeys([1 1;3 1;1 2;3 4],[3 4]),uint64([1;3;4;12]));
end
function testOverflowRetainsDistinctCoordinateKeys(testCase)
actual=spCellKeys([1 1 1;2 1 1;1 2 1;1 1 2],[2^32 2^32 2]);
verifyEqual(testCase,actual,{'1,1,1';'2,1,1';'1,2,1';'1,1,2'});
end
function testKeysBeyondDoublePrecisionRemainDistinct(testCase)
actual=spCellKeys([1 1 2;2 1 2],[2^27 2^27 2]);
verifyEqual(testCase,actual,[uint64(2^54)+uint64(1);uint64(2^54)+uint64(2)]);
end
function testEmptyQuery(testCase)
verifyEqual(testCase,spCellKeys(zeros(0,3),[3 4 5]),zeros(0,1,'uint64'));
end
function testRejectsInvalidCoordinateBeforeUnsignedCast(testCase)
verifyError(testCase,@() spCellKeys([0 1 1],[3 4 5]),'SpherePacking:InvalidCellIndex');
verifyError(testCase,@() spCellKeys([1 5 1],[3 4 5]),'SpherePacking:InvalidCellIndex');
end
function testSingleCoordinatesRetainIntegerPrecision(testCase)
indices=single([511 256 256;512 256 256]);
expected=[uint64(2^25)-uint64(1);uint64(2^25)];
verifyEqual(testCase,spCellKeys(indices,[512 256 256]),expected);
[~,spec]=spCellKeys([1 1 1],[512 256 256]);
verifyEqual(testCase,spCellKeys(indices,spec),expected);
end
