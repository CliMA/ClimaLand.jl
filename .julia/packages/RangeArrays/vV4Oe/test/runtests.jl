using RangeArrays
using Test

# Regular RangeMatrix tests

R = RangeMatrix(UnitRange{Int}[1:10, 11:20, 21:30, 31:40])

@test vec(collect(R)) == collect(1:40)
@test R[1,1] == 1

for i=1:40
    @test R[i] == i
end
for (i,j) in zip(eachindex(R), 1:length(R))
    @test R[i] == j
end

@test_throws BoundsError R[0]
@test_throws BoundsError R[41]

@test_throws BoundsError R[:, 0]
@test R[:, 1] == 1:10
@test R[:, 2] == 11:20
@test R[:, 3] == 21:30
@test R[:, 4] == 31:40
@test_throws BoundsError R[:, 5]

@test R[:, 1:2] == RangeMatrix(1:10, 11:20)
@test R[:, 1:4] == R[:,:] == RangeMatrix(1:10, 11:20, 21:30, 31:40)
@test R[:, 3:4] == RangeMatrix(21:30, 31:40)

# RepeatedRangeMatrix tests
R = RepeatedRangeMatrix(1:10, 0:10:30)

@test vec(collect(R)) == collect(1:40)
@test R[1,1] == 1

for i=1:40
    @test R[i] == i
end
for (i,j) in zip(eachindex(R), 1:length(R))
    @test R[i] == j
end

@test_throws BoundsError R[0]
@test_throws BoundsError R[41]

@test_throws BoundsError R[:, 0]
@test R[:, 1] == 1:10
@test R[:, 2] == 11:20
@test R[:, 3] == 21:30
@test R[:, 4] == 31:40
@test_throws BoundsError R[:, 5]

@test R[:, 1:2] == RangeMatrix(1:10, 11:20)
@test R[:, 1:4] == R[:,:] == RangeMatrix(1:10, 11:20, 21:30, 31:40)
@test R[:, 3:4] == RangeMatrix(21:30, 31:40)

A = 1:100
@test_throws BoundsError A[RangeMatrix(0:9, 20:29)]
@test_throws BoundsError A[RangeMatrix(20:29, 0:9)]
@test_throws BoundsError A[RangeMatrix(20:29, 92:101)]
@test_throws BoundsError A[RangeMatrix(92:101, 20:29)]
@test_throws BoundsError A[RangeMatrix(-100:100)]
@test A[RangeMatrix(1:10, 91:100)] == [1:10 91:100]

@test_throws BoundsError A[RepeatedRangeMatrix(1:10, [-1,33])]
@test_throws BoundsError A[RepeatedRangeMatrix(1:10, [33,-1])]
@test_throws BoundsError A[RepeatedRangeMatrix(1:10, [91,33])]
@test_throws BoundsError A[RepeatedRangeMatrix(1:10, [33,91])]
@test_throws BoundsError A[RepeatedRangeMatrix(-100:100, [0])]
@test A[RepeatedRangeMatrix(1:10, [0,90])] == [1:10 91:100]
