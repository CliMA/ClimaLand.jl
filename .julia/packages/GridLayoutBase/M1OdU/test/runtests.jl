using GridLayoutBase
using GridLayoutBase: GridSubposition
using GridLayoutBase: offsets
using GridLayoutBase: RectSides
using Test
using Observables

include("debugrect.jl")

# run this to make bboxes easier to read: Base.show(io::IO, b::GridLayoutBase.GeometryBasics.Rect2) = print(io, "BBox($(left(b)), $(right(b)), $(bottom(b)), $(top(b)))")

@testset "GridLayout Zero Outside AlignMode" begin
    bbox = BBox(0, 1000, 0, 1000)
    layout = GridLayout(bbox = bbox, alignmode = Outside(0))
    dr = layout[1, 1] = DebugRect()

    @test computedbboxobservable(dr)[] == bbox

    dr.topprot[] = 100
    @test computedbboxobservable(dr)[] == BBox(0, 1000, 0, 900)
    dr.bottomprot[] = 100
    @test computedbboxobservable(dr)[] == BBox(0, 1000, 100, 900)
    dr.leftprot[] = 100
    @test computedbboxobservable(dr)[] == BBox(100, 1000, 100, 900)
    dr.rightprot[] = 100
    @test computedbboxobservable(dr)[] == BBox(100, 900, 100, 900)

    dr2 = layout[1, 2] = DebugRect()
    @test nrows(layout) == 1 && ncols(layout) == 2
    colgap!(layout, 1, Fixed(0))

    @test computedbboxobservable(dr)[].widths == computedbboxobservable(dr2)[].widths == Float32[400.0, 800.0]
end

@testset "GridLayout Outside AlignMode" begin
    bbox = Observable(BBox(0, 1000, 0, 1000)) # just as observable for bbox conversion test
    layout = GridLayout(bbox = bbox, alignmode = Outside(100, 200, 50, 150))
    dr = layout[1, 1] = DebugRect()

    @test computedbboxobservable(dr)[] == BBox(100, 800, 50, 850)

    dr.topprot[] = 100
    @test computedbboxobservable(dr)[] == BBox(100, 800, 50, 750)
    dr.bottomprot[] = 100
    @test computedbboxobservable(dr)[] == BBox(100, 800, 150, 750)
    dr.leftprot[] = 100
    @test computedbboxobservable(dr)[] == BBox(200, 800, 150, 750)
    dr.rightprot[] = 100
    @test computedbboxobservable(dr)[] == BBox(200, 700, 150, 750)
end

@testset "GridLayout Inside AlignMode" begin
    bbox = BBox(0, 1000, 0, 1000)
    layout = GridLayout(bbox = bbox, alignmode = Inside())
    dr = layout[1, 1] = DebugRect()

    @test computedbboxobservable(dr)[] == BBox(0, 1000, 0, 1000)

    dr.topprot[] = 100
    @test computedbboxobservable(dr)[] == BBox(0, 1000, 0, 1000)
    dr.bottomprot[] = 100
    @test computedbboxobservable(dr)[] == BBox(0, 1000, 0, 1000)
    dr.leftprot[] = 100
    @test computedbboxobservable(dr)[] == BBox(0, 1000, 0, 1000)
    dr.rightprot[] = 100
    @test computedbboxobservable(dr)[] == BBox(0, 1000, 0, 1000)
end

@testset "GridLayout Mixed AlignMode" begin
    bbox = BBox(0, 1000, 0, 1000)
    layout = GridLayout(bbox = bbox, alignmode = Mixed(left = 0, top = 100))
    dr = layout[1, 1] = DebugRect()

    @test GridLayoutBase.protrusion(layout, Left()) == 0
    @test GridLayoutBase.protrusion(layout, Right()) == 0
    @test GridLayoutBase.protrusion(layout, Bottom()) == 0
    @test GridLayoutBase.protrusion(layout, Top()) == 0

    @test computedbboxobservable(dr)[] == BBox(0, 1000, 0, 900)

    dr.topprot[] = 100
    @test GridLayoutBase.protrusion(layout, Top()) == 0
    @test computedbboxobservable(dr)[] == BBox(0, 1000, 0, 800)

    dr.bottomprot[] = 100
    @test GridLayoutBase.protrusion(layout, Bottom()) == 100
    @test computedbboxobservable(dr)[] == BBox(0, 1000, 0, 800)

    dr.leftprot[] = 100
    @test GridLayoutBase.protrusion(layout, Left()) == 0
    @test computedbboxobservable(dr)[] == BBox(100, 1000, 0, 800)

    dr.rightprot[] = 100
    @test GridLayoutBase.protrusion(layout, Right()) == 100
    @test computedbboxobservable(dr)[] == BBox(100, 1000, 0, 800)

    layout.alignmode[] = Mixed(
        left = Protrusion(50),
        right = Protrusion(60),
        bottom = Protrusion(70),
        top = Protrusion(80),
    )

    # set layout to forced protrusion alignment
    @test GridLayoutBase.protrusion(layout, Left()) == 50
    @test GridLayoutBase.protrusion(layout, Right()) == 60
    @test GridLayoutBase.protrusion(layout, Bottom()) == 70
    @test GridLayoutBase.protrusion(layout, Top()) == 80
    # bb of dr only depends on its protrusions now
    @test computedbboxobservable(dr)[] == BBox(100, 1000, 0, 800)

    # force different protrusions for dr
    dr.alignmode[] = Mixed(
        left = Protrusion(50),
        right = Protrusion(60),
        bottom = Protrusion(70),
        top = Protrusion(80),
    )
    # bb should follow
    @test computedbboxobservable(dr)[] == BBox(50, 940, 70, 920)
end

@testset "alignmodes and reported dimensions" begin
    gl = GridLayout()
    dr = gl[1, 1] = DebugRect(width = 200, height = 300, alignmode = Mixed(left = 1, right = 2, bottom = 3, top = 4))
    @test reporteddimensionsobservable(dr)[].inner == (203, 307)
    @test reporteddimensionsobservable(dr)[].outer == RectSides{Float32}(0, 0, 0, 0)
    dr.alignmode[] = Inside()
    @test reporteddimensionsobservable(dr)[].inner == (200, 300)
    @test reporteddimensionsobservable(dr)[].outer == RectSides{Float32}(0, 0, 0, 0)
    dr.alignmode[] = Mixed(left = Protrusion(10), right = Protrusion(10), bottom = Protrusion(10), top = Protrusion(10))
    @test reporteddimensionsobservable(dr)[].inner == (200, 300)
    @test reporteddimensionsobservable(dr)[].outer == RectSides{Float32}(10, 10, 10, 10)
end

@testset "assigning content to protrusions" begin

    bbox = BBox(0, 1000, 0, 1000)
    layout = GridLayout(bbox = bbox, alignmode = Outside(0))
    subgl = layout[1, 1] = GridLayout()

    subgl[1, 1, Left()] = DebugRect(width = Fixed(100))
    @test GridLayoutBase.protrusion(subgl, Left()) == 100

    subgl[1, 1, Top()] = DebugRect(height = 50)
    @test GridLayoutBase.protrusion(subgl, Top()) == 50

    subgl[1, 1, Right()] = DebugRect(width = 120)
    @test GridLayoutBase.protrusion(subgl, Right()) == 120

    subgl[1, 1, Bottom()] = DebugRect(height = 40)
    @test GridLayoutBase.protrusion(subgl, Bottom()) == 40

    subgl[1, 1, TopLeft()] = DebugRect(width = 200, height = 200)
    @test GridLayoutBase.protrusion(subgl, Left()) == 200
    @test GridLayoutBase.protrusion(subgl, Top()) == 200

    subgl[1, 1, TopRight()] = DebugRect(width = 210, height = 210)
    @test GridLayoutBase.protrusion(subgl, Right()) == 210
    @test GridLayoutBase.protrusion(subgl, Top()) == 210

    subgl[1, 1, BottomRight()] = DebugRect(width = 220, height = 220)
    @test GridLayoutBase.protrusion(subgl, Right()) == 220
    @test GridLayoutBase.protrusion(subgl, Bottom()) == 220

    subgl[1, 1, BottomLeft()] = DebugRect(width = 230, height = 230)
    @test GridLayoutBase.protrusion(subgl, Left()) == 230
    @test GridLayoutBase.protrusion(subgl, Bottom()) == 230

    # issue 66
    bbox = BBox(0, 1000, 0, 1000)
    layout = GridLayout(bbox = bbox, alignmode = Outside(0))
    subgl = layout[1, 1] = GridLayout()
    dr = subgl[1, 1, GridLayoutBase.Outer()] = DebugRect()
    @test GridLayoutBase.protrusion(subgl, Left()) == 0
    @test GridLayoutBase.protrusion(subgl, Right()) == 0
    @test GridLayoutBase.protrusion(subgl, Bottom()) == 0
    @test GridLayoutBase.protrusion(subgl, Top()) == 0

    @test width(computedbboxobservable(dr)[]) ≈ 1000
    @test height(computedbboxobservable(dr)[]) ≈ 1000

    subgl[1, 1] = DebugRect(width = 200, height = 300, leftprot = 10, rightprot = 20, bottomprot = 30, topprot = 40)
    @test GridLayoutBase.protrusion(subgl, Left()) == 10
    @test GridLayoutBase.protrusion(subgl, Right()) == 20
    @test GridLayoutBase.protrusion(subgl, Bottom()) == 30
    @test GridLayoutBase.protrusion(subgl, Top()) == 40

    @test width(computedbboxobservable(dr)[]) ≈ 200 + 10 + 20
    @test height(computedbboxobservable(dr)[]) ≈ 300 + 30 + 40
end


@testset "resizing through indexing out of range and trim!" begin

    bbox = BBox(0, 1000, 0, 1000)
    layout = GridLayout(bbox = bbox, alignmode = Outside(0))

    dr = layout[1, 1] = DebugRect()
    @test size(layout) == (1, 1)
    @test offsets(layout) == (0, 0)

    layout[1, 2] = dr
    @test size(layout) == (1, 2)
    @test offsets(layout) == (0, 0)

    layout[3, 2] = dr
    @test size(layout) == (3, 2)
    @test offsets(layout) == (0, 0)

    layout[4, 4] = dr
    @test size(layout) == (4, 4)
    @test offsets(layout) == (0, 0)

    layout[0, 1] = dr
    @test size(layout) == (5, 4)
    @test offsets(layout) == (-1, 0)

    layout[1, 0] = dr
    @test size(layout) == (5, 5)
    @test offsets(layout) == (-1, -1)

    layout[-1, -1] = dr
    @test size(layout) == (6, 6)
    @test offsets(layout) == (-2, -2)

    layout[3, 3] = dr
    trim!(layout)
    @test size(layout) == (1, 1)
    @test offsets(layout) == (-2, -2)
    # reset offsets to zero
    layout.offsets = (0, 0)

    layout[2:3, 4:5] = dr
    @test size(layout) == (3, 5)

    trim!(layout)
    @test size(layout) == (2, 2)
end


@testset "manually deleting rows / cols" begin
    layout = GridLayout(4, 4)

    deletecol!(layout, 2)
    @test size(layout) == (4, 3)

    deleterow!(layout, 3)
    @test size(layout) == (3, 3)

    deleterow!(layout, 1)
    @test size(layout) == (2, 3)

    deletecol!(layout, 1)
    @test size(layout) == (2, 2)

    deleterow!(layout, 2)
    @test size(layout) == (1, 2)

    deletecol!(layout, 2)
    @test size(layout) == (1, 1)

    @test_throws ErrorException deletecol!(layout, 2)
    @test_throws ErrorException deleterow!(layout, 2)
    @test_throws ErrorException deletecol!(layout, 1)
    @test_throws ErrorException deleterow!(layout, 1)

    dr = layout[1, 2] = DebugRect()
    @test length(layout.content) == 1
    deletecol!(layout, 2)
    @test isempty(layout.content)

    dr = layout[2, 1] = DebugRect()
    @test length(layout.content) == 1
    deleterow!(layout, 2)
    @test isempty(layout.content)
end

@testset "inserting rows / cols" begin

    layout = GridLayout(4, 4)
    dr1 = layout[1, 1] = DebugRect()
    dr2 = layout[1, 2:3] = DebugRect()
    dr3 = layout[1, 4] = DebugRect()
    dr4 = layout[2:3, 1] = DebugRect()
    dr5 = layout[4, 1] = DebugRect()
    dr6 = layout[2:3, 2:3] = DebugRect()
    dr7 = layout[1, 3] = DebugRect()
    dr8 = layout[3, 1] = DebugRect()

    @test_throws ErrorException insertcols!(layout, 0, 1)
    @test_throws ErrorException insertcols!(layout, 5, 1)
    @test_throws ErrorException insertrows!(layout, 0, 1)
    @test_throws ErrorException insertrows!(layout, 5, 1)

    insertcols!(layout, 3, 2)
    @test gridcontent(dr1).span.cols == 1:1
    @test gridcontent(dr2).span.cols == 2:5
    @test gridcontent(dr3).span.cols == 6:6

    insertrows!(layout, 3, 2)
    @test gridcontent(dr1).span.rows == 1:1
    @test gridcontent(dr4).span.rows == 2:5
    @test gridcontent(dr5).span.rows == 6:6

    @test gridcontent(dr6).span.cols == 2:5
    @test gridcontent(dr6).span.rows == 2:5

    @test gridcontent(dr7).span.cols == 5:5
    @test gridcontent(dr8).span.rows == 5:5
end

@testset "setting col and row sizes and gaps" begin
    bbox = BBox(0, 1000, 0, 1000)
    layout = GridLayout(3, 3, bbox = bbox, alignmode = Outside(0))

    colsize!(layout, 1, Fixed(10))
    @test layout.colsizes[1] == Fixed(10)

    colsize!(layout, 2, 20)
    @test layout.colsizes[2] == Fixed(20)

    rowsize!(layout, 1, Relative(0.2))
    @test layout.rowsizes[1] == Relative(0.2)

    rowsize!(layout, 2, 15.3)
    @test layout.rowsizes[2] == Fixed(15.3)

    @test_throws ErrorException colsize!(layout, 4, Auto())
    @test_throws ErrorException rowsize!(layout, 0, Auto())


    colgap!(layout, 1, Fixed(10))
    @test layout.addedcolgaps[1] == Fixed(10)
    rowgap!(layout, 2, Relative(0.3))
    @test layout.addedrowgaps[2] == Relative(0.3)

    colgap!(layout, 10)
    @test all(layout.addedcolgaps .== Ref(Fixed(10)))
    rowgap!(layout, 20)
    @test all(layout.addedrowgaps .== Ref(Fixed(20)))

    colgap!(layout, Fixed(30))
    @test all(layout.addedcolgaps .== Ref(Fixed(30)))
    rowgap!(layout, Fixed(40))
    @test all(layout.addedrowgaps .== Ref(Fixed(40)))

    @test_throws ErrorException colgap!(layout, 10, Fixed(10))
    @test_throws ErrorException rowgap!(layout, 10, Fixed(10))
end

@testset "some constructors" begin
    @test Outside() == Outside(0f0)
    @test Auto(2).ratio == 2.0
end

@testset "gridlayout constructor colsizes" begin
    gl = GridLayout(2, 2; colsizes = Fixed(10), rowsizes = Relative(0.5))

    @test gl.colsizes == GridLayoutBase.ContentSize[Fixed(10), Fixed(10)]
    @test gl.rowsizes == GridLayoutBase.ContentSize[Relative(0.5), Relative(0.5)]

    gl2 = GridLayout(2, 2;
        colsizes = GridLayoutBase.ContentSize[Fixed(10), Relative(0.3)],
        rowsizes = GridLayoutBase.ContentSize[Auto(false), Auto(true)])
    @test gl2.colsizes == GridLayoutBase.ContentSize[Fixed(10), Relative(0.3)]
    @test gl2.rowsizes == GridLayoutBase.ContentSize[Auto(false), Auto(true)]

    @test_throws ErrorException GridLayout(; colsizes = "abc")
    @test_throws ErrorException GridLayout(; rowsizes = missing)

    gl3 = GridLayout(3, 3; addedcolgaps = Fixed(20), addedrowgaps = Fixed(0))
    @test gl3.addedcolgaps == GridLayoutBase.GapSize[Fixed(20), Fixed(20)]
    @test gl3.addedrowgaps == GridLayoutBase.GapSize[Fixed(0), Fixed(0)]

    @test_throws ErrorException GridLayout(3, 1; addedcolgaps = "abc")
    @test_throws ErrorException GridLayout(1, 3; addedrowgaps = "abc")

    @test_throws ErrorException GridLayout(3, 1; addedcolgaps = [Fixed(20)])
    @test_throws ErrorException GridLayout(3, 1; addedrowgaps = [Fixed(30)])

    @test_throws ErrorException GridLayout(3, 1; rowsizes = [Fixed(20)])
    @test_throws ErrorException GridLayout(3, 2; colsizes = [Relative(0.4)])

    @test_throws ErrorException GridLayout(0, 1)
    @test_throws ErrorException GridLayout(1, 0)
end


@testset "printing gridlayouts" begin
    gl = GridLayout(3, 3)
    gl[1, 1] = DebugRect()
    gl[2:3, 4:5] = DebugRect()

    text_long = repr(MIME"text/plain"(), gl)
    @test text_long == """
    GridLayout[1:3, 1:5] with 2 children
     ┣━ [1, 1] DebugRect
     ┗━ [2:3, 4:5] DebugRect
    """
    text_short = repr(gl)
    @test text_short == "GridLayout[3, 5] (2 children)"

    subgl = gl[1, 2] = GridLayout()
    subgl[1:5, 3] = DebugRect()

    text_longer = repr(MIME"text/plain"(), gl)
    # this is actually a bit buggy with the newline space space newline at the end
    @test text_longer == "GridLayout[1:3, 1:5] with 3 children\n ┣━ [1, 1] DebugRect\n ┣━ [2:3, 4:5] DebugRect\n ┗━ [1, 2] GridLayout[1:5, 1:3] with 1 children\n   ┗━ [1:5, 3] DebugRect\n  \n"


    gl3 = GridLayout()
    gl4 = gl3[1, 1] = GridLayout()
    gl4[1, 1] = DebugRect()
    gl3[2, 2] = DebugRect()
    text_long_downconnection = repr(MIME"text/plain"(), gl3)

    # this is also a bit buggy for the same reason as above
    @test text_long_downconnection == "GridLayout[1:2, 1:2] with 2 children\n ┣━ [1, 1] GridLayout[1:1, 1:1] with 1 children\n ┃ ┗━ [1, 1] DebugRect\n ┃\n ┗━ [2, 2] DebugRect\n"
end

@testset "vector and array assigning" begin
    gl = GridLayout()
    gl[1, 1:3] = [DebugRect() for i in 1:3]
    @test size(gl) == (1, 3)

    gl2 = GridLayout()
    gl2[2:3, 2] = [DebugRect() for i in 1:2]
    @test size(gl2) == (3, 2)

    gl3 = GridLayout()
    gl3[1:3, 1:4] = [DebugRect() for i in 1:12]
    @test size(gl3) == (3, 4)

    gl4 = GridLayout()
    gl4[1:3, 1:4] = [DebugRect() for i in 1:3, j in 1:4]
    @test size(gl4) == (3, 4)

    @test_throws ErrorException gl[1, 1:3] = [DebugRect() for i in 1:2]
    @test_throws ErrorException gl[1:3, 2] = [DebugRect() for i in 1:2]
    @test_throws ErrorException gl[1:3, 1:3] = [DebugRect() for i in 1:10]
    @test_throws ErrorException gl[1:3, 1:3] = [DebugRect() for i in 1:3, j in 1:4]

    gl5 = GridLayout()
    gl5[] = [DebugRect() for i in 1:2, j in 1:3]
    @test size(gl5) == (2, 3)

    gl6 = GridLayout()
    gl6[:v] = [DebugRect() for i in 1:3]
    @test size(gl6) == (3, 1)

    gl7 = GridLayout()
    gl7[:h] = [DebugRect() for i in 1:3]
    @test size(gl7) == (1, 3)

    @test_throws ErrorException gl7[:abc] = [DebugRect() for i in 1:3]
    @test_throws ErrorException gl7[] = [DebugRect() for i in 1:3]
end

@testset "grid api" begin

    gl1 = grid!([1:2, 1:2] => DebugRect(), [3, :] => DebugRect())
    @test size(gl1) == (3, 2)
    @test gl1.content[1].span == GridLayoutBase.Span(1:2, 1:2)
    @test gl1.content[2].span == GridLayoutBase.Span(3:3, 1:2)

    gl2 = grid!([DebugRect() for i in 1:3, j in 1:2])
    @test size(gl2) == (3, 2)
    for i in 1:3, j in 1:2
        n = (i - 1) * 2 + j
        @test gl2.content[n].span == GridLayoutBase.Span(i:i, j:j)
    end

    gl3 = vbox!(DebugRect(), DebugRect())
    @test size(gl3) == (2, 1)
    for i in 1:2
        @test gl3.content[i].span == GridLayoutBase.Span(i:i, 1:1)
    end

    gl4 = hbox!(DebugRect(), DebugRect())
    @test size(gl4) == (1, 2)
    for i in 1:2
        @test gl4.content[i].span == GridLayoutBase.Span(1:1, i:i)
    end
end

@testset "gridnest" begin
    layout = GridLayout()
    dr = layout[1:2, 3:4] = DebugRect()
    subgl = gridnest!(layout, 1:2, 3:4)

    @test size(subgl) == (2, 2)
    @test subgl.content[1].span == GridLayoutBase.Span(1:2, 1:2)
end

@testset "invalid removal" begin
    layout = GridLayout()
    dr = layout[1, 1] = DebugRect()
    # remove the item outside of the normal path
    deleteat!(layout.content, 1)
    # place the item somewhere else, this should error now
    @test_throws ErrorException layout[1, 2] = dr
end

@testset "equal protrusion gaps" begin
    bbox = BBox(0, 1000, 0, 700)
    layout = GridLayout(bbox = bbox, alignmode = Outside(0))
    subgl = layout[1, 1] = GridLayout(3, 3, equalprotrusiongaps = (true, true),
        addedcolgaps = Fixed(0), addedrowgaps = Fixed(0))
    subgl[1, 1, BottomRight()] = DebugRect(width = 100, height = 100)

    dr1 = subgl[1, 1] = DebugRect()
    dr2 = subgl[2, 2] = DebugRect()
    dr3 = subgl[3, 3] = DebugRect()

    @test width(computedbboxobservable(dr1)[]) ≈ (1000 - 2 * 100) / 3.0f0
    @test width(computedbboxobservable(dr2)[]) ≈ (1000 - 2 * 100) / 3.0f0
    @test width(computedbboxobservable(dr3)[]) ≈ (1000 - 2 * 100) / 3.0f0
    @test height(computedbboxobservable(dr1)[]) ≈ (700 - 2 * 100) / 3.0f0
    @test height(computedbboxobservable(dr2)[]) ≈ (700 - 2 * 100) / 3.0f0
    @test height(computedbboxobservable(dr3)[]) ≈ (700 - 2 * 100) / 3.0f0
end

@testset "getindex gridposition" begin
    layout = GridLayout()
    dr = layout[2, 2] = DebugRect()

    @test layout[2, 2] == GridPosition(layout, 2, 2)
    @test layout[end, end] == GridPosition(layout, 2, 2)
    @test_throws ErrorException layout[2, 2, end]

    gp = GridPosition(layout, 1, 1)
    gp[] = dr
    @test gridcontent(dr).span == GridLayoutBase.Span(1:1, 1:1)
end

@testset "gridposition contents" begin
    layout = GridLayout()
    dr1 = layout[1, 1] = DebugRect()
    dr2 = layout[1, 2] = DebugRect()
    dr3 = layout[2, 1] = DebugRect()
    dr4 = layout[2, 2] = DebugRect()

    co = contents(layout[1, 1])
    @test co == [dr1]

    co = contents(layout[1:2, 1])
    @test co == [dr1, dr3]

    co = contents(layout[2, 1:2])
    @test co == [dr3, dr4]

    co = contents(layout[:, :])
    @test co == [dr1, dr2, dr3, dr4]

    dr5 = layout[2, 2, Right()] = DebugRect()

    co = contents(layout[:, :]) # implicit Inner() side
    @test co == [dr1, dr2, dr3, dr4]

    @test contents(layout[2, 2, Right()]) == [dr5]

    dr6 = layout[1:2, 2] = DebugRect()
    @test contents(layout[1:2, 2]) == [dr2, dr4, dr6]
    @test contents(layout[1:2, 2], exact = true) == [dr6]

    @test contents(layout) == [dr1, dr2, dr3, dr4, dr5, dr6]
end

@testset "span containment" begin
    Span = GridLayoutBase.Span
    @test Span(1:2, 2:3) in Span(1:2, 2:3)
    @test Span(1:2, 2:3) in Span(0:3, 1:4)
    @test !(Span(0:3, 1:4) in Span(1:2, 2:3))
    @test !(Span(1:2, 1:4) in Span(1:2, 2:3))
    @test !(Span(0:3, 2:3) in Span(1:2, 2:3))
end

@testset "layoutobservables undefined" begin
    struct MyType end
    @test_throws ErrorException GridLayoutBase.layoutobservables(MyType())
end

@testset "determine gridlayout size" begin
    gl = GridLayout(;alignmode = Outside(0))
    gl[1, 1] = DebugRect(width = 800, height = 600)

    @test GridLayoutBase.determinedirsize(gl, GridLayoutBase.Row()) == 600
    @test GridLayoutBase.determinedirsize(gl, GridLayoutBase.Col()) == 800

    gl[1, 1, Left()] = DebugRect(width = 200)
    @test GridLayoutBase.determinedirsize(gl, GridLayoutBase.Col()) == 800 + 200
    gl[1, 1, Right()] = DebugRect(width = 200)
    @test GridLayoutBase.determinedirsize(gl, GridLayoutBase.Col()) == 800 + 200 + 200

    gl[1, 1, Top()] = DebugRect(height = 100)
    @test GridLayoutBase.determinedirsize(gl, GridLayoutBase.Row()) == 600 + 100
    gl[1, 1, Bottom()] = DebugRect(height = 100)
    @test GridLayoutBase.determinedirsize(gl, GridLayoutBase.Row()) == 600 + 100 + 100

    gl.alignmode[] = Outside(50)
    @test GridLayoutBase.determinedirsize(gl, GridLayoutBase.Col()) == 800 + 200 + 200 + 2 * 50
    @test GridLayoutBase.determinedirsize(gl, GridLayoutBase.Row()) == 600 + 100 + 100 + 2 * 50

    gl.alignmode[] = Inside()
    @test GridLayoutBase.determinedirsize(gl, GridLayoutBase.Col()) == 800
    @test GridLayoutBase.determinedirsize(gl, GridLayoutBase.Row()) == 600
end

@testset "alignment" begin
    boxwidth = 1000
    boxheight = 800

    gl = GridLayout(; bbox = (0, boxwidth, 0, boxheight), alignmode = Outside(0))
    gl[1, 1] = DebugRect(width = boxwidth, height = boxheight)
    halign = Observable{Any}(:center)
    valign = Observable{Any}(:center)

    rectwidth = 400
    rectheight = 200

    dr = gl[1, 1] = DebugRect(width = rectwidth, height = 200, halign = halign, valign = valign)

    @test computedbboxobservable(dr)[] == BBox(boxwidth/2 - rectwidth/2, boxwidth/2 + rectwidth/2, boxheight/2 - rectheight/2, boxheight/2 + rectheight/2)
    halign[] = :left
    @test computedbboxobservable(dr)[] == BBox(0, rectwidth, boxheight/2 - rectheight/2, boxheight/2 + rectheight/2)
    halign[] = :right
    @test computedbboxobservable(dr)[] == BBox(boxwidth - rectwidth, boxwidth, boxheight/2 - rectheight/2, boxheight/2 + rectheight/2)
    valign[] = :top
    @test computedbboxobservable(dr)[] == BBox(boxwidth - rectwidth, boxwidth, boxheight - rectheight, boxheight)
    valign[] = :bottom
    @test computedbboxobservable(dr)[] == BBox(boxwidth - rectwidth, boxwidth, 0, rectheight)

    valign[] = 0.5
    @test computedbboxobservable(dr)[] == BBox(boxwidth - rectwidth, boxwidth,  boxheight/2 - rectheight/2, boxheight/2 + rectheight/2)
    halign[] = 0.5
    @test computedbboxobservable(dr)[] == BBox(boxwidth/2 - rectwidth/2, boxwidth/2 + rectwidth/2,  boxheight/2 - rectheight/2, boxheight/2 + rectheight/2)

    @test_throws ErrorException valign[] = :abc
    valign[] = 0.5
    @test_throws ErrorException halign[] = :abc
end


@testset "invalid gridlayout construction" begin
    @test_throws ErrorException GridLayout(2, 2; rowsizes = [Fixed(10)])
    @test_throws ErrorException GridLayout(2, 2; colsizes = [Fixed(10)])
    @test_throws ErrorException GridLayout(2, 2; addedrowgaps = [Fixed(10), Fixed(10)])
    @test_throws ErrorException GridLayout(2, 2; addedcolgaps = [Fixed(10), Fixed(10)])
    @test_throws ErrorException GridLayout(2, 2; rowsizes = [Fixed(10), Fixed(10)], addedrowgaps = [Fixed(10), Fixed(10)])
    @test_throws ErrorException GridLayout(2, 2; colsizes = [Fixed(10), Fixed(10)], addedcolgaps = [Fixed(10), Fixed(10)])
end

@testset "tellsize" begin
    dr = DebugRect(width = 100, height = 200, tellwidth = true, tellheight = false)
    @test reporteddimensionsobservable(dr)[].inner == (100, nothing)
    dr = DebugRect(width = 100, height = 200, tellwidth = false, tellheight = true)
    @test reporteddimensionsobservable(dr)[].inner == (nothing, 200)
end


@testset "aspect sizes referencing auto sizes" begin
    bbox = BBox(0, 1000, 0, 1000)
    layout = GridLayout(2, 2, bbox = bbox, alignmode = Outside(0))

    colgap!(layout, 0)
    rowgap!(layout, 0)

    dr = layout[1, 1] = DebugRect()

    colsize!(layout, 1, Aspect(1, 1.5))
    @test computedbboxobservable(dr)[] == BBox(0, 750, 500, 1000)

    colsize!(layout, 1, Auto())

    rowsize!(layout, 1, Aspect(1, 1.5))
    @test computedbboxobservable(dr)[] == BBox(0, 500, 250, 1000)

    @test_throws ErrorException colsize!(layout, 1, Aspect(1, 1))
end

@testset "integer rect2 suggestedbbox" begin
    layout = GridLayout(bbox = Observable(GridLayoutBase.GeometryBasics.Rect(0, 10, 20, 30)))
    @test suggestedbboxobservable(layout)[] == GridLayoutBase.GeometryBasics.Rect2f(0, 10, 20, 30)
end

@testset "GridSubpositions" begin
    l = GridLayout()
    gp = l[1, 1]
    @test gp isa GridPosition
    gsp = gp[1, 1]
    @test gsp isa GridSubposition
    @test isempty(contents(gp))
    r = gsp[] = DebugRect()
    @test content(gp) isa GridLayout
    @test content(gsp) == r

    r2 = l[1, 1][1, 2] = DebugRect()
    @test content(l[1, 1][1, 2]) == r2
    r3 = l[1, 1][1, 3][1, 1] = DebugRect()
    @test content(l[1, 1][1, 3][1, 1]) == r3
end

@testset "Parents" begin
    g = GridLayout()
    @test GridLayoutBase.parent(g) === nothing
    @test GridLayoutBase.top_parent(g) === nothing
    g2 = GridLayout()
    @test GridLayoutBase.parent(g2) === nothing
    g[1, 1] = g2
    @test GridLayoutBase.parent(g2) === g
    @test GridLayoutBase.top_parent(g2) === nothing
    @test GridLayoutBase.top_parent_grid(g2) === g

    g3 = GridLayout()
    @test GridLayoutBase.parent(g3) === nothing
    g2[1, 1] = g3
    @test GridLayoutBase.parent(g3) === g2
    @test GridLayoutBase.top_parent(g3) === nothing
    @test GridLayoutBase.top_parent_grid(g3) === g

    # make some arbitrary parent object with an id
    p = Dict()
    g.parent = p
    @test GridLayoutBase.parent(g) === p
    @test GridLayoutBase.top_parent(g) === p
    @test GridLayoutBase.top_parent(g2) === p
    @test GridLayoutBase.top_parent(g3) === p
    @test GridLayoutBase.top_parent_grid(g2) === g
    @test GridLayoutBase.top_parent_grid(g3) === g
end

@testset "Number of updates" begin
    g = GridLayout()
    n = Ref(0)
    on(g.layoutobservables.suggestedbbox) do _
        n[] += 1
    end

    g2 = GridLayout()
    g[1, 1] = g2
    @test n[] == 0

    m = Ref(0)
    on(g2.layoutobservables.suggestedbbox) do _
        m[] += 1
    end
    for i in 1:10
        g2[1, i] = GridLayout()
    end
    # one update for each gridlayout
    @test m[] == 10
    # g shouldn't have changed
    @test n[] == 0
    with_updates_suspended(g2) do
        for i in 1:10
            g2[1, i] = GridLayout()
        end
    end
    # only one update should have happened at the end
    @test m[] == 11
    # still nothing for g
    @test n[] == 0
end

@testset "Number of updates 2" begin
    gl = GridLayout()
    gl2 = GridLayout(gl[1, 1])
    dr = gl2[1, 1] = DebugRect(width = 100, height = 100)
    dr2 = gl2[1, 2] = DebugRect(width =  100, height = 200)
    n = Ref(0)
    on(dr2.layoutobservables.suggestedbbox) do bb
        n[] += 1
    end
    dr2.rightprot[] = 10
    @test n[] == 1
end

@testset "GridLayout alignment when solved size mismatches suggested bbox" begin
    let
        bbox = BBox(0, 2000, 0, 1000)
        outer = GridLayout(bbox = bbox, alignmode = Outside(0), halign = :center)
        dr1 = outer[1, 1] = DebugRect()
        dr2 = outer[1, 2] = DebugRect()
        colsize!(outer, 2, 500)
        colsize!(outer, 1, Aspect(1, 1))
        colgap!(outer, 0)
        @test computedbboxobservable(dr1)[] == BBox(250, 1250, 0, 1000)
        @test computedbboxobservable(dr2)[] == BBox(1250, 1750, 0, 1000)

        outer.halign = :left
        @test computedbboxobservable(dr1)[] == BBox(0, 1000, 0, 1000)
        @test computedbboxobservable(dr2)[] == BBox(1000, 1500, 0, 1000)

        outer.halign = :right
        @test computedbboxobservable(dr1)[] == BBox(500, 1500, 0, 1000)
        @test computedbboxobservable(dr2)[] == BBox(1500, 2000, 0, 1000)
    end
    let
        bbox = BBox(0, 1000, 0, 2000)
        outer = GridLayout(bbox = bbox, alignmode = Outside(0), valign = :center)
        dr1 = outer[1, 1] = DebugRect()
        dr2 = outer[2, 1] = DebugRect()
        rowsize!(outer, 2, 500)
        rowsize!(outer, 1, Aspect(1, 1))
        rowgap!(outer, 0)
        @test computedbboxobservable(dr1)[] == BBox(0, 1000, 750, 1750)
        @test computedbboxobservable(dr2)[] == BBox(0, 1000, 250, 750)

        outer.valign = :top
        @test computedbboxobservable(dr1)[] == BBox(0, 1000, 1000, 2000)
        @test computedbboxobservable(dr2)[] == BBox(0, 1000, 500, 1000)

        outer.valign = :bottom
        @test computedbboxobservable(dr1)[] == BBox(0, 1000, 500, 1500)
        @test computedbboxobservable(dr2)[] == BBox(0, 1000, 0, 500)
    end
end

@testset "tight_bbox" begin
    bbox = BBox(0, 2000, 0, 1000)

    gl = GridLayout(bbox = bbox, alignmode = Outside(0), halign = :center)
    dr1 = gl[1, 1] = DebugRect()
    dr2 = gl[1, 2] = DebugRect()
    colsize!(gl, 2, 500)
    colsize!(gl, 1, Aspect(1, 1))
    colgap!(gl, 0)

    @test tight_bbox(gl) == BBox(250, 1750, 0, 1000)
    gl.alignmode[] = Outside(0, 100, 0, 200)
    # inner height now 800, so colwidth 1 also 800, plus 500 = 1300, plus padding 1400
    # but shifted because of unequal padding
    @test tight_bbox(gl) == BBox(300, 1700, 0, 1000)

    dr1.topprot[] = 100
    dr1.leftprot[] = 50
    # inner height now 700, so colwidth 1 also 700, plus 500 = 1200, plus padding + leftprot = 1350
    @test tight_bbox(gl) == BBox(325, 1675, 0, 1000)

    gl.alignmode[] = Inside()
    # padding and protrusions are gone, so height is again 1000
    @test tight_bbox(gl) == BBox(250, 1750, 0, 1000)

    gl.alignmode[] = Mixed()
    # not implemented yet
    @test_throws ErrorException tight_bbox(gl)
end

@testset "offset row/col auto size" begin
    bbox = BBox(0, 1000, 0, 1000)
    gl = GridLayout(bbox = bbox, alignmode = Outside(0), halign = :center)
    dr1 = gl[0, 1] = DebugRect(height = 200)
    dr2 = gl[1, 0] = DebugRect(width = 100)
    dr3 = gl[1, 1] = DebugRect()
    colgap!(gl, 0)
    rowgap!(gl, 0)

    @test GridLayoutBase.determinedirsize(0, gl, GridLayoutBase.Col()) == 100
    @test GridLayoutBase.determinedirsize(0, gl, GridLayoutBase.Row()) == 200

    @test suggestedbboxobservable(dr1)[] == BBox(100, 1000, 800, 1000)
    @test suggestedbboxobservable(dr2)[] == BBox(0, 100, 0, 800)
    @test suggestedbboxobservable(dr3)[] == BBox(100, 1000, 0, 800)
end

@testset "offset aspect" begin
    bbox = BBox(0, 1000, 0, 1000)
    gl = GridLayout(bbox = bbox, alignmode = Outside(0), halign = :center)
    # make columns and rows offset
    gl[0, 0] = DebugRect()
    colgap!(gl, 0)
    rowgap!(gl, 0)
    # aspect refers to row 0 
    colsize!(gl, 1, Aspect(0, 1.5))

    maxgrid, gridboxes = GridLayoutBase.compute_rowcols(gl, suggestedbboxobservable(gl)[])
    @test gridboxes.lefts == [0, 250]
    @test gridboxes.rights == [250, 1000]
    @test gridboxes.tops == [1000, 500]
    @test gridboxes.bottoms == [500, 0]

    bbox = BBox(0, 1000, 0, 1000)
    gl = GridLayout(bbox = bbox, alignmode = Outside(0), halign = :center)
    # make columns and rows offset
    gl[0, 0] = DebugRect()
    colgap!(gl, 0)
    rowgap!(gl, 0)
    # aspect refers to column 0
    rowsize!(gl, 1, Aspect(0, 1.5))

    maxgrid, gridboxes = GridLayoutBase.compute_rowcols(gl, suggestedbboxobservable(gl)[])
    @test gridboxes.lefts == [0, 500]
    @test gridboxes.rights == [500, 1000]
    @test gridboxes.tops == [1000, 750]
    @test gridboxes.bottoms == [750, 0]
end

# issue 37
@testset "Aspect indexing bug" begin
    gl = GridLayout()
    gl[0, 1] = DebugRect(height = 200)
    @test_nowarn colsize!(gl, 1, Aspect(0, 1.0))

    gl = GridLayout()
    gl[1, 0] = DebugRect(width = 200)
    @test_nowarn rowsize!(gl, 1, Aspect(0, 1.0))
end

@testset "colon rows/cols with offsets" begin
    gl = GridLayout()
    gl[0, 0] = GridLayout()
    g1 = gl[-1, :] = GridLayout()
    g2 = gl[:, -1] = GridLayout()
    @test gridcontent(g1).span == GridLayoutBase.Span(-1:-1, 0:1)
    @test gridcontent(g2).span == GridLayoutBase.Span(-1:1, -1:-1)
end

# https://github.com/JuliaPlots/Makie.jl/issues/2018
@testset "fixed size dirgaps" begin
    gl = GridLayout(width = 800, height = 600)
    gl[1, 1] = DebugRect(width = 600, height = 600)
    @test begin
        gl[0, 1] = DebugRect(height = 40, width = 100)
        true
    end
end
@testset "Effective protrusions" begin
    bbox = BBox(0, 1000, 0, 1000)
    gl = GridLayout(bbox = bbox, alignmode = Outside(0))

    dr1 = gl[1, 1] = DebugRect(rightprot = 50, bottomprot = 100, alignmode = Inside())
    @test GridLayoutBase.protrusionsobservable(dr1)[] == RectSides{Float32}(0, 50, 100, 0)
    @test GridLayoutBase.suggestedbboxobservable(dr1)[] == BBox(0, 950, 100, 1000)

    dr1.alignmode[] = Mixed(right = Protrusion(0), bottom = Protrusion(0))
    @test GridLayoutBase.protrusionsobservable(dr1)[] == RectSides{Float32}(0, 50, 100, 0)
    @test reporteddimensionsobservable(dr1)[].outer == RectSides{Float32}(0, 0, 0, 0)
    @test GridLayoutBase.suggestedbboxobservable(dr1)[] == BBox(0, 1000, 0, 1000)
    @test GridLayoutBase.computedbboxobservable(dr1)[] == BBox(0, 1000, 0, 1000)

    dr1.alignmode[] = Mixed(right = 50, bottom = 50)
    @test GridLayoutBase.protrusionsobservable(dr1)[] == RectSides{Float32}(0, 50, 100, 0)
    @test reporteddimensionsobservable(dr1)[].outer == RectSides{Float32}(0, 0, 0, 0)
    @test GridLayoutBase.suggestedbboxobservable(dr1)[] == BBox(0, 1000, 0, 1000)
    @test GridLayoutBase.computedbboxobservable(dr1)[] == BBox(0, 900, 150, 1000)

    dr1.alignmode[] = Mixed(right = Protrusion(50), bottom = Protrusion(50), left = Protrusion(50), top = Protrusion(50))
    @test GridLayoutBase.protrusionsobservable(dr1)[] == RectSides{Float32}(0, 50, 100, 0)
    @test reporteddimensionsobservable(dr1)[].outer == RectSides{Float32}(50, 50, 50, 50)
    @test GridLayoutBase.suggestedbboxobservable(dr1)[] == BBox(50, 950, 50, 950)
    @test GridLayoutBase.computedbboxobservable(dr1)[] == BBox(50, 950, 50, 950)

    dr1.layoutobservables.protrusions[] = RectSides{Float32}(0, 0, 0, 0)
    @test GridLayoutBase.protrusionsobservable(dr1)[] == RectSides{Float32}(0, 0, 0, 0)
    @test reporteddimensionsobservable(dr1)[].outer == RectSides{Float32}(50, 50, 50, 50)

    dr1.alignmode[] = Mixed(left = Protrusion(1), right = Protrusion(2), bottom = Protrusion(3), top = Protrusion(4))
    @test GridLayoutBase.protrusionsobservable(dr1)[] == RectSides{Float32}(0, 0, 0, 0)
    @test reporteddimensionsobservable(dr1)[].outer == RectSides{Float32}(1, 2, 3, 4)

    gl = GridLayout()
    @test GridLayoutBase.protrusionsobservable(gl)[] == RectSides{Float32}(0, 0, 0, 0)
    gl[1, 1, Right()] = DebugRect(width = 200)
    @test GridLayoutBase.protrusionsobservable(gl)[] == RectSides{Float32}(0, 200, 0, 0)
    @test reporteddimensionsobservable(gl)[].outer == RectSides{Float32}(0, 200, 0, 0)
end

@testset "GridLayout GridPosition/GridSubposition constructor" begin
    gl = GridLayout()
    gl2 = GridLayout(gl[1, 1])
    @test gl2.parent == gl
    @test gl2.layoutobservables.gridcontent[].span == GridLayoutBase.Span(1:1, 1:1)
    gl3 = GridLayout(gl[1, 2][1, 1])
    @test gl3.parent != gl
    @test gl3.parent == gl.content[2].content
    @test gl3.layoutobservables.gridcontent[].span == GridLayoutBase.Span(1:1, 1:1)
    @test gl.content[2].content.layoutobservables.gridcontent[].span == GridLayoutBase.Span(1:1, 2:2)
end

@testset "mixed alignmode bug" begin
    gl = GridLayout(bbox = BBox(0, 1000, 0, 500), alignmode = Outside(0))
    gl2 = GridLayout(gl[1, 1])
    dr1 = gl2[1, 1] = DebugRect()
    dr2 = gl2[1, 2] = DebugRect()
    dr3 = gl2[2, 1] = DebugRect()
    dr4 = gl2[2, 2] = DebugRect(width = 200)
    colgap!(gl2, 0)
    rowgap!(gl2, 0)
    @test reporteddimensionsobservable(dr4)[].inner == (200, nothing)
    @test computedbboxobservable(dr4)[] == BBox(800, 1000, 0, 250)
    @test suggestedbboxobservable(dr4)[] == BBox(800, 1000, 0, 250)
    @test reporteddimensionsobservable(dr4)[].outer == RectSides(0f0, 0f0, 0f0, 0f0)
    dr4.rightprot[] = 100
    @test reporteddimensionsobservable(dr4)[].inner == (200, nothing)
    @test reporteddimensionsobservable(dr4)[].outer == RectSides(0f0, 100f0, 0f0, 0f0)
    @test suggestedbboxobservable(dr4)[] == BBox(700, 900, 0, 250)
    @test computedbboxobservable(dr4)[] == BBox(700, 900, 0, 250)
    dr4.alignmode[] = Outside(0)
    @test reporteddimensionsobservable(dr4)[].inner == (300, nothing)
    @test reporteddimensionsobservable(dr4)[].outer == RectSides(0f0, 0f0, 0f0, 0f0)
    @test suggestedbboxobservable(dr4)[] == BBox(700, 1000, 0, 250)
    @test computedbboxobservable(dr4)[] == BBox(700, 900, 0, 250)
end
