function validate_poly(p)
    if first(p) != last(p)
        error("Polygon should have first and last elements equal")
    end
end
