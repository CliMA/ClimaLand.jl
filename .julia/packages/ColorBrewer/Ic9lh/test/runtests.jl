using ColorBrewer
using Test

### Test 1: ArgumentError thrown when the scheme is invalid.
@test_throws ArgumentError palette("NOTASCHEME", 5)

### Test 2: ArgumentError thrown when the number of categories for a scheme
#   is invalid.
@test_throws ArgumentError palette("Greens", 2)

### Test 3: Function runs correctly.
palette("Pastel1", 6)

### Test 4: Validate 9-color YlOrRd palette
palette("YlOrRd", 9)