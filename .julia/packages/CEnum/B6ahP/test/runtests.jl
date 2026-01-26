using CEnum
using Test

@cenum(ZeroOne, zero, one)
@test zero == 0
@test one == 1

@cenum(Fruit, apple=1, orange=2, kiwi=2)
@test orange == kiwi
@test apple | orange == 3
@test apple & orange == 0
@test apple + 1 == 2
@test kiwi - 1 == 1
@test kiwi ⊻ kiwi == 0
@test ~orange == -3
@test_nowarn print(devnull, Fruit(apple | orange))

@cenum(Boolean::Bool, alternativefact, fact)
@test alternativefact == false

@cenum(Day, Mon=1, Tue, Wed=3, suiyoubi=3, Fri=5, Sat)
@test Mon == 1
@test Tue == 2
@test Wed == 3 == suiyoubi
@test Fri == 5
@test Sat == 6

# issue #11
@test 0 - one == -1

# These tests are derived from https://en.cppreference.com/w/c/language/operator_arithmetic section "Shift operators"
# Here we test to ensure our promotion behavior matches that of C
@cenum var"##SHIFT_TEST_SIGNED#1"::Int64 begin
    SHIFT_TEST_N1000 = -1000
end
@cenum var"##SHIFT_TEST_UNSIGNED#1"::UInt64 begin
    SHIFT_TEST_123 = 0x123
    SHIFT_TEST_10 = 0x10
end
@test SHIFT_TEST_123 << 1 == 0x246
@test SHIFT_TEST_123 << 63 == 0x8000000000000000
@test SHIFT_TEST_10 << 10 == 0x4000
@test SHIFT_TEST_N1000 >> 1 == -500
