export is_commutative_semiring

# this function is used for testing
function is_commutative_semiring(a::T, b::T, c::T) where T
    # +
    (a + b) + c == a + (b + c) &&
    a + zero(T) == zero(T) + a == a &&
    a + b == b + a &&
    # *
    (a * b) * c == a * (b * c) &&
    a * one(T) == one(T) * a == a &&
    a * b == b * a &&
    # more
    a * (b+c) == a*b + a*c &&
    (a+b) * c == a*c + b*c &&
    a * zero(T) == zero(T) * a == zero(T)
end