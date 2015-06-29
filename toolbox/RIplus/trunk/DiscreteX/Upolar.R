Upolar <- function(rt1, rt2) {
    # just ||x - y]]^2 in polar coordinates
    -.5 * ( rt1[1]^2 + rt2[1]^2 - 2 * rt1[1] * rt2[1] * cos(rt1[2] - rt2[2]) )
}
