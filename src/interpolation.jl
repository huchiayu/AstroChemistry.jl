using Parameters
@with_kw struct InterpRange{T}
    xmin::T
    xmax::T
    dx_inv::T
    ymin::T
    ymax::T
    dy_inv::T
end

function interp(x::T,y::T, tbl::Matrix{T}, par::InterpRange{T}) where {T}
    @unpack xmin,xmax,dx_inv,ymin,ymax,dy_inv = par
    x = x >= xmax ? (1-1e-7)*xmax : x
    x = x <= xmin ? xmin : x
    y = y >= ymax ? (1-1e-7)*ymax : y
    y = y <= ymin ? ymin : y
    xi = (x - xmin) * dx_inv
    x0 = Int64(floor(xi))
    yi = (y - ymin) * dy_inv
    y0 = Int64(floor(yi))
    wx1 = xi - x0
    wy1 = yi - y0
    wx0 = 1.0 - wx1
    wy0 = 1.0 - wy1
    f00 = tbl[x0+1,y0+1] #julia array is 1-based
    f10 = tbl[x0+2,y0+1]
    f01 = tbl[x0+1,y0+2]
    f11 = tbl[x0+2,y0+2]

    fi = wx0*wy0 * f00 + wx0*wy1 * f01 + wx1*wy0 * f10 + wx1*wy1 * f11
end
