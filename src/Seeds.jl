"""
    uniform_seed1D(Nx; direction=(0.0, 0.0, 1.0))

Create a uniformly magnetized 1D seed with all spins aligned along `direction`.
The direction is normalized before being written into the returned `(3, Nx)`
array.
"""
function uniform_seed1D(Nx::Integer; direction = (0.0, 0.0, 1.0))
    Nx > 0 || throw(ArgumentError("Nx must be positive."))
    ux, uy, uz = _normalized_direction(direction)

    spins = Array{Float64}(undef, 3, Nx)
    @inbounds for i in 1:Nx
        spins[1, i] = ux
        spins[2, i] = uy
        spins[3, i] = uz
    end

    return spins
end

"""
    uniform_seed2D(Nx, Ny; direction=(0.0, 0.0, 1.0))

Create a uniformly magnetized 2D seed with all spins aligned along `direction`.
The direction is normalized before being written into the returned `(3, Nx, Ny)`
array.
"""
function uniform_seed2D(Nx::Integer, Ny::Integer; direction = (0.0, 0.0, 1.0))
    Nx > 0 || throw(ArgumentError("Nx must be positive."))
    Ny > 0 || throw(ArgumentError("Ny must be positive."))
    ux, uy, uz = _normalized_direction(direction)

    spins = Array{Float64}(undef, 3, Nx, Ny)
    @inbounds for j in 1:Ny, i in 1:Nx
        spins[1, i, j] = ux
        spins[2, i, j] = uy
        spins[3, i, j] = uz
    end

    return spins
end

"""
    paint_domain_wall1D!(
        spins;
        center=(size(spins, 2) + 1) / 2,
        width=5.0,
        domain_direction=(0.0, 0.0, 1.0),
        wall_normal=(1.0, 0.0, 0.0),
        background_direction=nothing,
    )

Paint a 1D domain wall in-place on top of a uniformly magnetized background.

The two domains point along `±domain_direction`, while `wall_normal` sets the
magnetization direction at the wall center after projection perpendicular to the
domain direction. For example, with `domain_direction = (0, 0, 1)`,
`wall_normal = (1, 0, 0)` gives a Neel wall and `wall_normal = (0, 1, 0)` gives
a Bloch wall.

This modifier is additive, so it is intended to be applied to a mostly uniform
state, possibly multiple times, followed by [`normalize_spins!`](@ref). The
painted wall profile is corrected by subtracting `background_direction`, which
defaults to `domain_direction`. When painting multiple walls on the same seed,
set `background_direction` to the original uniform seed direction.
"""
function paint_domain_wall1D!(
    spins::AbstractArray{<:AbstractFloat, 2};
    center::Real = (size(spins, 2) + 1) / 2,
    width::Real = 5.0,
    domain_direction = (0.0, 0.0, 1.0),
    wall_normal = (1.0, 0.0, 0.0),
    background_direction = nothing,
)
    size(spins, 1) == 3 || throw(ArgumentError("Expected a spin array of size (3, Nx)."))
    width > 0 || throw(ArgumentError("width must be positive."))

    dx, dy, dz = _normalized_direction(domain_direction)
    bx, by, bz = isnothing(background_direction) ? (dx, dy, dz) : _normalized_direction(background_direction)
    wx, wy, wz = _transverse_direction(wall_normal, (dx, dy, dz))

    @inbounds for i in axes(spins, 2)
        ξ = (i - center) / width
        longitudinal = tanh(ξ)
        transverse = sech(ξ)

        spins[1, i] += longitudinal * dx + transverse * wx - bx
        spins[2, i] += longitudinal * dy + transverse * wy - by
        spins[3, i] += longitudinal * dz + transverse * wz - bz
    end

    return spins
end

"""
    paint_domain_wall2D!(
        spins;
        point=((size(spins, 2) + 1) / 2, (size(spins, 3) + 1) / 2),
        slope=0.0,
        width=5.0,
        domain_direction=(0.0, 0.0, 1.0),
        wall_normal=(1.0, 0.0, 0.0),
        background_direction=nothing,
    )

Paint a 2D domain wall in-place on top of a uniformly magnetized background.

The wall front is the line through `point = (x0, y0)` with equation
`y - y0 = slope * (x - x0)`. Use `slope = Inf` for a vertical front `x = x0`.
The two domains point along `±domain_direction`.

`wall_normal` is specified in the local frame of the wall front: local `x` is
the direction normal to the front, local `y` is the direction
along the wall, and local `z` is the global out-of-plane direction. For example,
with `domain_direction = (0, 0, 1)`, `wall_normal = (1, 0, 0)` gives a Neel
wall and `wall_normal = (0, 1, 0)` gives a Bloch wall.

This modifier is additive, so it is intended to be applied to a mostly uniform
state, possibly multiple times, followed by [`normalize_spins!`](@ref). The
painted wall profile is corrected by subtracting `background_direction`, which
defaults to `domain_direction`. When painting multiple walls on the same seed,
set `background_direction` to the original uniform seed direction.
"""
function paint_domain_wall2D!(
    spins::AbstractArray{<:AbstractFloat, 3};
    point = ((size(spins, 2) + 1) / 2, (size(spins, 3) + 1) / 2),
    slope::Real = 0.0,
    width::Real = 5.0,
    domain_direction = (0.0, 0.0, 1.0),
    wall_normal = (1.0, 0.0, 0.0),
    background_direction = nothing,
)
    size(spins, 1) == 3 || throw(ArgumentError("Expected a spin array of size (3, Nx, Ny)."))
    width > 0 || throw(ArgumentError("width must be positive."))
    isnan(slope) && throw(ArgumentError("slope must be finite or Inf."))

    dx, dy, dz = _normalized_direction(domain_direction)
    bx, by, bz = isnothing(background_direction) ? (dx, dy, dz) : _normalized_direction(background_direction)
    front_normal, front_tangent = _domain_wall_frame(slope)
    wall_global = _local_wall_direction(wall_normal, front_normal, front_tangent)
    wx, wy, wz = _transverse_direction(wall_global, (dx, dy, dz))

    x0, y0 = point
    @inbounds for j in axes(spins, 3), i in axes(spins, 2)
        xi = _domain_wall_coordinate(i, j, x0, y0, slope) / width
        longitudinal = tanh(xi)
        transverse = sech(xi)

        spins[1, i, j] += longitudinal * dx + transverse * wx - bx
        spins[2, i, j] += longitudinal * dy + transverse * wy - by
        spins[3, i, j] += longitudinal * dz + transverse * wz - bz
    end

    return spins
end

"""
    paint_skyrmion2D!(
        spins;
        center=((size(spins, 2) + 1) / 2, (size(spins, 3) + 1) / 2),
        radius=10.0,
        width=3.0,
        center_direction=-1.0,
        helicity=0.0,
        vorticity=1.0,
    )

Paint a skyrmion in-place on top of a uniformly magnetized 2D background.

The skyrmion is centered at `center = (x0, y0)` with domain-wall radius
`radius` and wall width `width`. `center_direction` sets the spin direction at
the skyrmion core and must be either `+1` or `-1`, corresponding to `+z` or
`-z`. The in-plane
winding is set by `vorticity * phi + helicity`, where `phi = atan(y - y0, x - x0)`.

For a background aligned with `+z`, `helicity = 0` or `π` gives a Neel-type
skyrmion, while `helicity = ±π/2` gives a Bloch-type skyrmion. The sign of
`vorticity` sets the sense of rotation.

This modifier is additive, so it is intended to be applied to a mostly uniform
state, possibly multiple times, followed by [`normalize_spins!`](@ref). The
uniform background is taken to point opposite to `center_direction`, so the
intended use is to start from a seed aligned with that far-field direction.
"""
function paint_skyrmion2D!(
    spins::AbstractArray{<:AbstractFloat, 3};
    center = ((size(spins, 2) + 1) / 2, (size(spins, 3) + 1) / 2),
    radius::Real = 10.0,
    width::Real = 3.0,
    center_direction::Real = -1.0,
    helicity::Real = 0.0,
    vorticity::Real = 1.0,
)
    size(spins, 1) == 3 || throw(ArgumentError("Expected a spin array of size (3, Nx, Ny)."))
    radius >= 0 || throw(ArgumentError("radius must be nonnegative."))
    width > 0 || throw(ArgumentError("width must be positive."))
    abs(center_direction) == 1 || throw(ArgumentError("center_direction must be +1 or -1."))

    dx = 0.0
    dy = 0.0
    dz = -Float64(center_direction)
    x0, y0 = center

    @inbounds for j in axes(spins, 3), i in axes(spins, 2)
        rx = i - x0
        ry = j - y0
        r = hypot(rx, ry)
        phi = atan(ry, rx)
        psi = vorticity * phi + helicity

        qx = cos(psi)
        qy = sin(psi)
        qz = 0.0

        xi = (r - radius) / width
        longitudinal = tanh(xi)
        transverse = sech(xi)

        spins[1, i, j] += longitudinal * dx + transverse * qx - dx
        spins[2, i, j] += longitudinal * dy + transverse * qy - dy
        spins[3, i, j] += longitudinal * dz + transverse * qz - dz
    end

    return spins
end

"""
    normalize_spins!(spins)

Normalize each spin vector in-place for a `(3, Nx)` or `(3, Nx, Ny)` array.
"""
function normalize_spins!(spins::AbstractArray{<:AbstractFloat, 2})
    size(spins, 1) == 3 || throw(ArgumentError("Expected a spin array of size (3, Nx)."))

    @inbounds for i in axes(spins, 2)
        norm_i = sqrt(spins[1, i]^2 + spins[2, i]^2 + spins[3, i]^2)
        norm_i > 0 || throw(ArgumentError("Encountered a zero-length spin at site $i."))
        spins[1, i] /= norm_i
        spins[2, i] /= norm_i
        spins[3, i] /= norm_i
    end

    return spins
end

function normalize_spins!(spins::AbstractArray{<:AbstractFloat, 3})
    size(spins, 1) == 3 || throw(ArgumentError("Expected a spin array of size (3, Nx, Ny)."))

    @inbounds for j in axes(spins, 3), i in axes(spins, 2)
        norm_ij = sqrt(spins[1, i, j]^2 + spins[2, i, j]^2 + spins[3, i, j]^2)
        norm_ij > 0 || throw(ArgumentError("Encountered a zero-length spin at site ($i, $j)."))
        spins[1, i, j] /= norm_ij
        spins[2, i, j] /= norm_ij
        spins[3, i, j] /= norm_ij
    end

    return spins
end

function _normalized_direction(direction)
    length(direction) == 3 || throw(ArgumentError("direction must have exactly three components."))

    ux = Float64(direction[1])
    uy = Float64(direction[2])
    uz = Float64(direction[3])
    norm_u = sqrt(ux^2 + uy^2 + uz^2)
    norm_u > 0 || throw(ArgumentError("direction must be nonzero."))

    return ux / norm_u, uy / norm_u, uz / norm_u
end

function _transverse_direction(direction, axis)
    vx, vy, vz = _normalized_direction(direction)
    ax, ay, az = axis

    projection = vx * ax + vy * ay + vz * az
    tx = vx - projection * ax
    ty = vy - projection * ay
    tz = vz - projection * az
    norm_t = sqrt(tx^2 + ty^2 + tz^2)
    norm_t > 0 || throw(ArgumentError("wall_normal must not be parallel to domain_direction."))

    return tx / norm_t, ty / norm_t, tz / norm_t
end

function _domain_wall_frame(slope)
    if isinf(slope)
        return (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)
    end

    denom = sqrt(1 + slope^2)
    return (-slope / denom, 1.0 / denom, 0.0), (1.0 / denom, slope / denom, 0.0)
end

function _local_wall_direction(direction, front_normal, front_tangent)
    lx, ly, lz = _normalized_direction(direction)
    nx, ny, nz = front_normal
    tx, ty, tz = front_tangent

    return (
        lx * nx + ly * tx,
        lx * ny + ly * ty,
        lx * nz + ly * tz + lz,
    )
end

function _domain_wall_coordinate(i, j, x0, y0, slope)
    if isinf(slope)
        return i - x0
    end

    return (j - y0 - slope * (i - x0)) / sqrt(1 + slope^2)
end
