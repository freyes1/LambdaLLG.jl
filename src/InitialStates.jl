"""
    uniform_state1D(Nx; direction=(0.0, 0.0, 1.0))

Create a 1D spin state with all spins aligned along `direction`.
"""
function uniform_state1D(Nx::Integer; direction = (0.0, 0.0, 1.0))
    Nx > 0 || throw(ArgumentError("Nx must be positive."))
    unit_direction = _normalize_vector(direction, "direction")

    spins = zeros(Float64, 3, Nx)
    @inbounds for i in 1:Nx
        spins[:, i] .= unit_direction
    end

    return spins
end

"""
    uniform_state2D(Nx, Ny; direction=(0.0, 0.0, 1.0))

Create a 2D spin state with all spins aligned along `direction`.
"""
function uniform_state2D(Nx::Integer, Ny::Integer; direction = (0.0, 0.0, 1.0))
    Nx > 0 || throw(ArgumentError("Nx must be positive."))
    Ny > 0 || throw(ArgumentError("Ny must be positive."))
    unit_direction = _normalize_vector(direction, "direction")

    spins = zeros(Float64, 3, Nx, Ny)
    @inbounds for j in 1:Ny, i in 1:Nx
        spins[:, i, j] .= unit_direction
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
        for c in 1:3
            spins[c, i] /= norm_i
        end
    end

    return spins
end

function normalize_spins!(spins::AbstractArray{<:AbstractFloat, 3})
    size(spins, 1) == 3 || throw(ArgumentError("Expected a spin array of size (3, Nx, Ny)."))

    @inbounds for j in axes(spins, 3), i in axes(spins, 2)
        norm_ij = sqrt(spins[1, i, j]^2 + spins[2, i, j]^2 + spins[3, i, j]^2)
        norm_ij > 0 || throw(ArgumentError("Encountered a zero-length spin at site ($i, $j)."))
        for c in 1:3
            spins[c, i, j] /= norm_ij
        end
    end

    return spins
end

"""
    paint_domain_wall!(
        spins;
        point=((Nx + 1) / 2, (Ny + 1) / 2),
        slope=0.0,
        width=5.0,
        wall_type=:neel,
        easy_axis=(0.0, 0.0, 1.0),
        reference_axis=(1.0, 0.0, 0.0),
        chirality=1.0,
        polarity=1.0,
        amplitude=1.0,
    )

Paint a straight 2D domain wall in-place.

The wall front follows the line `y - y0 = slope * (x - x0)` passing through
`point = (x0, y0)`. Set `slope = Inf` for a vertical wall front `x = x0`.
This modifier is additive, so it is intended to be applied to a mostly uniform
background and followed by [`normalize_spins!`](@ref).
"""
function paint_domain_wall!(
    spins::AbstractArray{<:AbstractFloat, 3};
    point = ((size(spins, 2) + 1) / 2, (size(spins, 3) + 1) / 2),
    slope::Real = 0.0,
    width::Real = 5.0,
    wall_type::Symbol = :neel,
    easy_axis = (0.0, 0.0, 1.0),
    reference_axis = (1.0, 0.0, 0.0),
    chirality::Real = 1.0,
    polarity::Real = 1.0,
    amplitude::Real = 1.0,
)
    size(spins, 1) == 3 || throw(ArgumentError("Expected a spin array of size (3, Nx, Ny)."))
    width > 0 || throw(ArgumentError("width must be positive."))

    easy = _normalize_vector(easy_axis, "easy_axis")
    neel_axis, bloch_axis = _transverse_axes(easy, reference_axis)
    transverse_axis = _wall_axis(wall_type, neel_axis, bloch_axis)

    x0, y0 = point
    @inbounds for j in axes(spins, 3), i in axes(spins, 2)
        xi = _line_coordinate(i, j, x0, y0, slope) / width
        profile = polarity * tanh(xi) .* easy .+ chirality * sech(xi) .* transverse_axis
        spins[:, i, j] .+= amplitude .* (profile .- easy)
    end

    return spins
end

function _normalize_vector(vec, name::AbstractString)
    length(vec) == 3 || throw(ArgumentError("$name must have exactly three components."))

    out = Float64[vec[1], vec[2], vec[3]]
    norm_out = norm(out)
    norm_out > 0 || throw(ArgumentError("$name must be nonzero."))
    out ./= norm_out
    return out
end

function _transverse_axes(easy_axis, reference_axis)
    ref = _normalize_vector(reference_axis, "reference_axis")
    neel_axis = ref .- dot(ref, easy_axis) .* easy_axis
    neel_norm = norm(neel_axis)
    neel_norm > sqrt(eps(Float64)) || throw(ArgumentError("reference_axis must not be parallel to easy_axis."))
    neel_axis ./= neel_norm

    bloch_axis = cross(easy_axis, neel_axis)
    bloch_axis ./= norm(bloch_axis)

    return neel_axis, bloch_axis
end

function _wall_axis(wall_type::Symbol, neel_axis, bloch_axis)
    if wall_type == :neel
        return neel_axis
    elseif wall_type == :bloch
        return bloch_axis
    end

    throw(ArgumentError("Unsupported wall type `$wall_type`. Use `:neel` or `:bloch`."))
end

function _line_coordinate(i, j, x0, y0, slope)
    if isinf(slope)
        return i - x0
    end

    return (j - y0 - slope * (i - x0)) / sqrt(1 + slope^2)
end
