"""
    uniform_state1D(Nx; direction=(0.0, 0.0, 1.0), magnitude=1.0)

Create a 1D spin state with every site aligned along `direction`.

The input `direction` is normalized internally, and each spin is returned with
length `magnitude`.
"""
function uniform_state1D(
    Nx::Integer;
    direction = (0.0, 0.0, 1.0),
    magnitude::Real = 1.0,
)
    Nx > 0 || throw(ArgumentError("Nx must be positive."))
    state_direction = magnitude .* _normalize_vector(direction, "direction")

    spins = zeros(Float64, 3, Nx)
    @inbounds for i in 1:Nx
        spins[:, i] .= state_direction
    end

    return spins
end

"""
    uniform_state2D(Nx, Ny; direction=(0.0, 0.0, 1.0), magnitude=1.0)

Create a 2D spin state with every site aligned along `direction`.

The input `direction` is normalized internally, and each spin is returned with
length `magnitude`.
"""
function uniform_state2D(
    Nx::Integer,
    Ny::Integer;
    direction = (0.0, 0.0, 1.0),
    magnitude::Real = 1.0,
)
    Nx > 0 || throw(ArgumentError("Nx must be positive."))
    Ny > 0 || throw(ArgumentError("Ny must be positive."))
    state_direction = magnitude .* _normalize_vector(direction, "direction")

    spins = zeros(Float64, 3, Nx, Ny)
    @inbounds for j in 1:Ny
        for i in 1:Nx
            spins[:, i, j] .= state_direction
        end
    end

    return spins
end

"""
    normalize_spins!(spins)

Normalize each spin vector in-place for either a 1D `(3, Nx)` state or a 2D
`(3, Nx, Ny)` state.
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

    @inbounds for j in axes(spins, 3)
        for i in axes(spins, 2)
            norm_ij = sqrt(spins[1, i, j]^2 + spins[2, i, j]^2 + spins[3, i, j]^2)
            norm_ij > 0 || throw(ArgumentError("Encountered a zero-length spin at site ($i, $j)."))
            for c in 1:3
                spins[c, i, j] /= norm_ij
            end
        end
    end

    return spins
end

"""
    paint_domain_wall!(spins; kwargs...)

Paint a 1D or 2D domain-wall ansatz additively onto an existing spin state.

This modifier is intended to be layered on top of a mostly uniform background
and followed by [`normalize_spins!`](@ref). Supported wall types are `:neel`
and `:bloch`.
"""
function paint_domain_wall!(
    spins::AbstractArray{<:AbstractFloat, 2};
    center::Real = (size(spins, 2) + 1) / 2,
    width::Real = 5.0,
    wall_type::Symbol = :neel,
    easy_axis = nothing,
    reference_axis = nothing,
    chirality::Real = 1.0,
    polarity::Real = 1.0,
    background = nothing,
    amplitude::Real = 1.0,
)
    size(spins, 1) == 3 || throw(ArgumentError("Expected a spin array of size (3, Nx)."))
    width > 0 || throw(ArgumentError("width must be positive."))

    background_axis, wall_easy_axis = _resolve_texture_axes(spins, easy_axis, background)
    neel_axis, bloch_axis = _transverse_axes(wall_easy_axis, reference_axis)
    transverse_axis = _wall_axis(wall_type, neel_axis, bloch_axis)

    @inbounds for i in axes(spins, 2)
        xi = (i - center) / width
        profile = polarity * tanh(xi) .* wall_easy_axis .+ chirality * sech(xi) .* transverse_axis
        spins[:, i] .+= amplitude .* (profile .- background_axis)
    end

    return spins
end

function paint_domain_wall!(
    spins::AbstractArray{<:AbstractFloat, 3};
    center = ((size(spins, 2) + 1) / 2, (size(spins, 3) + 1) / 2),
    width::Real = 5.0,
    normal = (1.0, 0.0),
    wall_type::Symbol = :neel,
    easy_axis = nothing,
    reference_axis = nothing,
    chirality::Real = 1.0,
    polarity::Real = 1.0,
    background = nothing,
    amplitude::Real = 1.0,
)
    size(spins, 1) == 3 || throw(ArgumentError("Expected a spin array of size (3, Nx, Ny)."))
    width > 0 || throw(ArgumentError("width must be positive."))

    center_x, center_y = center
    wall_normal = _normalize_spatial_vector(normal, "normal")

    background_axis, wall_easy_axis = _resolve_texture_axes(spins, easy_axis, background)
    neel_axis, bloch_axis = _transverse_axes(wall_easy_axis, reference_axis)
    transverse_axis = _wall_axis(wall_type, neel_axis, bloch_axis)

    @inbounds for j in axes(spins, 3)
        for i in axes(spins, 2)
            xi = ((i - center_x) * wall_normal[1] + (j - center_y) * wall_normal[2]) / width
            profile = polarity * tanh(xi) .* wall_easy_axis .+ chirality * sech(xi) .* transverse_axis
            spins[:, i, j] .+= amplitude .* (profile .- background_axis)
        end
    end

    return spins
end

"""
    paint_skyrmion!(spins; kwargs...)

Paint a 2D skyrmion ansatz additively onto an existing spin state.

The far-field direction follows `background` when provided, otherwise a boundary
average of the current state is used. Supported skyrmion types are `:neel` and
`:bloch`.
"""
function paint_skyrmion!(
    spins::AbstractArray{<:AbstractFloat, 3};
    center = ((size(spins, 2) + 1) / 2, (size(spins, 3) + 1) / 2),
    radius::Real = min(size(spins, 2), size(spins, 3)) / 6,
    width::Real = max(radius / 4, 1.0),
    skyrmion_type::Symbol = :neel,
    easy_axis = nothing,
    reference_axis = nothing,
    chirality::Real = 1.0,
    background = nothing,
    amplitude::Real = 1.0,
)
    size(spins, 1) == 3 || throw(ArgumentError("Expected a spin array of size (3, Nx, Ny)."))
    radius > 0 || throw(ArgumentError("radius must be positive."))
    width > 0 || throw(ArgumentError("width must be positive."))

    center_x, center_y = center
    background_axis, skyrmion_axis = _resolve_texture_axes(spins, easy_axis, background)
    neel_axis, bloch_axis = _transverse_axes(skyrmion_axis, reference_axis)
    use_neel = _is_neel(skyrmion_type)

    @inbounds for j in axes(spins, 3)
        for i in axes(spins, 2)
            dx = i - center_x
            dy = j - center_y
            r = hypot(dx, dy)
            phi = r == 0 ? 0.0 : atan(dy, dx)

            radial_axis = cos(phi) .* neel_axis .+ sin(phi) .* bloch_axis
            tangential_axis = -sin(phi) .* neel_axis .+ cos(phi) .* bloch_axis
            transverse_axis = use_neel ? radial_axis : tangential_axis

            theta = r == 0 ? pi : 2 * atan(exp((radius - r) / width))
            profile = cos(theta) .* skyrmion_axis .+ chirality * sin(theta) .* transverse_axis
            spins[:, i, j] .+= amplitude .* (profile .- background_axis)
        end
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

function _normalize_spatial_vector(vec, name::AbstractString)
    length(vec) == 2 || throw(ArgumentError("$name must have exactly two components."))

    out = Float64[vec[1], vec[2]]
    norm_out = norm(out)
    norm_out > 0 || throw(ArgumentError("$name must be nonzero."))

    out ./= norm_out
    return out
end

function _resolve_texture_axes(spins, easy_axis, background)
    background_axis = if background === nothing
        if easy_axis === nothing
            _boundary_background(spins)
        else
            _normalize_vector(easy_axis, "easy_axis")
        end
    else
        _normalize_vector(background, "background")
    end

    texture_axis = easy_axis === nothing ? background_axis : _normalize_vector(easy_axis, "easy_axis")
    return background_axis, texture_axis
end

function _boundary_background(spins::AbstractArray{<:AbstractFloat, 2})
    guess = Float64[
        spins[1, begin] + spins[1, end],
        spins[2, begin] + spins[2, end],
        spins[3, begin] + spins[3, end],
    ]
    return _fallback_background(guess, vec(sum(spins, dims = 2)))
end

function _boundary_background(spins::AbstractArray{<:AbstractFloat, 3})
    guess = Float64[
        spins[1, begin, begin] + spins[1, end, begin] + spins[1, begin, end] + spins[1, end, end],
        spins[2, begin, begin] + spins[2, end, begin] + spins[2, begin, end] + spins[2, end, end],
        spins[3, begin, begin] + spins[3, end, begin] + spins[3, begin, end] + spins[3, end, end],
    ]
    return _fallback_background(guess, vec(sum(spins, dims = (2, 3))))
end

function _fallback_background(guess, fallback)
    guess_norm = norm(guess)
    if guess_norm > sqrt(eps(Float64))
        return guess ./ guess_norm
    end

    fallback_norm = norm(fallback)
    fallback_norm > 0 || throw(ArgumentError("Could not infer a background direction. Pass `background=` explicitly."))

    return fallback ./ fallback_norm
end

function _transverse_axes(easy_axis, reference_axis)
    ref_axis = _reference_axis(easy_axis, reference_axis)
    neel_axis = ref_axis .- dot(ref_axis, easy_axis) .* easy_axis
    neel_norm = norm(neel_axis)
    neel_norm > sqrt(eps(Float64)) || throw(ArgumentError("reference_axis must not be parallel to easy_axis."))
    neel_axis ./= neel_norm

    bloch_axis = cross(easy_axis, neel_axis)
    bloch_axis ./= norm(bloch_axis)

    return neel_axis, bloch_axis
end

function _reference_axis(easy_axis, reference_axis)
    if reference_axis !== nothing
        return _normalize_vector(reference_axis, "reference_axis")
    end

    candidates = (
        Float64[1.0, 0.0, 0.0],
        Float64[0.0, 1.0, 0.0],
        Float64[0.0, 0.0, 1.0],
    )
    dots = (abs(dot(candidates[1], easy_axis)), abs(dot(candidates[2], easy_axis)), abs(dot(candidates[3], easy_axis)))
    return copy(candidates[argmin(dots)])
end

function _wall_axis(wall_type::Symbol, neel_axis, bloch_axis)
    if _is_neel(wall_type)
        return neel_axis
    elseif _is_bloch(wall_type)
        return bloch_axis
    end

    throw(ArgumentError("Unsupported wall type `$wall_type`. Use `:neel` or `:bloch`."))
end

_is_neel(texture_type::Symbol) = texture_type == :neel
_is_bloch(texture_type::Symbol) = texture_type == :bloch
