using Test
using LinearAlgebra
using LambdaLLG

@testset "1D bulk DMI" begin
    spins = [
        1.0  0.3  -0.4  0.8;
        0.2 -0.6   0.5 -0.1;
       -0.3  0.7   0.1  0.4
    ]

    D = (1.2, -0.7, 0.3)
    p = LLGParams1D(
        size(spins, 2),
        0.0,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        0.0;
        D = D,
    )

    fill!(p.fields.Beff, 0.0)
    LambdaLLG.add_bulk_dmi1D!(spins, p)

    expected = zeros(3, size(spins, 2))
    Dvec = collect(D)

    for i in 1:size(spins, 2)
        left = i > 1 ? spins[:, i - 1] : zeros(3)
        right = i < size(spins, 2) ? spins[:, i + 1] : zeros(3)
        expected[:, i] .= cross(Dvec, right - left)
    end

    @test isapprox(p.fields.Beff, expected)
end

@testset "2D bulk DMI" begin
    spins = zeros(3, 3, 4)
    spins[:, 1, 1] .= (0.1, -0.2, 0.7)
    spins[:, 2, 1] .= (-0.4, 0.3, 0.5)
    spins[:, 3, 1] .= (0.6, 0.1, -0.2)
    spins[:, 1, 2] .= (-0.3, 0.8, -0.1)
    spins[:, 2, 2] .= (0.5, -0.7, 0.2)
    spins[:, 3, 2] .= (-0.6, 0.4, 0.3)
    spins[:, 1, 3] .= (0.2, 0.5, -0.4)
    spins[:, 2, 3] .= (-0.1, -0.3, 0.9)
    spins[:, 3, 3] .= (0.7, -0.5, -0.2)
    spins[:, 1, 4] .= (-0.8, 0.2, 0.1)
    spins[:, 2, 4] .= (0.4, 0.6, -0.3)
    spins[:, 3, 4] .= (-0.2, -0.4, 0.8)

    D = ((1.1, 0.0, 0.0), (0.0, 1.1, 0.0))
    p = LLGParams2D(
        size(spins, 2),
        size(spins, 3),
        0.0,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        0.0;
        D = D,
    )

    fill!(p.fields.Beff, 0.0)
    LambdaLLG.add_dmi2D!(spins, p)

    expected = zeros(size(spins))
    Dx = collect(D[1])
    Dy = collect(D[2])

    for j in 1:size(spins, 3), i in 1:size(spins, 2)
        left = i > 1 ? spins[:, i - 1, j] : zeros(3)
        right = i < size(spins, 2) ? spins[:, i + 1, j] : zeros(3)
        down = j > 1 ? spins[:, i, j - 1] : zeros(3)
        up = j < size(spins, 3) ? spins[:, i, j + 1] : zeros(3)
        expected[:, i, j] .= cross(Dx, right - left) + cross(Dy, up - down)
    end

    @test isapprox(p.fields.Beff, expected)
end

@testset "2D interface DMI" begin
    spins = zeros(3, 4, 3)
    spins[:, 1, 1] .= (0.2, 0.1, -0.5)
    spins[:, 2, 1] .= (-0.7, 0.4, 0.3)
    spins[:, 3, 1] .= (0.6, -0.2, 0.1)
    spins[:, 4, 1] .= (-0.1, 0.8, -0.4)
    spins[:, 1, 2] .= (0.5, -0.6, 0.2)
    spins[:, 2, 2] .= (-0.3, 0.7, -0.8)
    spins[:, 3, 2] .= (0.4, 0.2, 0.9)
    spins[:, 4, 2] .= (-0.5, -0.1, 0.6)
    spins[:, 1, 3] .= (0.7, -0.4, 0.3)
    spins[:, 2, 3] .= (-0.2, 0.5, -0.7)
    spins[:, 3, 3] .= (0.1, -0.8, 0.4)
    spins[:, 4, 3] .= (0.6, 0.3, -0.2)

    D0 = 0.9
    D = ((0.0, D0, 0.0), (-D0, 0.0, 0.0))
    p = LLGParams2D(
        size(spins, 2),
        size(spins, 3),
        0.0,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        0.0;
        D = D,
    )

    fill!(p.fields.Beff, 0.0)
    LambdaLLG.add_dmi2D!(spins, p)

    expected = zeros(size(spins))
    Dx = collect(D[1])
    Dy = collect(D[2])

    for j in 1:size(spins, 3), i in 1:size(spins, 2)
        left = i > 1 ? spins[:, i - 1, j] : zeros(3)
        right = i < size(spins, 2) ? spins[:, i + 1, j] : zeros(3)
        down = j > 1 ? spins[:, i, j - 1] : zeros(3)
        up = j < size(spins, 3) ? spins[:, i, j + 1] : zeros(3)
        expected[:, i, j] .= cross(Dx, right - left) + cross(Dy, up - down)
    end

    @test isapprox(p.fields.Beff, expected)
end

@testset "Uniform magnetized seeds" begin
    spins1D = uniform_seed1D(5; direction = (1.0, 2.0, -2.0))
    expected1D = [1.0, 2.0, -2.0] ./ norm([1.0, 2.0, -2.0])
    for i in 1:5
        @test isapprox(spins1D[:, i], expected1D)
    end

    spins2D = uniform_seed2D(4, 3; direction = (0.0, 3.0, 4.0))
    expected2D = [0.0, 3.0, 4.0] ./ norm([0.0, 3.0, 4.0])
    for j in 1:3, i in 1:4
        @test isapprox(spins2D[:, i, j], expected2D)
    end

    @test size(spins1D) == (3, 5)
    @test size(spins2D) == (3, 4, 3)

    @test_throws ArgumentError uniform_seed1D(0)
    @test_throws ArgumentError uniform_seed2D(4, 0)
    @test_throws ArgumentError uniform_seed1D(4; direction = (0.0, 0.0, 0.0))
end

@testset "1D domain-wall painter" begin
    neel = uniform_seed1D(13; direction = (0.0, 0.0, 1.0))
    paint_domain_wall1D!(
        neel;
        center = 7.0,
        width = 1.5,
        domain_direction = (0.0, 0.0, 1.0),
        wall_normal = (1.0, 0.0, 0.0),
    )
    normalize_spins!(neel)

    @test neel[3, 3] < -0.95
    @test neel[3, 11] > 0.95
    @test neel[1, 7] > 0.99
    @test abs(neel[2, 7]) < 1e-12

    bloch = uniform_seed1D(13; direction = (0.0, 0.0, 1.0))
    paint_domain_wall1D!(
        bloch;
        center = 7.0,
        width = 1.5,
        domain_direction = (0.0, 0.0, 1.0),
        wall_normal = (0.0, 1.0, 0.0),
    )
    normalize_spins!(bloch)

    @test bloch[3, 3] < -0.95
    @test bloch[3, 11] > 0.95
    @test abs(bloch[1, 7]) < 1e-12
    @test bloch[2, 7] > 0.99

    tilted = uniform_seed1D(11; direction = (1.0, 1.0, 0.0))
    paint_domain_wall1D!(
        tilted;
        center = 6.0,
        width = 1.2,
        domain_direction = (1.0, 1.0, 0.0),
        wall_normal = (1.0, 0.0, 1.0),
    )
    normalize_spins!(tilted)

    domain = [1.0, 1.0, 0.0] ./ norm([1.0, 1.0, 0.0])
    @test isapprox(dot(tilted[:, 2], domain), -1.0; atol = 5e-2)
    @test isapprox(dot(tilted[:, 10], domain), 1.0; atol = 5e-2)

    layered = uniform_seed1D(21; direction = (0.0, 0.0, 1.0))
    paint_domain_wall1D!(layered; center = 7.0, width = 1.5, wall_normal = (1.0, 0.0, 0.0))
    paint_domain_wall1D!(
        layered;
        center = 15.0,
        width = 1.5,
        domain_direction = (0.0, 0.0, -1.0),
        wall_normal = (1.0, 0.0, 0.0),
        background_direction = (0.0, 0.0, 1.0),
    )
    layered_norm_error = maximum(abs(norm(layered[:, i]) - 1.0) for i in axes(layered, 2))
    @test layered_norm_error > 1e-3
    normalize_spins!(layered)
    @test all(isapprox(norm(layered[:, i]), 1.0; atol = 1e-12) for i in axes(layered, 2))
    @test layered[3, 3] < -0.95
    @test layered[3, 11] > 0.95
    @test layered[3, 19] < -0.95

    bad = uniform_seed1D(5)
    @test_throws ArgumentError paint_domain_wall1D!(bad; width = 0.0)
    @test_throws ArgumentError paint_domain_wall1D!(bad; wall_normal = (0.0, 0.0, 2.0))
end

@testset "STT activation uses parameter flag" begin
    spins1D = zeros(3, 9)
    for i in 1:9
        theta = 0.45 + 0.03 * i
        phi = 0.35 * i
        spins1D[:, i] .= (sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta))
    end

    p1_off = LLGParams1D(
        size(spins1D, 2),
        0.0,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        0.0;
        u_stt = 0.2,
        beta_stt = 0.3,
        stt_active = false,
    )
    p1_on = LLGParams1D(
        size(spins1D, 2),
        0.0,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        0.0;
        u_stt = 0.2,
        beta_stt = 0.3,
        stt_active = true,
    )

    rhs1_off = copy(LambdaLLG.rhs1D!(copy(spins1D), p1_off, 0.0))
    rhs1_on = copy(LambdaLLG.rhs1D!(copy(spins1D), p1_on, 0.0))

    @test p1_off.stt_active == false
    @test p1_on.stt_active == true
    @test norm(rhs1_on - rhs1_off) > 1e-10

    spins2D = zeros(3, 4, 5)
    for j in 1:5, i in 1:4
        theta = 0.35 + 0.04 * i
        phi = 0.25 * i + 0.15 * j
        spins2D[:, i, j] .= (sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta))
    end

    p2_off = LLGParams2D(
        size(spins2D, 2),
        size(spins2D, 3),
        0.0,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        0.0;
        u_stt = (0.15, -0.05),
        beta_stt = 0.25,
        stt_active = false,
    )
    p2_on = LLGParams2D(
        size(spins2D, 2),
        size(spins2D, 3),
        0.0,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        0.0;
        u_stt = (0.15, -0.05),
        beta_stt = 0.25,
        stt_active = true,
    )

    rhs2_off = copy(LambdaLLG.rhs2D!(copy(spins2D), p2_off, 0.0))
    rhs2_on = copy(LambdaLLG.rhs2D!(copy(spins2D), p2_on, 0.0))

    @test p2_off.stt_active == false
    @test p2_on.stt_active == true
    @test norm(rhs2_on - rhs2_off) > 1e-10

    wall = uniform_seed1D(21; direction = (0.0, 0.0, 1.0))
    paint_domain_wall1D!(wall; center = 11.0, width = 2.0, wall_normal = (0.0, 1.0, 0.0))
    normalize_spins!(wall)

    psolve_off = LLGParams1D(
        21,
        0.0,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        0.0;
        u_stt = 0.2,
        beta_stt = 0.3,
        stt_active = false,
    )
    psolve_on = LLGParams1D(
        21,
        0.0,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        0.0;
        u_stt = 0.2,
        beta_stt = 0.3,
        stt_active = true,
    )

    sol_off = evolve1D(copy(wall), (0.0, 0.5), psolve_off)
    sol_on = evolve1D(copy(wall), (0.0, 0.5), psolve_on)

    @test psolve_off.stt_active == false
    @test psolve_on.stt_active == true
    @test !isapprox(sol_off.u[end], sol_on.u[end]; atol = 1e-8, rtol = 1e-8)
end

@testset "2D domain-wall painter" begin
    horizontal_neel = uniform_seed2D(13, 11; direction = (0.0, 0.0, 1.0))
    paint_domain_wall2D!(
        horizontal_neel;
        point = (7.0, 6.0),
        slope = 0.0,
        width = 1.5,
        domain_direction = (0.0, 0.0, 1.0),
        wall_normal = (1.0, 0.0, 0.0),
    )
    normalize_spins!(horizontal_neel)

    @test horizontal_neel[3, 7, 2] < -0.95
    @test horizontal_neel[3, 7, 10] > 0.95
    @test abs(horizontal_neel[1, 7, 6]) < 1e-12
    @test horizontal_neel[2, 7, 6] > 0.99

    horizontal_bloch = uniform_seed2D(13, 11; direction = (0.0, 0.0, 1.0))
    paint_domain_wall2D!(
        horizontal_bloch;
        point = (7.0, 6.0),
        slope = 0.0,
        width = 1.5,
        domain_direction = (0.0, 0.0, 1.0),
        wall_normal = (0.0, 1.0, 0.0),
    )
    normalize_spins!(horizontal_bloch)

    @test horizontal_bloch[1, 7, 6] > 0.99
    @test abs(horizontal_bloch[2, 7, 6]) < 1e-12

    vertical_neel = uniform_seed2D(13, 11; direction = (0.0, 0.0, 1.0))
    paint_domain_wall2D!(
        vertical_neel;
        point = (7.0, 5.0),
        slope = Inf,
        width = 1.5,
        domain_direction = (0.0, 0.0, 1.0),
        wall_normal = (1.0, 0.0, 0.0),
    )
    normalize_spins!(vertical_neel)

    @test vertical_neel[3, 3, 5] < -0.95
    @test vertical_neel[3, 11, 5] > 0.95
    @test vertical_neel[1, 7, 5] > 0.99
    @test abs(vertical_neel[2, 7, 5]) < 1e-12

    sloped = uniform_seed2D(13, 13; direction = (0.0, 0.0, 1.0))
    paint_domain_wall2D!(
        sloped;
        point = (6.0, 6.0),
        slope = 1.0,
        width = 1.5,
        domain_direction = (0.0, 0.0, 1.0),
        wall_normal = (1.0, 0.0, 0.0),
    )
    normalize_spins!(sloped)

    expected_center = [-1.0, 1.0, 0.0] ./ sqrt(2.0)
    @test isapprox(sloped[:, 6, 6], expected_center; atol = 1e-12)
    @test sloped[3, 4, 8] > 0.95
    @test sloped[3, 8, 4] < -0.95

    layered = uniform_seed2D(21, 21; direction = (0.0, 0.0, 1.0))
    paint_domain_wall2D!(layered; point = (11.0, 6.0), slope = 0.0, width = 1.5, wall_normal = (1.0, 0.0, 0.0))
    paint_domain_wall2D!(
        layered;
        point = (11.0, 16.0),
        slope = 0.0,
        width = 1.5,
        domain_direction = (0.0, 0.0, -1.0),
        wall_normal = (1.0, 0.0, 0.0),
        background_direction = (0.0, 0.0, 1.0),
    )
    layered_norm_error = maximum(abs(norm(layered[:, i, j]) - 1.0) for j in axes(layered, 3), i in axes(layered, 2))
    @test layered_norm_error > 1e-3
    normalize_spins!(layered)
    @test all(isapprox(norm(layered[:, i, j]), 1.0; atol = 1e-12) for j in axes(layered, 3), i in axes(layered, 2))
    @test layered[3, 11, 2] < -0.95
    @test layered[3, 11, 11] > 0.95
    @test layered[3, 11, 20] < -0.95

    bad = uniform_seed2D(5, 5)
    @test_throws ArgumentError paint_domain_wall2D!(bad; width = 0.0)
    @test_throws ArgumentError paint_domain_wall2D!(bad; wall_normal = (0.0, 0.0, 2.0))
end

@testset "2D skyrmion painter" begin
    neel = uniform_seed2D(41, 41; direction = (0.0, 0.0, 1.0))
    paint_skyrmion2D!(
        neel;
        center = (21.0, 21.0),
        radius = 8.0,
        width = 2.0,
        center_direction = -1.0,
        helicity = 0.0,
        vorticity = 1.0,
    )
    normalize_spins!(neel)

    @test neel[3, 21, 21] < -0.99
    @test neel[3, 21, 1] > 0.99
    @test neel[1, 29, 21] > 0.95
    @test abs(neel[2, 29, 21]) < 1e-12
    @test neel[2, 21, 29] > 0.95
    @test abs(neel[1, 21, 29]) < 1e-12

    bloch = uniform_seed2D(41, 41; direction = (0.0, 0.0, 1.0))
    paint_skyrmion2D!(
        bloch;
        center = (21.0, 21.0),
        radius = 8.0,
        width = 2.0,
        center_direction = -1.0,
        helicity = pi / 2,
        vorticity = 1.0,
    )
    normalize_spins!(bloch)

    @test abs(bloch[1, 29, 21]) < 1e-12
    @test bloch[2, 29, 21] > 0.95
    @test bloch[1, 21, 29] < -0.95
    @test abs(bloch[2, 21, 29]) < 1e-12

    reverse_bg = uniform_seed2D(41, 41; direction = (0.0, 0.0, -1.0))
    paint_skyrmion2D!(
        reverse_bg;
        center = (21.0, 21.0),
        radius = 8.0,
        width = 2.0,
        center_direction = 1.0,
        helicity = 0.0,
        vorticity = 1.0,
    )
    normalize_spins!(reverse_bg)

    @test reverse_bg[3, 21, 21] > 0.99
    @test reverse_bg[3, 21, 1] < -0.99

    layered = uniform_seed2D(61, 61; direction = (0.0, 0.0, 1.0))
    paint_skyrmion2D!(layered; center = (21.0, 31.0), radius = 7.0, width = 2.0, center_direction = -1.0, helicity = 0.0, vorticity = 1.0)
    paint_skyrmion2D!(layered; center = (41.0, 31.0), radius = 7.0, width = 2.0, center_direction = -1.0, helicity = pi / 2, vorticity = -1.0)
    layered_norm_error = maximum(abs(norm(layered[:, i, j]) - 1.0) for j in axes(layered, 3), i in axes(layered, 2))
    @test layered_norm_error > 1e-3
    normalize_spins!(layered)
    @test all(isapprox(norm(layered[:, i, j]), 1.0; atol = 1e-12) for j in axes(layered, 3), i in axes(layered, 2))

    bad = uniform_seed2D(5, 5)
    @test_throws ArgumentError paint_skyrmion2D!(bad; radius = -1.0)
    @test_throws ArgumentError paint_skyrmion2D!(bad; width = 0.0)
    @test_throws ArgumentError paint_skyrmion2D!(bad; center_direction = 0.0)
end

@testset "2D spin-transfer torque" begin
    spins = zeros(3, 4, 3)
    spins[:, 1, 1] .= (0.8, 0.1, 0.6)
    spins[:, 2, 1] .= (0.2, 0.9, 0.4)
    spins[:, 3, 1] .= (-0.4, 0.3, 0.85)
    spins[:, 4, 1] .= (-0.6, -0.2, 0.77)
    spins[:, 1, 2] .= (0.7, -0.5, 0.5)
    spins[:, 2, 2] .= (0.1, 0.8, 0.6)
    spins[:, 3, 2] .= (-0.3, 0.4, 0.85)
    spins[:, 4, 2] .= (-0.5, -0.1, 0.8)
    spins[:, 1, 3] .= (0.6, -0.4, 0.65)
    spins[:, 2, 3] .= (0.0, 0.7, 0.7)
    spins[:, 3, 3] .= (-0.2, 0.5, 0.84)
    spins[:, 4, 3] .= (-0.4, 0.2, 0.89)
    normalize_spins!(spins)

    p = LLGParams2D(
        size(spins, 2),
        size(spins, 3),
        0.0,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        0.0;
        u_stt = (0.35, -0.2),
        beta_stt = 0.4,
    )

    fill!(p.fields.dS_1, 0.0)
    LambdaLLG.add_stt2D!(spins, p)

    expected = zeros(size(spins))
    ux, uy = p.u_stt

    for j in axes(spins, 3), i in axes(spins, 2)
        grad_x = if size(spins, 2) == 1
            zeros(3)
        elseif i == 1
            spins[:, 2, j] - spins[:, 1, j]
        elseif i == size(spins, 2)
            spins[:, end, j] - spins[:, end - 1, j]
        else
            0.5 .* (spins[:, i + 1, j] - spins[:, i - 1, j])
        end

        grad_y = if size(spins, 3) == 1
            zeros(3)
        elseif j == 1
            spins[:, i, 2] - spins[:, i, 1]
        elseif j == size(spins, 3)
            spins[:, i, end] - spins[:, i, end - 1]
        else
            0.5 .* (spins[:, i, j + 1] - spins[:, i, j - 1])
        end

        directional_grad = ux .* grad_x .+ uy .* grad_y
        expected[:, i, j] .= -directional_grad .+ p.beta_stt .* cross(spins[:, i, j], directional_grad)
    end

    @test isapprox(p.fields.dS_1, expected)
end

@testset "2D STT parameter flag" begin
    spins = zeros(3, 5, 5)
    x0 = 3.0
    y0 = 3.0
    for j in 1:5, i in 1:5
        dx = i - x0
        dy = j - y0
        r = hypot(dx, dy)
        mz = -tanh((1.8 - r) / 0.9)
        mxy = sqrt(max(0.0, 1.0 - mz^2))
        if r < 1e-12
            spins[:, i, j] .= (0.0, 0.0, -1.0)
        else
            spins[:, i, j] .= (mxy * dx / r, mxy * dy / r, mz)
        end
    end
    normalize_spins!(spins)

    p_off = LLGParams2D(
        5,
        5,
        0.0,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        0.0;
        u_stt = (0.2, 0.1),
        beta_stt = 0.15,
        stt_active = false,
    )
    p_on = LLGParams2D(
        5,
        5,
        0.0,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        0.0;
        u_stt = (0.2, 0.1),
        beta_stt = 0.15,
        stt_active = true,
    )

    sol_off = evolve2D(copy(spins), (0.0, 0.2), p_off; reltol = 1e-7, abstol = 1e-7)
    @test p_off.stt_active == false

    sol_on = evolve2D(copy(spins), (0.0, 0.2), p_on; reltol = 1e-7, abstol = 1e-7)
    @test p_on.stt_active == true

    final_off = reshape(sol_off.u[end], size(spins))
    final_on = reshape(sol_on.u[end], size(spins))

    @test isapprox(final_off, spins; atol = 1e-10, rtol = 1e-10)
    @test maximum(abs.(final_on .- final_off)) > 1e-5
end
