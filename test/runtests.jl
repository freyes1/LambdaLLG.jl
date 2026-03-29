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

@testset "Uniform states and normalization" begin
    spins1D = uniform_state1D(5; direction = (1.0, 2.0, -2.0))
    expected1D = [1.0, 2.0, -2.0] ./ norm([1.0, 2.0, -2.0])
    for i in 1:5
        @test isapprox(spins1D[:, i], expected1D)
    end

    spins2D = uniform_state2D(4, 3; direction = (0.0, 3.0, 4.0))
    expected2D = [0.0, 3.0, 4.0] ./ norm([0.0, 3.0, 4.0])
    for j in 1:3, i in 1:4
        @test isapprox(spins2D[:, i, j], expected2D)
    end

    raw = 3.0 .* copy(spins2D)
    normalize_spins!(raw)
    for j in 1:3, i in 1:4
        @test isapprox(norm(raw[:, i, j]), 1.0)
        @test isapprox(raw[:, i, j], expected2D)
    end
end

@testset "2D straight domain wall state" begin
    horizontal = uniform_state2D(13, 11; direction = (0.0, 0.0, 1.0))
    paint_domain_wall!(
        horizontal;
        point = (7.0, 6.0),
        slope = 0.0,
        width = 1.5,
        wall_type = :bloch,
        reference_axis = (1.0, 0.0, 0.0),
    )
    normalize_spins!(horizontal)

    @test horizontal[3, 7, 2] < -0.95
    @test horizontal[3, 7, 10] > 0.95
    @test abs(horizontal[1, 7, 6]) < 1e-10
    @test horizontal[2, 7, 6] > 0.99

    sloped = uniform_state2D(13, 13; direction = (0.0, 0.0, 1.0))
    paint_domain_wall!(
        sloped;
        point = (6.0, 6.0),
        slope = 1.0,
        width = 1.5,
        wall_type = :neel,
        reference_axis = (1.0, 0.0, 0.0),
    )
    normalize_spins!(sloped)

    @test sloped[1, 6, 6] > 0.99
    @test sloped[3, 4, 8] > 0.95
    @test sloped[3, 8, 4] < -0.95

    vertical = uniform_state2D(13, 11; direction = (0.0, 0.0, 1.0))
    paint_domain_wall!(
        vertical;
        point = (7.0, 5.0),
        slope = Inf,
        width = 1.5,
        wall_type = :neel,
        reference_axis = (1.0, 0.0, 0.0),
    )
    normalize_spins!(vertical)

    @test vertical[3, 3, 5] < -0.95
    @test vertical[3, 11, 5] > 0.95
    @test vertical[1, 7, 5] > 0.99

    layered = uniform_state2D(13, 11; direction = (0.0, 0.0, 1.0))
    paint_domain_wall!(
        layered;
        point = (7.0, 4.0),
        slope = 0.0,
        width = 1.5,
        wall_type = :bloch,
        reference_axis = (1.0, 0.0, 0.0),
    )
    paint_domain_wall!(
        layered;
        point = (7.0, 8.0),
        slope = 0.0,
        width = 1.5,
        wall_type = :bloch,
        reference_axis = (1.0, 0.0, 0.0),
        chirality = -1.0,
    )

    layered_norm_error = maximum(abs(norm(layered[:, i, j]) - 1.0) for j in axes(layered, 3), i in axes(layered, 2))
    @test layered_norm_error > 1e-3
    normalize_spins!(layered)
    @test all(isapprox(norm(layered[:, i, j]), 1.0; atol = 1e-12) for j in axes(layered, 3), i in axes(layered, 2))
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

@testset "2D STT integrator flag" begin
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

    p = LLGParams2D(
        5,
        5,
        0.0,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        0.0;
        u_stt = (0.2, 0.1),
        beta_stt = 0.15,
    )

    sol_off = evolve2D(copy(spins), (0.0, 0.2), p; reltol = 1e-7, abstol = 1e-7, stt = false)
    @test p.stt_active == false

    sol_on = evolve2D(copy(spins), (0.0, 0.2), p; reltol = 1e-7, abstol = 1e-7, stt = true)
    @test p.stt_active == false

    final_off = reshape(sol_off.u[end], size(spins))
    final_on = reshape(sol_on.u[end], size(spins))

    @test isapprox(final_off, spins; atol = 1e-10, rtol = 1e-10)
    @test maximum(abs.(final_on .- final_off)) > 1e-5
end
