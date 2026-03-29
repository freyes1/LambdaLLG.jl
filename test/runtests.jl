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

@testset "1D spin-transfer torque" begin
    spins = [
        0.8  0.2  -0.4  -0.6;
        0.1  0.9   0.3  -0.2;
        0.6  0.4   0.85  0.77
    ]
    for i in 1:size(spins, 2)
        spins[:, i] ./= norm(spins[:, i])
    end

    p = LLGParams1D(
        size(spins, 2),
        0.0,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        0.0;
        u_stt = 0.45,
        beta_stt = 0.3,
    )

    fill!(p.fields.dS_1, 0.0)
    LambdaLLG.add_stt1D!(spins, p)

    expected = zeros(size(spins))

    for i in 1:size(spins, 2)
        grad = if i == 1
            spins[:, 2] - spins[:, 1]
        elseif i == size(spins, 2)
            spins[:, end] - spins[:, end - 1]
        else
            0.5 .* (spins[:, i + 1] - spins[:, i - 1])
        end

        expected[:, i] .= -p.u_stt .* grad .+ p.beta_stt * p.u_stt .* cross(spins[:, i], grad)
    end

    @test isapprox(p.fields.dS_1, expected)
end

@testset "1D STT integrator flag" begin
    spins = [
        0.0   0.0   0.0   0.0   0.0;
        0.2   0.6   1.0   0.6   0.2;
       -0.98 -0.8   0.0   0.8   0.98
    ]
    for i in 1:size(spins, 2)
        spins[:, i] ./= norm(spins[:, i])
    end

    p = LLGParams1D(
        size(spins, 2),
        0.0,
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        0.0;
        u_stt = 0.25,
        beta_stt = 0.2,
    )

    sol_off = evolve1D(copy(spins), (0.0, 0.2), p; reltol = 1e-7, abstol = 1e-7, stt = false)
    @test p.stt_active == false

    sol_on = evolve1D(copy(spins), (0.0, 0.2), p; reltol = 1e-7, abstol = 1e-7, stt = true)
    @test p.stt_active == false

    final_off = reshape(sol_off.u[end], size(spins))
    final_on = reshape(sol_on.u[end], size(spins))

    @test isapprox(final_off, spins; atol = 1e-10, rtol = 1e-10)
    @test maximum(abs.(final_on .- final_off)) > 1e-5
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

    for j in 1:size(spins, 3)
        for i in 1:size(spins, 2)
            left = i > 1 ? spins[:, i - 1, j] : zeros(3)
            right = i < size(spins, 2) ? spins[:, i + 1, j] : zeros(3)
            down = j > 1 ? spins[:, i, j - 1] : zeros(3)
            up = j < size(spins, 3) ? spins[:, i, j + 1] : zeros(3)
            expected[:, i, j] .= cross(Dx, right - left) + cross(Dy, up - down)
        end
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

    for j in 1:size(spins, 3)
        for i in 1:size(spins, 2)
            left = i > 1 ? spins[:, i - 1, j] : zeros(3)
            right = i < size(spins, 2) ? spins[:, i + 1, j] : zeros(3)
            down = j > 1 ? spins[:, i, j - 1] : zeros(3)
            up = j < size(spins, 3) ? spins[:, i, j + 1] : zeros(3)
            expected[:, i, j] .= cross(Dx, right - left) + cross(Dy, up - down)
        end
    end

    @test isapprox(p.fields.Beff, expected)
end

@testset "Initial-state constructors and normalization" begin
    direction = [1.0, 2.0, -2.0]
    direction ./= norm(direction)

    spins1D = uniform_state1D(5; direction = (1.0, 2.0, -2.0), magnitude = 2.5)
    @test size(spins1D) == (3, 5)
    for i in 1:size(spins1D, 2)
        @test isapprox(spins1D[:, i], 2.5 .* direction)
    end

    spins2D = uniform_state2D(4, 3; direction = (1.0, 2.0, -2.0))
    @test size(spins2D) == (3, 4, 3)
    for j in 1:size(spins2D, 3), i in 1:size(spins2D, 2)
        @test isapprox(spins2D[:, i, j], direction)
    end

    normalize_spins!(spins1D)
    normalize_spins!(spins2D)

    for i in 1:size(spins1D, 2)
        @test isapprox(norm(spins1D[:, i]), 1.0)
    end
    for j in 1:size(spins2D, 3), i in 1:size(spins2D, 2)
        @test isapprox(norm(spins2D[:, i, j]), 1.0)
    end
end

@testset "1D domain-wall painter" begin
    spins = uniform_state1D(41; direction = (0.0, 0.0, 1.0))

    paint_domain_wall!(
        spins;
        center = 21.0,
        width = 3.0,
        wall_type = :neel,
        reference_axis = (1.0, 0.0, 0.0),
    )
    normalize_spins!(spins)

    @test spins[3, 1] < -0.999
    @test abs(spins[1, 1]) < 0.01
    @test abs(spins[2, 1]) < 0.01
    @test spins[3, end] > 0.999
    @test abs(spins[1, end]) < 0.01
    @test abs(spins[2, end]) < 0.01
    @test isapprox(spins[:, 21], [1.0, 0.0, 0.0]; atol = 1e-6)

    for i in 1:size(spins, 2)
        @test isapprox(norm(spins[:, i]), 1.0; atol = 1e-12)
    end
end

@testset "2D painters" begin
    wall_state = uniform_state2D(31, 15; direction = (0.0, 0.0, 1.0))

    paint_domain_wall!(
        wall_state;
        center = (16.0, 8.0),
        width = 2.5,
        normal = (1.0, 0.0),
        wall_type = :bloch,
        reference_axis = (1.0, 0.0, 0.0),
    )
    normalize_spins!(wall_state)

    @test wall_state[3, 1, 8] < -0.999
    @test abs(wall_state[1, 1, 8]) < 0.01
    @test abs(wall_state[2, 1, 8]) < 0.01
    @test wall_state[3, end, 8] > 0.999
    @test abs(wall_state[1, end, 8]) < 0.01
    @test abs(wall_state[2, end, 8]) < 0.01
    @test isapprox(wall_state[:, 16, 8], [0.0, 1.0, 0.0]; atol = 1e-6)

    skyrmion_state = uniform_state2D(41, 41; direction = (0.0, 0.0, 1.0))
    paint_skyrmion!(
        skyrmion_state;
        center = (15.0, 21.0),
        radius = 6.0,
        width = 1.5,
        skyrmion_type = :neel,
        reference_axis = (1.0, 0.0, 0.0),
    )
    paint_skyrmion!(
        skyrmion_state;
        center = (29.0, 21.0),
        radius = 6.0,
        width = 1.5,
        skyrmion_type = :bloch,
        reference_axis = (1.0, 0.0, 0.0),
        chirality = -1.0,
    )
    normalize_spins!(skyrmion_state)

    @test skyrmion_state[3, 15, 21] < -0.999
    @test abs(skyrmion_state[1, 15, 21]) < 0.02
    @test abs(skyrmion_state[2, 15, 21]) < 0.02
    @test skyrmion_state[3, 29, 21] < -0.999
    @test abs(skyrmion_state[1, 29, 21]) < 0.02
    @test abs(skyrmion_state[2, 29, 21]) < 0.02
    @test skyrmion_state[3, 1, 1] > 0.99

    for j in 1:size(skyrmion_state, 3), i in 1:size(skyrmion_state, 2)
        @test isapprox(norm(skyrmion_state[:, i, j]), 1.0; atol = 1e-12)
    end
end
