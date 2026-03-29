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

    @test p.fields.Beff ≈ expected
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

    sol_off = evolve1D(copy(spins), (0.0, 0.2), p; reltol=1e-7, abstol=1e-7, stt=false)
    @test p.stt_active == false

    sol_on = evolve1D(copy(spins), (0.0, 0.2), p; reltol=1e-7, abstol=1e-7, stt=true)
    @test p.stt_active == false

    final_off = reshape(sol_off.u[end], size(spins))
    final_on = reshape(sol_on.u[end], size(spins))

    @test isapprox(final_off, spins; atol=1e-10, rtol=1e-10)
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

    @test p.fields.Beff ≈ expected
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

    @test p.fields.Beff ≈ expected
end
