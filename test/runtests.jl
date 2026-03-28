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
