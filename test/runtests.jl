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
