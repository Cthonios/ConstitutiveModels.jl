function test_safe_log()
    for _ in 1:100
        input = 5. * rand()
        @test ConstitutiveModels.safe_log(input) ≈ log(input)
        @test isnan(ConstitutiveModels.safe_log(-input))
    end
end

test_safe_log()
