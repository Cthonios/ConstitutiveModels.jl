function test_state_variables()
    a = zero(SymmetricTensor{2, 3, Float64, 6})
    names = state_variable_names(a, "some_long_var_name")
    @test names[1] == "some_long_var_name_xx"
    @test names[2] == "some_long_var_name_xy"
    @test names[3] == "some_long_var_name_xz"
    @test names[4] == "some_long_var_name_yy"
    @test names[5] == "some_long_var_name_yz"
    @test names[6] == "some_long_var_name_zz"

    a = zero(Tensor{2, 3, Float64, 9})
    names = state_variable_names(a, "some_other_long_var_name")
    @test names[1] == "some_other_long_var_name_xx"
    @test names[2] == "some_other_long_var_name_yx"
    @test names[3] == "some_other_long_var_name_zx"
    @test names[4] == "some_other_long_var_name_xy"
    @test names[5] == "some_other_long_var_name_yy"
    @test names[6] == "some_other_long_var_name_zy"
    @test names[7] == "some_other_long_var_name_xz"
    @test names[8] == "some_other_long_var_name_yz"
    @test names[9] == "some_other_long_var_name_zz"
end

test_state_variables()
